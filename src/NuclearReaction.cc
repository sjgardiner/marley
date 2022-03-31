/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/MatrixElement.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/Reaction.hh"

using ME_Type = marley::MatrixElement::TransitionType;
using ProcType = marley::Reaction::ProcessType;
using CMode = marley::NuclearReaction::CoulombMode;

std::map< CMode, std::string > marley::NuclearReaction
  ::coulomb_mode_string_map_ =
{
  { CMode::NO_CORRECTION, "none" },
  { CMode::FERMI_FUNCTION, "Fermi" },
  { CMode::EMA, "EMA" },
  { CMode::MEMA, "MEMA" },
  { CMode::FERMI_AND_EMA, "Fermi-EMA" },
  { CMode::FERMI_AND_MEMA, "Fermi-MEMA" },
};

namespace {
  constexpr int BOGUS_TWO_J_VALUE = -99999;

  // Helper function that computes an approximate nuclear radius in natural
  // units (1/MeV)
  double nuclear_radius_natural_units(int A) {
    // First compute the approximate nuclear radius in fm
    double R = marley_utils::r0 * std::pow(A, marley_utils::ONE_THIRD);

    // Convert to natural units (MeV^(-1))
    double R_nat = R / marley_utils::hbar_c;

    return R_nat;
  }

}

marley::NuclearReaction::NuclearReaction(ProcType pt, int pdg_a, int pdg_b,
  int pdg_c, int pdg_d, int q_d,
  const std::shared_ptr<std::vector<marley::MatrixElement> >& mat_els,
  double lambda_dmcc)
  : q_d_( q_d ), matrix_elements_( mat_els ), lambda_dmcc_( lambda_dmcc )
{
  // Initialize the process type (NC, neutrino/antineutrino CC)
  process_type_ = pt;

  // Initialize the PDG codes for the 2->2 scatter particles
  pdg_a_ = pdg_a;
  pdg_b_ = pdg_b;
  pdg_c_ = pdg_c;
  pdg_d_ = pdg_d;

  // Get initial and final values of the nuclear charge and mass number from
  // the PDG codes
  Zi_ = (pdg_b_ % 10000000) / 10000;
  Ai_ = (pdg_b_ % 10000) / 10;
  Zf_ = (pdg_d_ % 10000000) / 10000;
  Af_ = (pdg_d_ % 10000) / 10;

  const marley::MassTable& mt = marley::MassTable::Instance();

  // Get the particle masses from the mass table
  ma_ = mt.get_particle_mass( pdg_a_ );
  mc_ = mt.get_particle_mass( pdg_c_ );

  // If the target (particle b) or residue (particle d)
  // has a particle ID greater than 10^9, assume that it
  // is an atom rather than a bare nucleus
  if ( pdg_b_ > 1000000000 ) mb_ = mt.get_atomic_mass( pdg_b_ );
  else mb_ = mt.get_particle_mass( pdg_b_ );

  if ( pdg_d_ > 1000000000 ) {
    // If particle d is an atom and is ionized as a result of this reaction
    // (e.g., q_d_ != 0), then approximate its ground-state ionized mass by
    // subtracting the appropriate number of electron masses from its atomic
    // (i.e., neutral) ground state mass.
    md_gs_ = mt.get_atomic_mass(pdg_d_)
      - (q_d_ * mt.get_particle_mass(marley_utils::ELECTRON));
  }
  else {
    md_gs_ = mt.get_particle_mass(pdg_d_);
  }

  KEa_threshold_ = (std::pow(mc_ + md_gs_, 2)
    - std::pow(ma_ + mb_, 2))/(2.*mb_);

  this->set_description();
}

// Fermi function used to calculate cross-sections
// The form used here is based on http://en.wikipedia.org/wiki/Beta_decay
// but rewritten for convenient use inside this class.
// Input: beta_c (3-momentum magnitude of particle c / total energy of particle
// c), where we assume that the ejectile (particle c) is the light product from
// 2-2 scattering.
double marley::NuclearReaction::fermi_function(double beta_c) const {

  // If the PDG code for particle c is positive, then it is a
  // negatively-charged lepton.
  bool c_minus = (pdg_c_ > 0);

  // Lorentz factor gamma for particle c
  double gamma_c = std::pow(1. - beta_c*beta_c, -marley_utils::ONE_HALF);

  double s = std::sqrt(1. - std::pow(marley_utils::alpha * Zf_, 2));

  // Estimate the nuclear radius (in MeV^(-1))
  double rho = nuclear_radius_natural_units( Af_ );

  // Sommerfeld parameter
  double eta = marley_utils::alpha * Zf_ / beta_c;

  // Adjust the value of eta if the light product from this reaction is
  // an antilepton
  if ( !c_minus ) eta *= -1;

  // Complex variable for the gamma function
  std::complex<double> a(s, eta);
  double b = std::tgamma(1 + 2*s);

  return 2 * (1 + s) * std::pow(2*beta_c*gamma_c*rho*mc_, 2*s-2)
    * std::exp(marley_utils::pi*eta) * std::norm(marley_utils::gamma(a))
    / std::pow(b, 2);
}

// Return the maximum residue excitation energy E_level that can
// be achieved in the lab frame for a given projectile kinetic energy KEa
// (this corresponds to the final particles all being produced
// at rest in the CM frame). This maximum level energy is used
// to find the allowed levels when creating events.
double marley::NuclearReaction::max_level_energy(double KEa) const {
  // Calculate the total CM frame energy using known quantities
  // from the lab frame
  double E_CM = std::sqrt(std::pow(ma_ + mb_, 2) + 2*mb_*KEa);
  // The maximum level energy is achieved when the final state
  // particles are produced at rest in the CM frame. Subtracting
  // the ground-state rest masses of particles c and d from the
  // total CM energy leaves us with the energy available to create
  // an excited level in the residue (particle d).
  return E_CM - mc_ - md_gs_;
}

double marley::NuclearReaction::threshold_kinetic_energy() const {
  return KEa_threshold_;
}

// Creates an event object by sampling the appropriate quantities and
// performing kinematic calculations
marley::Event marley::NuclearReaction::create_event(int pdg_a, double KEa,
  marley::Generator& gen) const
{
  // Check that the projectile supplied to this event is correct. If not, alert
  // the user that this event does not use the requested projectile.
  if ( pdg_a != pdg_a_ ) throw marley::Error(std::string("Could")
    + " not create this event. The requested projectile particle ID, "
    + std::to_string(pdg_a) + ", does not match the projectile"
    + " particle ID, " + std::to_string(pdg_a_) + ", in the reaction dataset.");

  // Sample a final residue energy level. First, check to make sure the given
  // projectile energy is above threshold for this reaction.
  if ( KEa < KEa_threshold_ ) throw std::range_error(std::string("Could")
    + " not create this event. Projectile kinetic energy " + std::to_string(KEa)
    + " MeV is below the threshold value " + std::to_string(KEa_threshold_)
    + " MeV.");

  /// @todo Add more error checks to NuclearReaction::create_event as necessary

  // Create an empty vector of sampling weights (partial total cross
  // sections to each kinematically accessible final level)
  std::vector<double> level_weights;

  // Create a discrete distribution object for level sampling.
  // Its default constructor creates a single weight of 1.
  // We will always explicitly give it weights to use when sampling
  // levels, so we won't worry about its default behavior.
  static std::discrete_distribution<size_t> ldist;

  // Compute the total cross section for a transition to each individual nuclear
  // level, and save the results in the level_weights vector (which will be
  // cleared by summed_xs_helper() before being loaded with the cross sections).
  // The summed_xs_helper() method can also be used for differential
  // (d\sigma/d\cos\theta_c^{CM}) cross sections, so supply a dummy cos_theta_c_cm
  // value and request total cross sections by setting the last argument to false.
  double dummy = 0.;
  double sum_of_xsecs = summed_xs_helper(pdg_a, KEa, dummy,
    &level_weights, false);

  // Note that the elements in matrix_elements_ are given in order of
  // increasing excitation energy (this is currently enforced by the reaction
  // data format and is checked during parsing). This ensures that we can
  // sample a matrix element index from level_weights (which is populated in
  // the same order by summed_xs_helper()) and have it refer to the correct
  // object.

  // Complain if none of the levels we have data for are kinematically allowed
  if ( level_weights.empty() ) {
    throw marley::Error("Could not create this event. The DecayScheme object"
      " associated with this reaction does not contain data for any"
      " kinematically accessible levels for a projectile kinetic energy of "
      + std::to_string(KEa) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) ) + " MeV).");
  }

  // Complain if the total cross section (the sum of all partial level cross
  // sections) is zero or negative (the latter is just to cover all possibilities).
  if ( sum_of_xsecs <= 0. ) {
    throw marley::Error("Could not create this event. All kinematically"
      " accessible levels for a projectile kinetic energy of "
      + std::to_string(KEa) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) )
      + " MeV) have vanishing matrix elements.");
  }

  // Create a list of parameters used to supply the weights to our discrete
  // level sampling distribution
  std::discrete_distribution<size_t>::param_type params( level_weights.begin(),
    level_weights.end() );

  // Sample a matrix_element using our discrete distribution and the
  // current set of weights
  size_t me_index = gen.sample_from_distribution( ldist, params );

  const auto& sampled_matrix_el = matrix_elements_->at( me_index );

  // Get the energy of the selected level.
  double E_level = sampled_matrix_el.level_energy();

  // Update the residue mass based on its excitation energy for the current
  // event
  md_ = md_gs_ + E_level;

  // Compute Mandelstam s, the ejectile's CM frame total energy, the magnitude
  // of the ejectile's CM frame 3-momentum, and the residue's CM frame total
  // energy.
  double s, Ec_cm, pc_cm, Ed_cm;
  two_two_scatter( KEa, s, Ec_cm, pc_cm, Ed_cm );

  // Determine the CM frame velocity of the ejectile
  double beta_c_cm = pc_cm / Ec_cm;

  // For a massive projectile, the angular distributions depend on the product
  // of the two speeds rather than just the speed of the ejectile. To avoid
  // needing to refactor the call to sample_cos_theta_c_cm() below, we can
  // cheat by doing the multiplication here. For massless projectiles, nothing
  // changes since they have beta = 1.

  // Determine the CM frame speed of the projectile. We actually need to do the
  // calculation only for massive projectiles.
  double beta_a_cm = 1.;
  if ( ma_ > 0. ) {
    // Projectile total energy
    double Ea_cm = ( s + ma_*ma_ - mb_*mb_ ) / ( 2 * std::sqrt(s) );

    // Projectile momentum
    double pa_cm = marley_utils::real_sqrt( std::pow(Ea_cm, 2) - ma_*ma_ );

    // Projectile speed
    beta_a_cm = pa_cm / Ea_cm;
  }

  // Sample a CM frame scattering cosine for the ejectile.
  double cos_theta_c_cm = sample_cos_theta_c_cm( sampled_matrix_el,
    beta_a_cm * beta_c_cm, gen );

  // Sample a CM frame azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the matrix elements are azimuthally invariant
  double phi_c_cm = gen.uniform_random_double( 0.,
    marley_utils::two_pi, false );

  // Load the initial residue twoJ and parity values into twoJ and P. These
  // variables are included in the event record and used by NucleusDecayer
  // to start the Hauser-Feshbach decay cascade.
  int twoJ = BOGUS_TWO_J_VALUE;
  marley::Parity P; // defaults to positive parity

  // Get access to the nuclear structure database owned by the Generator
  auto& sdb = gen.get_structure_db();

  // Retrieve the ground-state spin-parity of the initial nucleus
  int twoJ_gs;
  marley::Parity P_gs;

  sdb.get_gs_spin_parity( pdg_b_, twoJ_gs, P_gs );

  // For transitions to discrete nuclear levels, all we need to do is retrieve
  // these values directly from the Level object
  const marley::Level* final_lev = sampled_matrix_el.level();
  if ( final_lev ) {
    twoJ = final_lev->twoJ();
    P = final_lev->parity();
  }
  // For transitions to the continuum, we rely on the spin-parity selection
  // rules to determine suitable values of twoJ and P. In cases where more than
  // one value is allowed, assume equipartition, and sample a spin-parity based
  // on the relative nuclear level densities at the excitation energy of
  // interest.
  /// @todo Revisit the equipartition assumption made here in favor of
  /// something better motivated.
  else {

    // For a Fermi transition, the final spin-parity is always the same as the
    // initial one
    if ( sampled_matrix_el.type() == ME_Type::FERMI ) {
      twoJ = twoJ_gs;
      P = P_gs;
    }
    else if ( sampled_matrix_el.type() == ME_Type::GAMOW_TELLER ) {
      // For a Gamow-Teller transition, the final parity is the same as the
      // initial parity
      P = P_gs;

      // For a spin-zero initial state, take a shortcut: the final spin will
      // always be one.
      if ( twoJ_gs == 0 ) twoJ = 2;
      else {

        // For an initial state with a non-zero spin, make a vector storing
        // all of the spin values allowed by the GT selection rules.
        // Sample an allowed value assuming equipartition of spin. Use the
        // relative final nuclear level densities as sampling weights.

        std::vector<int> allowed_twoJs;
        std::vector<double> ld_weights;

        auto& ldm = sdb.get_level_density_model( pdg_d_ );

        for ( int myTwoJ = std::abs(twoJ_gs - 2); myTwoJ <= twoJ_gs + 2;
          myTwoJ += 2 )
        {
          allowed_twoJs.push_back( myTwoJ );
          ld_weights.push_back( ldm.level_density(E_level, myTwoJ, P) );
        }

        std::discrete_distribution<size_t> my_twoJ_dist( ld_weights.begin(),
          ld_weights.end() );

        size_t my_index = gen.sample_from_distribution( my_twoJ_dist );
        twoJ = allowed_twoJs.at( my_index );
      }
    }
    else throw marley::Error( "Unrecognized matrix element type encountered"
      " in marley::NuclearReaction::create_event()" );
  }

  MARLEY_LOG_DEBUG() << "Sampled a " << sampled_matrix_el.type_str()
    << " transition from the " << marley::TargetAtom( pdg_b_ )
    << " ground state (with spin-parity " << static_cast<double>( twoJ_gs ) / 2.
    << P_gs << ") to the " << marley::TargetAtom( pdg_d_ )
    << " level with Ex = " << E_level << " MeV and spin-parity "
    << static_cast<double>( twoJ ) / 2. << P;

  // Create the preliminary event object (after 2-->2 scattering, but before
  // de-excitation of the residual nucleus)
  marley::Event event = make_event_object( KEa, pc_cm, cos_theta_c_cm, phi_c_cm,
    Ec_cm, Ed_cm, E_level, twoJ, P );

  // Return the preliminary event object (to be processed later by the
  // NucleusDecayer class)
  return event;
}

// Compute the total reaction cross section (summed over all final nuclear levels)
// in units of MeV^(-2) using the center of momentum frame.
double marley::NuclearReaction::total_xs(int pdg_a, double KEa) const {
  double dummy_cos_theta = 0.;
  return summed_xs_helper(pdg_a, KEa, dummy_cos_theta, nullptr, false);
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM}
// summed over all final nuclear levels. This is done in units of MeV^(-2)
// using the center of momentum frame.
double marley::NuclearReaction::diff_xs(int pdg_a, double KEa,
  double cos_theta_c_cm) const
{
  return summed_xs_helper(pdg_a, KEa, cos_theta_c_cm, nullptr, true);
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM} for a
// transition to a particular final nuclear level. This is done in units of
// MeV^(-2) using the center of momentum frame.
double marley::NuclearReaction::diff_xs(const marley::MatrixElement& mat_el,
  double KEa, double cos_theta_c_cm) const
{
  // Check that the scattering cosine is within the physically meaningful range
  if ( std::abs(cos_theta_c_cm) > 1. ) return 0.;
  double beta_c_cm;
  double xsec = total_xs(mat_el, KEa, beta_c_cm, true);
  xsec *= mat_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm);
  return xsec;
}

// Helper function for total_xs and diff_xs()
double marley::NuclearReaction::summed_xs_helper(int pdg_a, double KEa,
  double cos_theta_c_cm, std::vector<double>* level_xsecs, bool differential)
  const
{
  // Check that the projectile supplied to this event is correct. If not,
  // return a total cross section of zero since this reaction is not available
  // for the given projectile.
  /// @todo Consider whether you should use an exception if pdg_a != pdg_a_
  if (pdg_a != pdg_a_) return 0.;

  // If we're evaluating a differential cross section, check that the
  // scattering cosine is within the physically meaningful range. If it's
  // not, then just return 0.
  if ( differential && std::abs(cos_theta_c_cm) > 1. ) return 0.;

  // If the projectile kinetic energy is zero (or negative), then
  // just return zero.
  if ( KEa <= 0. ) return 0.;

  // If we've been passed a vector to load with the partial cross sections
  // to each nuclear level, then clear it before storing them
  if ( level_xsecs ) level_xsecs->clear();

  double max_E_level = max_level_energy( KEa );
  double xsec = 0.;
  for ( const auto& mat_el : *matrix_elements_ ) {

    // Get the excitation energy for the current level
    double level_energy = mat_el.level_energy();

    // Exit the loop early if you reach a level with an energy that's too high
    if ( level_energy > max_E_level ) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // computing the total xs.
    if ( mat_el.strength() != 0. ) {

      // Set the check_max_E_level flag to false when calculating the total
      // cross section for this level (we've already verified that the
      // current level is kinematically accessible in the check against
      // max_E_level above)
      double beta_c_cm = 0.;
      double partial_xsec = total_xs(mat_el, KEa, beta_c_cm, false);

      // If a differential cross section (d\sigma / d\cos\theta_{CM})
      // is desired, then multiply by the appropriate angular factor
      if ( differential ) {
        partial_xsec *= mat_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm);
      }

      if ( std::isnan(partial_xsec) ) {
        MARLEY_LOG_WARNING() << "Partial cross section for reaction "
          << description_ << " gave NaN result.";
        MARLEY_LOG_DEBUG() << "Parameters were level energy = "
          << mat_el.level_energy() << " MeV, projectile kinetic energy = "
          << KEa << " MeV, and reduced matrix element = " << mat_el.strength();
        MARLEY_LOG_DEBUG() << "The partial cross section to this level"
          << " will be set to zero.";
        partial_xsec = 0.;
      }

      xsec += partial_xsec;

      // Store the partial cross section to the current individual nuclear
      // level if needed (i.e., if level_xsecs is not nullptr)
      if ( level_xsecs ) level_xsecs->push_back( partial_xsec );
    }
  }

  return xsec;
}

// Compute the total reaction cross section (in MeV^(-2)) for a transition to a
// particular nuclear level using the center of momentum frame
double marley::NuclearReaction::total_xs(const marley::MatrixElement& me,
  double KEa, double& beta_c_cm, bool check_max_E_level) const
{
  // Don't bother to compute anything if the matrix element vanishes for this
  // level
  if ( me.strength() == 0. ) return 0.;

  // Also don't proceed further if the reaction is below threshold (equivalently,
  // if the requested level excitation energy E_level exceeds that maximum
  // kinematically-allowed value). To avoid redundant checks of the threshold,
  // skip this check if check_max_E_level is set to false.
  if ( check_max_E_level ) {
    double max_E_level = max_level_energy( KEa );
    if ( me.level_energy() > max_E_level ) return 0.;
  }

  // The final nuclear mass (before nuclear de-excitations) is the sum of the
  // ground state residue mass plus the excitation energy of the accessed level
  double md2 = std::pow(md_gs_ + me.level_energy(), 2);

  // Compute Mandelstam s (the square of the total CM frame energy)
  double s = std::pow(ma_ + mb_, 2) + 2.*mb_*KEa;
  double sqrt_s = std::sqrt(s);

  // Compute CM frame total energies for two of the particles. Also
  // compute the magnitude of the ejectile CM frame momentum.
  double Eb_cm = (s + mb_*mb_ - ma_*ma_) / (2. * sqrt_s);
  double Ec_cm = (s + mc_*mc_ - md2) / (2. * sqrt_s);
  double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc_*mc_);

  // Compute the (dimensionless) speed of the ejectile in the CM frame
  beta_c_cm = pc_cm / Ec_cm;

  // CM frame total energy of the nuclear residue
  double Ed_cm = sqrt_s - Ec_cm;

  // Dot product of the four-momenta of particles c and d
  double pc_dot_pd = Ed_cm*Ec_cm + std::pow(pc_cm, 2);

  // Relative speed of particles c and d, computed with a manifestly
  // Lorentz-invariant expression
  double beta_rel_cd = marley_utils::real_sqrt(
    std::pow(pc_dot_pd, 2) - mc_*mc_*md2) / pc_dot_pd;

  // Common factors for the allowed approximation total cross sections
  // for both CC and NC reactions
  double total_xsec = ( Eb_cm * Ed_cm / s ) * Ec_cm
    * pc_cm * me.strength() / marley_utils::pi;

  // Apply extra factors based on the current process type
  if ( process_type_ == ProcessType::NeutrinoCC
    || process_type_ == ProcessType::AntiNeutrinoCC )
  {
    // Calculate a Coulomb correction factor using either a Fermi function
    // or the effective momentum approximation
    double factor_C = coulomb_correction_factor( beta_rel_cd );
    total_xsec *= marley_utils::GF2 * marley_utils::Vud2 * factor_C;
  }
  else if ( process_type_ == ProcessType::NC )
  {
    // For NC, extra factors are only needed for Fermi transitions (which
    // correspond to CEvNS since they can only access the nuclear ground state)
    if ( me.type() == ME_Type::FERMI ) {
      double Q_w = weak_nuclear_charge();
      total_xsec *= 0.25 * std::pow(Q_w, 2) * marley_utils::GF2;
    }
  }
  else if ( process_type_ == ProcessType::DMCC ) {

    // Compute CM frame projectile total energy, momentum, and speed. We need
    // the last of these when computing the cross section for a massive
    // projectile
    double Ea_cm = ( s + ma_*ma_ - mb_*mb_ ) / ( 2. * sqrt_s );
    double pa_cm = marley_utils::real_sqrt( std::pow(Ea_cm, 2) - ma_*ma_);
    double beta_a_cm = pa_cm / Ea_cm;

    double factor_C = coulomb_correction_factor( beta_rel_cd );

    total_xsec *= factor_C / ( 2. * beta_a_cm * std::pow(lambda_dmcc_, 4) );
  }
  else throw marley::Error("Unrecognized process type encountered in"
    " marley::NuclearReaction::total_xs()");

  return total_xsec;
}

// Sample an ejectile scattering cosine in the CM frame.
double marley::NuclearReaction::sample_cos_theta_c_cm(
  const marley::MatrixElement& matrix_el, double beta_c_cm,
  marley::Generator& gen) const
{
  // To avoid wasting time searching for the maximum of these distributions
  // for rejection sampling, set the maximum to the known value before
  // proceeding.
  double max;
  if ( matrix_el.type() == ME_Type::FERMI ) {
    // B(F)
    max = matrix_el.cos_theta_pdf(1., beta_c_cm);
  }
  else if (matrix_el.type() == ME_Type::GAMOW_TELLER) {
    // B(GT)
    max = matrix_el.cos_theta_pdf(-1., beta_c_cm);
  }
  else throw marley::Error("Unrecognized matrix element type "
    + std::to_string(matrix_el.type()) + " encountered while sampling a"
    " CM frame scattering angle");

  // Sample a CM frame scattering cosine using the appropriate distribution for
  // this matrix element.
  return gen.rejection_sample(
    [&matrix_el, &beta_c_cm](double cos_theta_c_cm)
    -> double { return matrix_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm); },
    -1., 1., max);
}

marley::Event marley::NuclearReaction::make_event_object(double KEa,
  double pc_cm, double cos_theta_c_cm, double phi_c_cm, double Ec_cm,
  double Ed_cm, double E_level, int twoJ, const marley::Parity& P) const
{
  marley::Event event = marley::Reaction::make_event_object(KEa, pc_cm,
    cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm, E_level, twoJ, P);

  // Assume that the target is a neutral atom (q_b = 0)
  event.target().set_charge(0);

  // Assign the correct charge to the residue
  event.residue().set_charge(q_d_);

  return event;
}

double marley::NuclearReaction::coulomb_correction_factor(double beta_rel_cd)
  const
{
  // Don't do anything if Coulomb corrections are switched off
  if ( coulomb_mode_ == CoulombMode::NO_CORRECTION ) return 1.;

  // Fermi function approach to the Coulomb correction
  double fermi_func = fermi_function( beta_rel_cd );

  // Unconditionally return the value of the Fermi function if the user
  // has configured things this way
  if ( coulomb_mode_ == CoulombMode::FERMI_FUNCTION ) return fermi_func;

  bool use_mema = false;
  if ( coulomb_mode_ == CoulombMode::MEMA
    || coulomb_mode_ == CoulombMode::FERMI_AND_MEMA )
  {
    use_mema = true;
  }

  // Effective momentum approximation for the Coulomb correction
  bool EMA_ok = false;
  double factor_EMA = ema_factor( beta_rel_cd, EMA_ok, use_mema );

  if ( coulomb_mode_ == CoulombMode::EMA || coulomb_mode_ == CoulombMode::MEMA ) {
    if ( EMA_ok ) return factor_EMA;
    else {
      std::string model_name( "EMA" );
      if ( coulomb_mode_ == CoulombMode::MEMA ) model_name = "MEMA";
      throw marley::Error( "Invalid " + model_name + " factor encountered"
        " in marley::NuclearReaction::coulomb_correction_factor()" );
    }
  }

  if ( coulomb_mode_ != CoulombMode::FERMI_AND_EMA
    && coulomb_mode_ != CoulombMode::FERMI_AND_MEMA )
  {
    throw marley::Error( "Unrecognized Coulomb correction mode encountered"
      " in marley::NuclearReaction::coulomb_correction_factor()" );
  }

  // If we've gotten this far, then we're interpolating between the Fermi
  // function and the (M)EMA correction factor. If the (modified) effective
  // momentum approximation is invalid because subtracting off the Coulomb
  // potential brings the reaction below threshold, then just use the Fermi
  // function
  if ( !EMA_ok ) return fermi_func;

  // Otherwise, choose the approach that yields the smaller correction (i.e., the
  // correction factor that is closest to unity).
  double diff_Fermi = std::abs( fermi_func - 1. );
  double diff_EMA = std::abs( factor_EMA - 1. );

  if ( diff_Fermi < diff_EMA ) return fermi_func;
  else return factor_EMA;
}

// Effective momentum approximation for the Coulomb correction factor
double marley::NuclearReaction::ema_factor(double beta_rel_cd, bool& ok,
  bool modified_ema) const
{
  // If particle c has a positive PDG code, then it is a negatively-charged
  // lepton
  bool minus_c = (pdg_c_ > 0);

  // Approximate nuclear radius (MeV^(-1))
  double R_nuc = nuclear_radius_natural_units( Af_ );

  // Approximate Coulomb potential
  double Vc = ( -3. * Zf_ * marley_utils::alpha ) / ( 2. * R_nuc );
  // Adjust if needed for a final-state charged antilepton
  if ( !minus_c ) Vc *= -1;

  // Like the Fermi function, this approximation uses a static nuclear Coulomb
  // potential (a sphere at the origin). Typically nuclear recoil is neglected,
  // allowing one to compute the effective lepton momentum in the lab frame. In
  // MARLEY's case, we do this by calculating it in the rest frame of the final
  // nucleus ("FNR" frame). We already have the relative speed of the final
  // nucleus and outgoing lepton, so this is easy.
  double gamma_rel_cd = std::pow( 1. - std::pow(beta_rel_cd, 2), -0.5 );

  // Check for numerical errors from the square root
  if ( !std::isfinite(gamma_rel_cd) ) {
    MARLEY_LOG_WARNING() << "Invalid beta_rel = " << beta_rel_cd
      << " encountered in marley::NuclearReaction::ema_factor()";
  }

  // Total energy of the outgoing lepton in the FNR frame
  double E_c_FNR = gamma_rel_cd * mc_;

  // Lepton momentum in FNR frame
  double p_c_FNR = beta_rel_cd * E_c_FNR;

  // Effective FNR frame total energy
  double E_c_FNR_eff = E_c_FNR - Vc;

  // If subtracting off the Coulomb potential drops the effective energy
  // below the lepton mass, then the expression for the effective momentum
  // will give an imaginary value. Signal this by setting the "ok" flag to
  // false.
  ok = ( E_c_FNR_eff >= mc_ );

  // Effective momentum in FNR frame
  double p_c_FNR_eff = marley_utils::real_sqrt(
    std::pow(E_c_FNR_eff, 2) - mc_*mc_ );

  // Coulomb correction factor for the original EMA
  double F_EMA = std::pow( p_c_FNR_eff / p_c_FNR, 2 );

  // Coulomb correction factor for the modified EMA
  double F_MEMA = ( p_c_FNR_eff * E_c_FNR_eff ) / ( p_c_FNR * E_c_FNR );

  if ( modified_ema ) return F_MEMA;
  else return F_EMA;
}

// Factor that appears in the cross section for coherent elastic
// neutrino-nucleus scattering (CEvNS), which corresponds to the Fermi
// component of NC scattering under the allowed approximation
double marley::NuclearReaction::weak_nuclear_charge() const
{
  int Ni = Ai_ - Zi_;
  double Qw = Ni - (1. - 4.*marley_utils::sin2thetaw)*Zi_;
  return Qw;
}

// Sets the description_ string based on the member PDG codes
void marley::NuclearReaction::set_description() {
  description_ = marley_utils::get_particle_symbol( pdg_a_ ) + " + ";
  description_ += std::to_string( Ai_ );
  description_ += marley_utils::element_symbols.at( Zi_ ) + " --> ";
  description_ += marley_utils::get_particle_symbol( pdg_c_ ) + " + ";
  description_ += std::to_string( Af_ );
  description_ += marley_utils::element_symbols.at( Zf_ );
  bool has_excited_state = false;
  for ( const auto& me : *matrix_elements_ ) {
    if ( me.level_energy() > 0. ) {
      has_excited_state = true;
      break;
    }
  }
  if ( has_excited_state ) description_ += '*';
  else description_ += " (g.s.)";
}

// Convert a string to a CoulombMode value
CMode marley::NuclearReaction::coulomb_mode_from_string(
  const std::string& str )
{
  for ( const auto& pair : coulomb_mode_string_map_ ) {
    if ( str == pair.second ) return pair.first;
  }
  throw marley::Error( "The string \"" + str + "\" was not recognized"
    " as a valid Coloumb mode setting" );
}

// Convert a CoulombMode value to a string
std::string marley::NuclearReaction::string_from_coulomb_mode( CMode mode ) {
  auto it = coulomb_mode_string_map_.find( mode );
  if ( it != coulomb_mode_string_map_.end() ) return it->second;
  else throw marley::Error( "Unrecognized CoulombMode value encountered in"
    " marley::NuclearReaction::string_from_coulomb_mode()" );
}

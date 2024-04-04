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
#include "marley/AllowedNuclearReaction.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/MatrixElement.hh"
#include "marley/NucleusDecayer.hh"

using ME_Type = marley::MatrixElement::TransitionType;
using ProcType = marley::Reaction::ProcessType;

namespace {
  constexpr int BOGUS_TWO_J_VALUE = -99999;
}

marley::AllowedNuclearReaction::AllowedNuclearReaction( ProcType pt,
  int pdg_a, int pdg_b, int pdg_c, int pdg_d, int q_d,
  const std::shared_ptr< std::vector<marley::MatrixElement> >& mat_els,
  marley::CoulombCorrector::CoulombMode mode )
  : marley::NuclearReaction( pt, pdg_a, pdg_b, pdg_c, pdg_d, q_d ),
  matrix_elements_( mat_els ), coulomb_corrector_( pdg_c, pdg_d, mode )
{
}

// Creates an event object by sampling the appropriate quantities and
// performing kinematic calculations
std::shared_ptr< HepMC3::GenEvent > marley::AllowedNuclearReaction
  ::create_event( int pdg_a, double KEa, marley::Generator& gen ) const
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

  /// @todo Add more error checks to AllowedNuclearReaction::create_event() as
  /// necessary

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
  // (d\sigma/d\cos\theta_c^{CM}) cross sections, so supply a dummy
  // cos_theta_c_cm value and request total cross sections by setting the last
  // argument to false.
  double dummy = 0.;
  double sum_of_xsecs = summed_xs_helper( pdg_a, KEa, dummy,
    &level_weights, false );

  // Note that the elements in matrix_elements_ are given in order of
  // increasing excitation energy (this is currently enforced by the reaction
  // data format and is checked during parsing). This ensures that we can
  // sample a matrix element index from level_weights (which is populated in
  // the same order by summed_xs_helper()) and have it refer to the correct
  // object.

  // Complain if none of the levels we have data for are kinematically allowed
  if ( level_weights.empty() ) {
    throw marley::Error( "Could not create this event. The DecayScheme object"
      " associated with this reaction does not contain data for any"
      " kinematically accessible levels for a projectile kinetic energy of "
      + std::to_string( KEa ) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) ) + " MeV)." );
  }

  // Complain if the total cross section (the sum of all partial level cross
  // sections) is zero or negative (the latter is just to cover all possibilities).
  if ( sum_of_xsecs <= 0. ) {
    throw marley::Error( "Could not create this event. All kinematically"
      " accessible levels for a projectile kinetic energy of "
      + std::to_string( KEa ) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) )
      + " MeV) have vanishing matrix elements." );
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

  // Sample a CM frame scattering cosine for the ejectile.
  double cos_theta_c_cm = sample_cos_theta_c_cm( sampled_matrix_el,
    beta_c_cm, gen );

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
      " in marley::AllowedNuclearReaction::create_event()" );
  }

  MARLEY_LOG_DEBUG() << "Sampled a " << sampled_matrix_el.type_str()
    << " transition from the " << marley::TargetAtom( pdg_b_ )
    << " ground state (with spin-parity " << static_cast<double>( twoJ_gs ) / 2.
    << P_gs << ") to the " << marley::TargetAtom( pdg_d_ )
    << " level with Ex = " << E_level << " MeV and spin-parity "
    << static_cast<double>( twoJ ) / 2. << P;

  // Create the preliminary event object (after 2-->2 scattering, but before
  // de-excitation of the residual nucleus)
  auto event = this->make_event_object( KEa, pc_cm, cos_theta_c_cm, phi_c_cm,
    Ec_cm, Ed_cm, E_level, twoJ, P );

  // Return the preliminary event object (to be processed later by the
  // NucleusDecayer class)
  return event;
}

// Compute the total reaction cross section (summed over all final nuclear
// levels) in units of MeV^(-2) using the center of momentum frame.
double marley::AllowedNuclearReaction::total_xs( int pdg_a, double KEa ) const
{
  double dummy_cos_theta = 0.;
  return summed_xs_helper( pdg_a, KEa, dummy_cos_theta, nullptr, false );
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM}
// summed over all final nuclear levels. This is done in units of MeV^(-2)
// using the center of momentum frame.
double marley::AllowedNuclearReaction::diff_xs( int pdg_a, double KEa,
  double cos_theta_c_cm ) const
{
  return summed_xs_helper( pdg_a, KEa, cos_theta_c_cm, nullptr, true );
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM} for a
// transition to a particular final nuclear level. This is done in units of
// MeV^(-2) using the center of momentum frame.
double marley::AllowedNuclearReaction::diff_xs(
  const marley::MatrixElement& mat_el, double KEa, double cos_theta_c_cm ) const
{
  // Check that the scattering cosine is within the physically meaningful range
  if ( std::abs(cos_theta_c_cm) > 1. ) return 0.;
  double beta_c_cm;
  double xsec = total_xs(mat_el, KEa, beta_c_cm, true);
  xsec *= mat_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm);
  return xsec;
}

// Helper function for total_xs and diff_xs()
double marley::AllowedNuclearReaction::summed_xs_helper( int pdg_a, double KEa,
  double cos_theta_c_cm, std::vector<double>* level_xsecs, bool differential )
  const
{
  // Check that the projectile supplied to this event is correct. If not,
  // return a total cross section of zero since this reaction is not available
  // for the given projectile.
  /// @todo Consider whether you should use an exception if pdg_a != pdg_a_
  if ( pdg_a != pdg_a_ ) return 0.;

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
double marley::AllowedNuclearReaction::total_xs(
  const marley::MatrixElement& me, double KEa, double& beta_c_cm,
  bool check_max_E_level ) const
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
  double md2 = std::pow( md_gs_ + me.level_energy(), 2 );

  // Compute Mandelstam s (the square of the total CM frame energy)
  double s = std::pow( ma_ + mb_, 2 ) + 2.*mb_*KEa;
  double sqrt_s = std::sqrt( s );

  // Compute CM frame total energies for two of the particles. Also
  // compute the magnitude of the ejectile CM frame momentum.
  double Eb_cm = ( s + mb_*mb_ - ma_*ma_ ) / ( 2. * sqrt_s );
  double Ec_cm = ( s + mc_*mc_ - md2 ) / ( 2. * sqrt_s );
  double pc_cm = marley_utils::real_sqrt( std::pow(Ec_cm, 2) - mc_*mc_ );

  // Compute the (dimensionless) speed of the ejectile in the CM frame
  beta_c_cm = pc_cm / Ec_cm;

  // CM frame total energy of the nuclear residue
  double Ed_cm = sqrt_s - Ec_cm;

  // Dot product of the four-momenta of particles c and d
  double pc_dot_pd = Ed_cm*Ec_cm + std::pow( pc_cm, 2 );

  // Relative speed of particles c and d, computed with a manifestly
  // Lorentz-invariant expression
  double beta_rel_cd = marley_utils::real_sqrt(
    std::pow(pc_dot_pd, 2) - mc_*mc_*md2 ) / pc_dot_pd;

  // Common factors for the allowed approximation total cross sections
  // for both CC and NC reactions
  double total_xsec = ( marley_utils::GF2 / marley_utils::pi )
    * ( Eb_cm * Ed_cm / s ) * Ec_cm * pc_cm * me.strength();

  // Apply extra factors based on the current process type
  if ( process_type_ == ProcessType::NeutrinoCC_Discrete
    || process_type_ == ProcessType::AntiNeutrinoCC_Discrete )
  {
    // Calculate a Coulomb correction factor using either a Fermi function
    // or the effective momentum approximation
    double factor_C = coulomb_corrector_.coulomb_correction_factor(
      beta_rel_cd );
    total_xsec *= marley_utils::Vud2 * factor_C;
  }
  else if ( process_type_ == ProcessType::NC_Discrete )
  {
    // For NC, extra factors are only needed for Fermi transitions (which
    // correspond to CEvNS since they can only access the nuclear ground state)
    if ( me.type() == ME_Type::FERMI ) {
      double Q_w = weak_nuclear_charge();
      total_xsec *= 0.25*std::pow( Q_w, 2 );
    }
  }
  else throw marley::Error( "Unrecognized process type encountered in"
    " marley::AllowedNuclearReaction::total_xs()" );

  return total_xsec;
}

// Sample an ejectile scattering cosine in the CM frame.
double marley::AllowedNuclearReaction::sample_cos_theta_c_cm(
  const marley::MatrixElement& matrix_el, double beta_c_cm,
  marley::Generator& gen ) const
{
  // To avoid wasting time searching for the maximum of these distributions
  // for rejection sampling, set the maximum to the known value before
  // proceeding.
  double max;
  if ( matrix_el.type() == ME_Type::FERMI ) {
    // B(F)
    max = matrix_el.cos_theta_pdf( 1., beta_c_cm );
  }
  else if ( matrix_el.type() == ME_Type::GAMOW_TELLER ) {
    // B(GT)
    max = matrix_el.cos_theta_pdf( -1., beta_c_cm );
  }
  else throw marley::Error("Unrecognized matrix element type "
    + std::to_string(matrix_el.type()) + " encountered while sampling a"
    " CM frame scattering angle");

  // Sample a CM frame scattering cosine using the appropriate distribution for
  // this matrix element.
  return gen.rejection_sample(
    [ &matrix_el, &beta_c_cm ]( double cos_theta_c_cm )
    -> double { return matrix_el.cos_theta_pdf( cos_theta_c_cm, beta_c_cm ); },
    -1., 1., max );
}

std::shared_ptr< HepMC3::GenEvent > marley::AllowedNuclearReaction
  ::make_event_object( double KEa, double pc_cm, double cos_theta_c_cm,
  double phi_c_cm, double Ec_cm, double Ed_cm, double E_level, int twoJ,
  const marley::Parity& P ) const
{
  auto event = marley::Reaction::make_event_object( KEa, pc_cm,
    cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm, E_level, twoJ, P );

  this->set_charge_attributes( event );

  return event;
}

// Adds an indication of whether the reaction populates excited levels of the
// daughter nucleus or only accesses the ground state (e.g., CEvNS)
void marley::AllowedNuclearReaction::set_description() {
  // Initialize the description_ member variable with the basic information
  // first
  marley::NuclearReaction::set_description();
  // Now decide whether this reaction can access excited levels in the daughter
  // nucleus
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

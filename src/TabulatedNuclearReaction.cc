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

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/MassTable.hh"
#include "marley/TabulatedNuclearReaction.hh"
#include "marley/marley_utils.hh"

marley::TabulatedNuclearReaction::TabulatedNuclearReaction(
  Reaction::ProcessType pt, int pdg_a, int pdg_b, int pdg_c, int pdg_d,
  int q_d, const std::shared_ptr<TabulatedXSec>& txsec )
  : marley::NuclearReaction( pt, pdg_a, pdg_b, pdg_c, pdg_d, q_d ),
  xsec_( txsec )
{
}

double marley::TabulatedNuclearReaction::total_xs( int pdg_a,
  double KEa ) const
{
  if ( pdg_a != pdg_a_ ) return 0.;
  return xsec_->integral( pdg_a, KEa );
}

marley::Event marley::TabulatedNuclearReaction::create_event( int pdg_a,
  double KEa, marley::Generator& gen ) const
{
  // TODO: reduce code duplication here with AllowedNuclearReaction using the
  // common base class NuclearReaction

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

  // Select a specific multipole to use for the current event using the
  // individual total cross sections
  std::vector<double> multipole_weights;
  std::vector<marley::TabulatedXSec::MultipoleLabel> multipoles;
  std::vector<double> diff_max_values;
  const auto& table_map = xsec_->get_table_map();
  double sum_of_xsecs = 0.;
  for ( const auto& pair : table_map ) {
    const auto& ml = pair.first;
    double diff_max;
    double total_xsec = xsec_->integral( pdg_a_, KEa, ml, diff_max );

    sum_of_xsecs += total_xsec;

    multipole_weights.push_back( total_xsec );
    multipoles.push_back( ml );
    diff_max_values.push_back( diff_max );
  }

  // If there are no multipole weights, we can't go on. Complain if this
  // is the case.
  if ( multipole_weights.empty() ) {
    throw marley::Error("Could not create this event. The TabulatedXSec object"
      " associated with this reaction does not own any nuclear response"
      " tables." );
  }

  // Complain if the total cross section (the sum of all partial cross
  // sections) is zero or negative (the latter is just to cover all
  // possibilities).
  if ( sum_of_xsecs <= 0. ) {
    throw marley::Error( "Could not create this event. All"
      " multipole total cross sections are nonpositive." );
  }

  // Create a discrete distribution based on the weights. This will be
  // used to choose a single multipole for the current event.
  std::discrete_distribution<size_t> multipole_dist( multipole_weights.begin(),
    multipole_weights.end() );

  // Sample a matrix_element using our discrete distribution and the
  // current set of weights
  size_t multipole_index = gen.sample_from_distribution( multipole_dist );

  // Label of the multipole chosen for this event
  const auto& sampled_ml = multipoles.at( multipole_index );
  // Maximum value of the differential cross section for this multipole.
  // This will be used for rejection sampling of inclusive kinematics below.
  double diff_max = diff_max_values.at( multipole_index );

  // ResponseTable object to use for computing the differential cross section
  // during kinematic sampling below
  const auto& rt = xsec_->get_table( sampled_ml );

  // Get the minimum energy transfer value that appears in the table of
  // nuclear responses
  double table_wmin = rt.w_min();
  double table_wmax = rt.w_max();

  // Choose a reasonable sampling interval for the energy transfer
  double Ea = KEa + ma_; // Projectile total energy

  // Set the lower bound for the energy transfer to be either zero or
  // the lowest tabulated value (whichever is larger)
  double wmin = std::max( 0., table_wmin );
  // Set the upper bound for the energy transfer to be either the projectile
  // energy minus the final lepton mass or the highest tabulated value
  // (whichever is lower)
  double wmax = std::min( Ea - mc_, table_wmax );

  // Sample values for the energy transfer and scattering cosine using the
  // differential cross section for the chosen multipole.
  // Use a simple rejection sampling technique.
  double w, ctl, diff, y;
  do {
    w = gen.uniform_random_double( wmin, wmax, true );
    ctl = gen.uniform_random_double( -1., 1., true );
    diff = xsec_->diff_xsec( pdg_a_, KEa, w, ctl, sampled_ml );
    y = gen.uniform_random_double( 0., diff_max, true );
  } while ( y > diff );

  // Sample a lab-frame azimuthal scattering angle uniformly
  double phi_c = gen.uniform_random_double( 0., marley_utils::two_pi, false );

  // Load the initial residue twoJ and parity values into twoJ and P. These
  // variables are included in the event record and used by NucleusDecayer to
  // start the Hauser-Feshbach decay cascade.
  // NOTE: right now, these are taken directly from the multipole involved in
  // the current event. This is only valid for scattering on a 0+ target
  // nucleus
  // TODO: revisit this assumption and do something better
  int twoJ = 2 * sampled_ml.J_;
  marley::Parity P = sampled_ml.Pi_; // defaults to positive parity

  // Sine of the ejectile scattering angle
  double stl = marley_utils::real_sqrt( 1. - std::pow(ctl, 2) );

  // Get the ejectile total energy and momentum using the sampled energy
  // transfer
  double Ec = Ea - w;
  double pc = marley_utils::real_sqrt( Ec*Ec - mc_*mc_ );

  // Determine the Cartesian components of the ejectile's lab-frame momentum
  double pc_x = stl * std::cos( phi_c ) * pc;
  double pc_y = stl * std::sin( phi_c ) * pc;
  double pc_z = ctl * pc;

  // Determine the magnitude of the lab-frame 3-momentum of the projectile
  double pa = marley_utils::real_sqrt( KEa*(KEa + 2*ma_) );

  // Create particle objects representing the projectile, target, and ejectile
  // in the lab frame
  marley::Particle projectile( pdg_a_, Ea, 0., 0., pa, ma_ );
  marley::Particle target( pdg_b_, mb_, 0., 0., 0., mb_ );
  marley::Particle ejectile( pdg_c_, Ec, pc_x, pc_y, pc_z, mc_ );

  // Get the 4-momentum of the residue in the lab frame using conservation
  double Ed = projectile.total_energy() + target.total_energy()
    - ejectile.total_energy();
  double pd_x = projectile.px() + target.px() - ejectile.px();
  double pd_y = projectile.py() + target.py() - ejectile.py();
  double pd_z = projectile.pz() + target.pz() - ejectile.pz();

  // Determine the residue mass from its 4-momentum
  md_ = marley_utils::real_sqrt( Ed*Ed - pd_x*pd_x - pd_y*pd_y - pd_z*pd_z );

  // The excitation energy is the mass difference between this mass and
  // the residue's ground-state mass
  double Ex = md_ - md_gs_;

  // TODO: add some sanity-checking here for the excitation energy

  // Create a particle to represent the residue
  marley::Particle residue( pdg_d_, Ed, pd_x, pd_y, pd_z, md_, q_d_ );

  // Create the event object and load it with the appropriate information
  marley::Event event( projectile, target, ejectile, residue, Ex, twoJ, P );
  return event;
}

// Adds an indication to the description that the daughter nucleus will always
// be left in an excited state.
// TODO: revisit this as needed. I assume here that the
// TabulatedNuclearReaction class will always be used for calculations in the
// unbound continuum of nuclear levels
void marley::TabulatedNuclearReaction::set_description() {
  marley::NuclearReaction::set_description();
  description_ += '*';
}

/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
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

#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/MatrixElement.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/Reaction.hh"

using ProcType = marley::Reaction::ProcessType;

marley::NuclearReaction::NuclearReaction( ProcType pt, int pdg_a, int pdg_b,
  int pdg_c, int pdg_d, int q_d )
  : q_d_( q_d ), coulomb_corrector_( pdg_c, pdg_d )
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
    md_gs_ = mt.get_atomic_mass( pdg_d_ )
      - ( q_d_ * mt.get_particle_mass(marley_utils::ELECTRON) );
  }
  else {
    md_gs_ = mt.get_particle_mass( pdg_d_ );
  }

  KEa_threshold_ = ( std::pow(mc_ + md_gs_, 2)
    - std::pow(ma_ + mb_, 2) ) / ( 2.*mb_ );

  this->set_description();
}

// Return the maximum residue excitation energy E_level that can
// be achieved in the lab frame for a given projectile kinetic energy KEa
// (this corresponds to the final particles all being produced
// at rest in the CM frame). This maximum level energy is used
// to find the allowed levels when creating events.
double marley::NuclearReaction::max_level_energy( double KEa ) const {
  // Calculate the total CM frame energy using known quantities
  // from the lab frame
  double E_CM = std::sqrt( std::pow(ma_ + mb_, 2) + 2*mb_*KEa );
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

// Factor that appears in the cross section for coherent elastic
// neutrino-nucleus scattering (CEvNS), which corresponds to the Fermi
// component of NC scattering under the allowed approximation
double marley::NuclearReaction::weak_nuclear_charge() const
{
  int Ni = Ai_ - Zi_;
  double Qw = Ni - ( 1. - 4.*marley_utils::sin2thetaw )*Zi_;
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
}

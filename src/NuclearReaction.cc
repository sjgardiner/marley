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
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

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
  // Helper function that computes an approximate nuclear radius in natural
  // units (1/MeV)
  double nuclear_radius_natural_units( int A ) {
    // First compute the approximate nuclear radius in fm
    double R = marley_utils::r0 * std::pow( A, marley_utils::ONE_THIRD );

    // Convert to natural units (MeV^(-1))
    double R_nat = R / marley_utils::hbar_c;

    return R_nat;
  }

}

marley::NuclearReaction::NuclearReaction( ProcType pt, int pdg_a, int pdg_b,
  int pdg_c, int pdg_d, int q_d )
  : q_d_( q_d )
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

// Fermi function used to calculate cross-sections
// The form used here is based on http://en.wikipedia.org/wiki/Beta_decay
// but rewritten for convenient use inside this class.
// Input: beta_c (3-momentum magnitude of particle c / total energy of particle
// c), where we assume that the ejectile (particle c) is the light product from
// 2-2 scattering.
double marley::NuclearReaction::fermi_function( double beta_c ) const {

  // If the PDG code for particle c is positive, then it is a
  // negatively-charged lepton.
  bool c_minus = ( pdg_c_ > 0 );

  // Lorentz factor gamma for particle c
  double gamma_c = std::pow( 1. - beta_c*beta_c, -marley_utils::ONE_HALF );

  double s = std::sqrt( 1. - std::pow(marley_utils::alpha * Zf_, 2) );

  // Estimate the nuclear radius (in MeV^(-1))
  double rho = nuclear_radius_natural_units( Af_ );

  // Sommerfeld parameter
  double eta = marley_utils::alpha * Zf_ / beta_c;

  // Adjust the value of eta if the light product from this reaction is
  // an antilepton
  if ( !c_minus ) eta *= -1;

  // Complex variable for the gamma function
  std::complex<double> a( s, eta );
  double b = std::tgamma( 1 + 2*s );

  return 2. * ( 1. + s ) * std::pow( 2.*beta_c*gamma_c*rho*mc_, 2.*s - 2. )
    * std::exp( marley_utils::pi*eta ) * std::norm( marley_utils::gamma(a) )
    / std::pow( b, 2 );
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

double marley::NuclearReaction::coulomb_correction_factor( double beta_rel_cd )
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

  if ( coulomb_mode_ == CoulombMode::EMA
    || coulomb_mode_ == CoulombMode::MEMA )
  {
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

  // Otherwise, choose the approach that yields the smaller correction (i.e.,
  // the correction factor that is closest to unity).
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
  bool minus_c = ( pdg_c_ > 0 );

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

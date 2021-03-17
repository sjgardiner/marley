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

#include "marley/marley_utils.hh"
#include "marley/CoulombCorrector.hh"
#include "marley/Error.hh"
#include "marley/Logger.hh"
#include "marley/MassTable.hh"

using CMode = marley::CoulombCorrector::CoulombMode;

std::map< CMode, std::string > marley::CoulombCorrector
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

  bool is_charged_lepton_or_antilepton( int pdg ) {
    int abs_pdg = std::abs( pdg );
    if ( abs_pdg != marley_utils::ELECTRON && abs_pdg != marley_utils::MUON
      && abs_pdg != marley_utils::TAU ) return false;
    return true;
  }

}

marley::CoulombCorrector::CoulombCorrector( int pdg_c,
  int pdg_d, CoulombMode mode ) : coulomb_mode_( mode ), pdg_c_( pdg_c )
{
  // Extract the mass number and proton number from the residue PDG code
  Af_ = marley_utils::get_particle_A( pdg_d );
  Zf_ = marley_utils::get_particle_Z( pdg_d );

  // Get the ejectile mass from the mass table
  const marley::MassTable& mt = marley::MassTable::Instance();
  mc_ = mt.get_particle_mass( pdg_c_ );

  // Do a sanity check that the ejectile is actually a charged (anti)lepton. If
  // it's not, then set the Coulomb mode to NO_CORRECTION to prevent misuse
  // of this object.
  if ( is_charged_lepton_or_antilepton(pdg_c_) ) {
    coulomb_mode_ = CoulombMode::NO_CORRECTION;
  }
}

// Fermi function used to calculate cross-sections
// The form used here is based on http://en.wikipedia.org/wiki/Beta_decay
// but rewritten for convenient use inside this class.
// Input: beta_c (3-momentum magnitude of particle c / total energy of particle
// c), where we assume that the ejectile (particle c) is the light product from
// 2-2 scattering.
double marley::CoulombCorrector::fermi_function( double beta_c ) const {

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

double marley::CoulombCorrector::coulomb_correction_factor( double beta_rel_cd )
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
        " in marley::CoulombCorrector::coulomb_correction_factor()" );
    }
  }

  if ( coulomb_mode_ != CoulombMode::FERMI_AND_EMA
    && coulomb_mode_ != CoulombMode::FERMI_AND_MEMA )
  {
    throw marley::Error( "Unrecognized Coulomb correction mode encountered"
      " in marley::CoulombCorrector::coulomb_correction_factor()" );
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
double marley::CoulombCorrector::ema_factor(double beta_rel_cd, bool& ok,
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
      << " encountered in marley::CoulombCorrector::ema_factor()";
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

// Convert a string to a CoulombMode value
CMode marley::CoulombCorrector::coulomb_mode_from_string(
  const std::string& str )
{
  for ( const auto& pair : coulomb_mode_string_map_ ) {
    if ( str == pair.second ) return pair.first;
  }
  throw marley::Error( "The string \"" + str + "\" was not recognized"
    " as a valid Coloumb mode setting" );
}

// Convert a CoulombMode value to a string
std::string marley::CoulombCorrector::string_from_coulomb_mode( CMode mode ) {
  auto it = coulomb_mode_string_map_.find( mode );
  if ( it != coulomb_mode_string_map_.end() ) return it->second;
  else throw marley::Error( "Unrecognized CoulombMode value encountered in"
    " marley::CoulombCorrector::string_from_coulomb_mode()" );
}

// Sets the method to use for computing Coulomb corrections. If the PDG
// code for particle c is not a charged (anti)lepton, then ignore the
// input and set the mode to "no correction"
void marley::CoulombCorrector::set_coulomb_mode( CoulombMode mode ) {
  if ( is_charged_lepton_or_antilepton(pdg_c_) ) coulomb_mode_ = mode;
  else coulomb_mode_ = CoulombMode::NO_CORRECTION;
}

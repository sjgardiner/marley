/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
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

// Standard library includes

// MARLEY includes
#include "marley/FormFactor.hh"
#include "marley/marley_utils.hh"
#include "marley/Error.hh"

// The default constructor sets the Q^2 scaling mode
marley::FormFactor::FormFactor( FFScalingMode mode ) : ff_scaling_mode_(mode) {}

// Helper map for converting strings to Q^2 scaling modes
std::map < marley::FormFactor::FFScalingMode, std::string > marley::FormFactor::ff_scaling_mode_map_ = {
  { marley::FormFactor::FFScalingMode::DIPOLE, "dipole" },
  { marley::FormFactor::FFScalingMode::FLAT, "flat" }
};

// Convert a string to a Q^2 scaling mode
marley::FormFactor::FFScalingMode marley::FormFactor::ff_scaling_mode_from_string(
  const std::string& str ) 
{
  for (const auto& pair : ff_scaling_mode_map_) {
    if (pair.second == str) return pair.first;
  }
  throw marley::Error("Invalid Q^2 scaling mode string: " + str);
}

// Convert a Q^2 scaling mode to a string
std::string marley::FormFactor::string_from_ff_scaling_mode( FFScalingMode mode )
{
  auto it = ff_scaling_mode_map_.find(mode);
  if (it != ff_scaling_mode_map_.end()) return it->second;
  else throw marley::Error("Invalid FFScalingMode value encountered in"
    " marley::FormFactor::string_from_q2_scaling_mode" );
}

// Dipole dependence form factor
double marley::FormFactor::dipole( double Q2, double M ) const {
  return 1.0 / ( ( 1.0 + Q2 / M / M ) * ( 1.0 + Q2 / M / M ) );
}

// F1_n form factor
double marley::FormFactor::F1_n( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) {
    return 1.0;
  } else if (ff_scaling_mode_ == FFScalingMode::DIPOLE) {

    double prefactor = ( marley_utils::mu_n * Q2 )
                    / (4 * marley_utils::m_nucleon2 + Q2);

    return prefactor * dipole(Q2, M);

  } else {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }
}

// F1_p form factor
double marley::FormFactor::F1_p( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) {
    return 1.0;
  } else if (ff_scaling_mode_ == FFScalingMode::DIPOLE) {

    double prefactor = ( 4 * marley_utils::m_nucleon2 + marley_utils::mu_p * Q2 )
                    / (4 * marley_utils::m_nucleon2 + Q2);

    return prefactor * dipole(Q2, M);

  } else {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }
}

// F2_n form factor
double marley::FormFactor::F2_n( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) {
    return 1.0;
  } else if (ff_scaling_mode_ == FFScalingMode::DIPOLE) {

    double prefactor = ( 4 * marley_utils::m_nucleon2 * marley_utils::mu_n )
                      / (4 * marley_utils::m_nucleon2 + Q2);

    return prefactor * dipole(Q2, M);

  } else {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }
}

// F2_p form factor
double marley::FormFactor::F2_p( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) {
    return 1.0;
  } else if (ff_scaling_mode_ == FFScalingMode::DIPOLE) {

    double prefactor = ( 4 * marley_utils::m_nucleon2 * (marley_utils::mu_p - 1) )
                    / (4 * marley_utils::m_nucleon2 + Q2);

    return prefactor * dipole(Q2, M);

  } else {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }
}

// F3 form factor
double marley::FormFactor::F3( double Q2, double M ) const { return 0.0; }

// GA form factor
double marley::FormFactor::GA( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) {
    return 1.0;
  } else if (ff_scaling_mode_ == FFScalingMode::DIPOLE) {

    return marley_utils::g_A * dipole(Q2, M);

  } else {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }
  
}

// GP form factor
double marley::FormFactor::GP( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) {
    return 1.0;
  } else if (ff_scaling_mode_ == FFScalingMode::DIPOLE) {

    double prefactor = (2 * marley_utils::m_nucleon2 * marley_utils::g_A) 
              / (marley_utils::m_pion * marley_utils::m_pion + Q2);
    
    return prefactor * dipole(Q2, M);

  } else {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }
}

// G3 form factor
double marley::FormFactor::G3( double Q2, double M ) const { return 0.0; }

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

#pragma once

// Standard library includes
#include <functional>
#include <map>
#include <string>

// MARLEY includes


namespace marley {

  /// @brief Class implementing the different nuclear form factors used in
  /// MARLEY

  class FormFactor {

    public:

      /// @brief Enumeration of the possible Q^2 scaling modes for the form factors
      enum class FFScalingMode { DIPOLE, FLAT };

      /// @brief Class constructor
      FormFactor( FFScalingMode mode = FFScalingMode::DIPOLE );

      /// @brief Get the current Q^2 scaling mode
      inline FFScalingMode ff_scaling_mode() const;

      /// @brief Set the Q^2 scaling mode
      inline void set_ff_scaling_mode( FFScalingMode mode );

      /// @brief Convert a string to a Q^2 scaling mode
      static FFScalingMode ff_scaling_mode_from_string( const std::string& str );

      /// @brief Convert a Q^2 scaling mode to a string
      static std::string string_from_ff_scaling_mode( FFScalingMode mode );

      // Q^2 scaling functions

      /// @brief Dipole dependence
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Axial or vector mass parameter
      double dipole( double Q2, double M ) const;

      // Form factors

      /// @brief F1_n form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      double F1_n( double Q2, double M ) const;

      /// @brief F1_p form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      double F1_p( double Q2, double M ) const;

      /// @brief F1 form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      inline double F1( double Q2, double M ) const;

      /// @brief F2_n form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      double F2_n( double Q2, double M ) const;

      /// @brief F2_p form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      double F2_p( double Q2, double M ) const;

      /// @brief F2 form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      inline double F2( double Q2, double M ) const;

      /// @brief FS form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Vector mass parameter
      double FS( double Q2, double M ) const;

      /// @brief FA form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Axial mass parameter
      double FA( double Q2, double M ) const;

      /// @brief FP form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Axial mass parameter
      double FP( double Q2, double M ) const;

      /// @brief FT form factor
      /// @param Q2 The (minus) squared four-momentum transfer
      /// @param M Axial mass parameter
      double FT( double Q2, double M ) const;

    private:

      /// @brief Helper map for converting strings to Q^2 scaling modes and vice versa
      static std::map<FFScalingMode, std::string> ff_scaling_mode_map_;

      /// @brief The current Q^2 scaling mode
      FFScalingMode ff_scaling_mode_;

  };

  // Inline function definitions
  inline FormFactor::FFScalingMode FormFactor::ff_scaling_mode() const
    { return ff_scaling_mode_; }

  inline void FormFactor::set_ff_scaling_mode( FFScalingMode mode )
    { ff_scaling_mode_ = mode; }

  inline double FormFactor::F1( double Q2, double M ) const
    { return ( F1_p( Q2, M ) - F1_n( Q2, M ) ); }

  inline double FormFactor::F2( double Q2, double M ) const
    { return ( F2_p( Q2, M ) - F2_n( Q2, M ) ); }

}

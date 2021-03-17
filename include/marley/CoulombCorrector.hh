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

#pragma once

#include <map>
#include <string>

namespace marley {

  /// @brief Computes Coulomb correction factors for neutrino-nucleus
  /// differential cross sections
  class CoulombCorrector {

    public:

      /// @brief Enumerated type used to set the method for handling Coulomb
      /// corrections for CC nuclear reactions
      enum class CoulombMode { NO_CORRECTION, FERMI_FUNCTION, EMA, MEMA,
        FERMI_AND_EMA, FERMI_AND_MEMA };

      /// @param pdg_c Ejectile PDG code
      /// @param pdg_d Residue PDG code
      CoulombCorrector( int pdg_c, int pdg_d,
        CoulombMode mode = CoulombMode::FERMI_AND_MEMA );

      /// @brief Compute the
      /// <a href="https://en.wikipedia.org/wiki/Beta_decay#Fermi_function">
      /// Fermi function</a>
      /// @param beta_c <a
      /// href="http://scienceworld.wolfram.com/physics/RelativisticBeta.html">
      /// Dimensionless speed</a> of the ejectile
      double fermi_function( double beta_c ) const;

      /// Computes an approximate correction factor to account for
      /// effects of the Coulomb potential when calculating cross sections
      /// @param beta_rel_cd The relative speed of the final particles c and d
      /// (dimensionless)
      double coulomb_correction_factor( double beta_rel_cd ) const;

      /// Computes a Coulomb correction factor according to the effective
      /// momentum approximation. See J. Engel, Phys. Rev. C 57, 2004 (1998)
      /// @param beta_rel_cd The relative speed of the final particles c and d
      /// (dimensionless)
      /// @param[out] ok Flag that is set to false if subtracting the Coulomb
      /// potential pulls the event below threshold
      /// @param modified_ema If true, the modified EMA correction factor will
      /// be returned instead of that specified by the original EMA
      double ema_factor( double beta_rel_cd, bool& ok,
        bool modified_ema ) const;

      /// Return the method used by this reaction for handling Coulomb
      /// corrections
      inline CoulombMode coulomb_mode() const
        { return coulomb_mode_; }

      /// Set the method for handling Coulomb corrections for this reaction
      void set_coulomb_mode( CoulombMode mode );

      /// Convert a string to a CoulombMode value
      static CoulombMode coulomb_mode_from_string( const std::string& str );

      /// Convert a CoulombMode value to a string
      static std::string string_from_coulomb_mode( CoulombMode mode );

    protected:

      /// Helper map used by the methods to convert a CoulombMode value
      /// to and from a std::string
      static std::map< CoulombMode, std::string > coulomb_mode_string_map_;

      int pdg_c_; ///< Ejectile PDG code
      double mc_; ///< Ejectile mass

      int Zf_; ///< Residue atomic number
      int Af_; ///< Residue mass number

      /// @brief The method to use when computing Coulomb corrections
      CoulombMode coulomb_mode_;
  };

}

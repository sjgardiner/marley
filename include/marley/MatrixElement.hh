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

#pragma once

#include "marley/Error.hh"
#include "marley/Level.hh"

namespace marley {

  /// @brief A reduced nuclear matrix element that represents a transition
  /// caused by a neutrino-nucleus reaction
  class MatrixElement {

    public:

      /// @brief Enumerated type that represents the possible kinds of nuclear
      /// transitions recognized by MARLEY
      enum TransitionType {
        FERMI = 0, ///< The reduced Fermi matrix element is defined by
                   ///< @f$ \frac{ g_V^2 \big| \big< J_f \big\lVert
                   /// \mathcal{O}_\mathrm{F} \big\rVert J_i \big> \big|^2 }
                   /// { 2J_i + 1 }, @f$ where @f$J_i@f$ (@f$J_f@f$) is the
                   ///< initial (final) nuclear spin and @f$g_V@f$ is the
                   ///< vector coupling constant of the nucleon. The Fermi
                   ///< operator @f$ \mathcal{O}_\mathrm{F} @f$
                   ///< depends on the interaction process type
                   ///< (charged- vs. neutral-current) and is given by
                   ///< @f$ \mathcal{O}_\mathrm{F} = \begin{cases}
                   ///   \sum_{k = 1}^A t_{-}(k)
                   /// & \text{CC, }\nu\text{ projectile}
                   /// \\  \sum_{k = 1}^A t_{+}(k)
                   /// & \text{CC, }\bar{\nu}\text{ projectile}
                   /// \\ I & \text{NC}
                   /// \\ \end{cases} @f$ <br> where @f$ t_{\pm} @f$ are the
                   ///< isospin raising and lowering operators and @f$ I @f$ is
                   ///< the identity operator in isospace. MARLEY uses the
                   ///< convention where @f$ t_{-}\big|n\big> = \big|p\big> @f$.

        GAMOW_TELLER = 1 ///< The reduced Gamow-Teller matrix element is
                   ///< defined by @f$ \frac{ g_A^2 \big| \big< J_f \big\lVert
                   /// \sum_{k = 1}^A \boldsymbol{\sigma}(k)\Theta(k)
                   /// \big\rVert J_i \big> \big|^2 }{ 2J_i + 1 } @f$ where
                   ///< @f$ g_A @f$ is the axial coupling constant of the
                   ///< nucleon and @f$ \Theta(k) @f$ is an operator in
                   ///< isospace that depends on the type of scattering process:
                   ///< @f$ \Theta = \begin{cases}
                   /// t_{-} & \text{CC, }\nu\text{ projectile}
                   /// \\ t_{+} & \text{CC, }\bar{\nu}\text{ projectile}
                   /// \\ t_3   & \text{NC}
                   /// \\ \end{cases} @f$
      };

      /// @param level_energy Excitation energy (MeV) of the final-state nuclear
      /// level accessed by the matrix element
      /// @param strength Numerical value (dimensionless) of the matrix element
      /// @param type Type of nuclear transition (e.g., Fermi, Gamow-Teller)
      /// represented by the matrix element
      /// @param final_level Pointer to the Level object that represents the
      /// final nuclear level
      inline MatrixElement(double level_energy, double strength, TransitionType
        type, marley::Level* final_level = nullptr)
        : level_energy_(level_energy), strength_(strength), type_(type),
        final_level_(final_level) {}

      /// @brief Get the excitation energy (MeV) of the final-state nuclear
      /// level accessed by the matrix element
      inline double level_energy() const {
        // If this matrix element accesses a discrete final nuclear level,
        // then return that level's excitation energy
        if ( final_level_ ) return final_level_->energy();
        // Otherwise, we're in the ubound continuum, and we can just
        // use the tabulated value from the reaction matrix element data file
        else return level_energy_;
      }

      /// @brief Get the excitation energy (MeV) listed for this level
      /// in the reaction matrix element data file
      /// @details This value may differ from that returned by level_energy()
      /// for discrete nuclear levels. When final_level_ is not nullptr,
      /// level_energy() returns the excitation energy owned by final_level_,
      /// while tabulated_level_energy() returns the value of the level_energy_
      /// member variable. The latter is initialized from the reaction data
      /// file and may not match the discrete level data. To achieve consistency
      /// between the two sets of level energies, MARLEY overrides the
      /// values tabulated in the reaction data files with their closest
      /// matching discrete level energies. The level_energy_ member
      /// variable retains the reaction dataset value for validation purposes.
      /// For any physics calculation in MARLEY, level_energy() should be
      /// used instead of tabulated_level_energy().
      inline double tabulated_level_energy() const { return level_energy_; }

      /// @brief Get the numerical value (dimensionless) of the matrix element
      inline double strength() const { return strength_; }

      /// @brief Get the kind of nuclear transition (e.g., Fermi, Gamow-Teller)
      /// represented by the matrix element
      inline TransitionType type() const { return type_; }

      /// @brief Get a pointer to the final-state nuclear Level accessed by the
      /// matrix element, or nullptr if it is a transition to the unbound
      /// continuum
      inline const marley::Level* level() const { return final_level_; }

      /// @brief Get a pointer to the final-state nuclear Level accessed by the
      /// matrix element, or nullptr if it is a transition to the unbound
      /// continuum
      inline marley::Level* level() { return final_level_; }

      /// @brief Set the excitation energy (MeV) of the final-state nuclear
      /// level accessed by the matrix element
      inline void set_level_energy(double energy) { level_energy_ = energy; }

      /// @brief Set the numerical value (dimensionless) of the matrix element
      inline void set_strength(double strength) { strength_ = strength; }

      /// @brief Set the kind of nuclear transition (e.g., Fermi, Gamow-Teller)
      /// represented by the matrix element
      inline void set_type(TransitionType type) { type_ = type; }

      /// @brief Set the pointer to the final nuclear Level object accessed
      /// by the matrix element
      inline void set_level(marley::Level* lev) { final_level_ = lev; }

      /// @brief Compute the PDF for the CM frame scattering cosine
      /// @param cos_theta_c_cm Ejectile scattering cosine in the CM frame
      /// @param beta_c_cm Ejectile speed (dimensionless) as measured
      /// in the CM frame
      inline double cos_theta_pdf(double cos_theta_c_cm, double beta_c_cm) const
      {
        double pdf = 0.;
        if ( type_ == TransitionType::FERMI ) {
          pdf = 1. + beta_c_cm * cos_theta_c_cm;
        }
        else if ( type_ == TransitionType::GAMOW_TELLER ) {
          pdf = (3. - beta_c_cm * cos_theta_c_cm) / 3.;
        }
        else throw marley::Error("Unrecognized transition type encountered"
          " in marley::MatrixElement::cos_theta_pdf()");

        // Normalize to unit integral
        pdf *= 0.5;
        return pdf;
      }

      /// @brief Returns a string representation of the transition type for
      /// this matrix element
      inline std::string type_str() const {
        if ( type_ == TransitionType::FERMI ) return "Fermi";
        else if ( type_ == TransitionType::GAMOW_TELLER ) return "Gamow-Teller";
        else throw marley::Error( "Unrecognized transition type encountered"
          " in marley::MatrixElement::type_str()" );
      }

    protected:

      /// @brief Energy (MeV) of the final-state nuclear level accessed by this
      /// matrix element
      double level_energy_;

      /// @brief Numerical value of the matrix element (dimensionless)
      double strength_;

      /// @brief The kind of transition represented by this matrix element
      /// (Fermi, Gamow-Teller, etc.)
      TransitionType type_;

      /// @brief Pointer to the final Level object for a transition to a
      /// discrete nuclear level, or nullptr for an unbound state
      marley::Level* final_level_;
  };
}

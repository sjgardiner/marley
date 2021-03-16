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
#include <functional>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "marley/DecayScheme.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/MatrixElement.hh"
#include "marley/Reaction.hh"
#include "marley/StructureDatabase.hh"

namespace marley {

  class Generator;

  /// @brief A neutrino-nucleus reaction
  class NuclearReaction : public Reaction {

    public:

      /// @param pt Type of scattering process represented by this Reaction
      /// @param pdg_a Projectile PDG code
      /// @param pdg_b Target PDG code
      /// @param pdg_c Ejectile PDG code
      /// @param pdg_d Residue PDG code
      /// @param q_d Charge of the residue after the prompt 2->2 scatter
      /// represented by this NuclearReaction object
      NuclearReaction( ProcessType pt, int pdg_a, int pdg_b, int pdg_c,
        int pdg_d, int q_d );

      /// @brief Enumerated type used to set the method for handling Coulomb
      /// corrections for CC nuclear reactions
      enum class CoulombMode { NO_CORRECTION, FERMI_FUNCTION, EMA, MEMA,
        FERMI_AND_EMA, FERMI_AND_MEMA };

      inline virtual marley::TargetAtom atomic_target() const override final
        { return marley::TargetAtom( pdg_b_ ); }

      /// Produces a two-two scattering Event that proceeds via this reaction
      virtual std::shared_ptr< HepMC3::GenEvent > create_event(
        int particle_id_a, double KEa, marley::Generator& gen) const override;

      /// @brief Compute the
      /// <a href="https://en.wikipedia.org/wiki/Beta_decay#Fermi_function">
      /// Fermi function</a>
      /// @param beta_c <a
      /// href="http://scienceworld.wolfram.com/physics/RelativisticBeta.html">
      /// Dimensionless speed</a> of the ejectile
      double fermi_function( double beta_c ) const;

      /// @brief Get the minimum lab-frame kinetic energy (MeV) of the
      /// projectile that allows this reaction to proceed via a transition to
      /// the residue's ground state
      double threshold_kinetic_energy() const override;

      /// @brief Get the maximum possible excitation energy (MeV) of the
      /// final-state residue that is kinematically allowed
      /// @param KEa Projectile lab-frame kinetic energy (MeV)
      double max_level_energy( double KEa ) const;

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

      /// Computes the weak nuclear charge @f$ Q_W = N - (1
      /// - 4\sin^2\theta_W)Z @f$ for the target nucleus
      /// @details In the expression above, @f$ N @f$ (@f$Z@f$) is the
      /// neutron (proton) number of the target nucleus and
      /// @f$ \theta_W @f$ is the weak mixing angle.
      double weak_nuclear_charge() const;

      /// Return the method used by this reaction for handling Coulomb
      /// corrections
      inline CoulombMode coulomb_mode() const
        { return coulomb_mode_; }

      /// Set the method for handling Coulomb corrections for this reaction
      inline void set_coulomb_mode( CoulombMode mode )
        { coulomb_mode_ = mode; }

      /// Convert a string to a CoulombMode value
      static CoulombMode coulomb_mode_from_string( const std::string& str );

      /// Convert a CoulombMode value to a string
      static std::string string_from_coulomb_mode( CoulombMode mode );

    protected:

      /// Helper map used by the methods to convert a CoulombMode value
      /// to and from a std::string
      static std::map< CoulombMode, std::string > coulomb_mode_string_map_;

      /// @brief Creates the description string based on the
      /// PDG code values for the initial and final particles
      virtual void set_description();

      double md_gs_; ///< Ground state mass (MeV) of the residue

      int Zi_; ///< Target atomic number
      int Ai_; ///< Target mass number
      int Zf_; ///< Residue atomic number
      int Af_; ///< Residue mass number

      /// @brief Net charge of the residue (in units of the proton charge)
      /// following this reaction
      int q_d_;

      /// @brief Lab-frame kinetic energy of the projectile at threshold for
      /// this reaction (i.e., the residue is produced in its ground state, and
      /// all final-state particles are at rest in the CM frame)
      double KEa_threshold_;

      /// @brief The method to use when computing Coulomb corrections (if
      /// needed) for this reaction
      CoulombMode coulomb_mode_ = CoulombMode::FERMI_AND_MEMA;

  };

}

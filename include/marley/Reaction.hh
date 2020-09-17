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
#include <memory>
#include <string>
#include <vector>

#include "marley/MassTable.hh"
#include "marley/TargetAtom.hh"

namespace HepMC3 {
  class GenEvent;
}

namespace marley {

  class Generator;
  class Parity;
  class StructureDatabase;

  /// @brief Abstract base class that represents a 2 &rarr; 2 scattering
  /// reaction
  /// @details This class models a reaction of the form a + b &rarr; c + d.
  /// The projectile (particle a) is taken to have lab-frame kinetic energy KEa
  /// and to be traveling toward the target along the positive z direction.
  /// The target (particle b) is taken to be at rest in the lab frame.
  /// @todo Add check that projectile kinetic energy KEa >= 0. for all
  /// relevant member functions of Reaction and NuclearReaction
  class Reaction {

    public:

      virtual ~Reaction() = default;

      /// @brief Enumerated type describing the kind of scattering process
      /// represented by a Reaction
      enum ProcessType {
        NeutrinoCC = 0,     ///< Nuclear matrix elements contain @f$ t_{-} @f$
        AntiNeutrinoCC = 1, ///< Nuclear matrix elements contain @f$ t_{+} @f$
        NC = 2,             ///< Nuclear matrix elements contain @f$ t_{3} @f$
        NuElectronElastic = 3 ///< Neutrino-electron elastic scattering
      };

      /// @brief Compute the reaction's total cross section (MeV<sup> -2</sup>)
      /// @param pdg_a Projectile's PDG code
      /// @param KEa Lab-frame kinetic energy of the incident projectile
      /// @return %Reaction total cross section (MeV<sup> -2</sup>)
      /// @note Functions that override total_xs() should always return zero
      /// if pdg_a != pdg_a_.
      virtual double total_xs( int pdg_a, double KEa ) const = 0;

      /// @brief Differential cross section
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$
      /// (MeV<sup> -2</sup>)
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param cos_theta_c_cm CM frame scattering cosine of the ejectile
      /// @note Functions that override diff_xs() should always return zero
      /// if pdg_a != pdg_a_.
      virtual double diff_xs( int pdg_a, double KEa,
        double cos_theta_c_cm ) const = 0;

      /// @brief Create an event object for this reaction
      /// @param pdg_a PDG code for the incident projectile
      /// @param KEa Lab-frame kinetic energy of the projectile
      /// @param gen Reference to the Generator to use for random sampling
      /// @note Functions that override create_event() should throw an
      /// Error if pdg_a != pdg_a_.
      virtual std::shared_ptr< HepMC3::GenEvent > create_event( int pdg_a,
        double KEa, marley::Generator& gen ) const = 0;

      /// @brief Get a string that contains the formula for this reaction
      inline const std::string& get_description() const { return description_; }

      /// @brief Get the process type for this reaction
      inline ProcessType process_type() const { return process_type_; }

      static std::string proc_type_to_string( const ProcessType& pt );

      /// @brief Get the projectile PDG code
      inline int pdg_a() const { return pdg_a_; }

      /// @brief Get the target PDG code
      inline int pdg_b() const { return pdg_b_; }

      /// @brief Get the minimum lab-frame kinetic energy (MeV) of the
      /// projectile that allows this reaction to proceed via a transition to
      /// the residue's ground state
      virtual double threshold_kinetic_energy() const = 0;

      /// @brief Returns the target atom involved in this reaction
      /// @details For nuclear reactions, this is identical to the pdg_b_
      /// member variable. For electron reactions, it is distinct (since
      /// particle b is the initial struck electron).
      virtual marley::TargetAtom atomic_target() const = 0;

      /// Factory method called by JSONConfig to build
      /// Reaction objects given a file with matrix element data
      static std::vector< std::unique_ptr<Reaction> >
        load_from_file( const std::string& filename,
        marley::StructureDatabase& db );

      /// Function that returns the ejectile PDG code given the projectile
      /// PDG code and the ProcessType
      static int get_ejectile_pdg( int pdg_a, ProcessType proc_type );

    protected:

      int pdg_a_; ///< PDG code for the projectile
      int pdg_b_; ///< PDG code for the target
      int pdg_c_; ///< PDG code for the ejectile
      int pdg_d_; ///< PDG code for the residue

      double ma_; ///< Projectile mass (MeV)
      double mb_; ///< Target mass (MeV)
      double mc_; ///< Ejectile mass (MeV)

      // Mutable because, for the nuclear reaction case, we
      // may modify the final nuclear mass based on a sampled
      // excitation energy
      mutable double md_; ///< Residue mass (MeV)

      /// @brief String that contains a formula describing the reaction
      std::string description_;

      /// @brief Type of scattering process (CC, NC) represented by
      /// this reaction
      ProcessType process_type_;

      /// @brief Helper function that handles CM frame kinematics for the
      /// reaction
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param[out] s <a
      /// href="https://en.wikipedia.org/wiki/Mandelstam_variables"> Mandelstam
      /// s</a> (MeV<sup>2</sup>)
      /// @param[out] Ec_cm Ejectile total energy (MeV) in the CM frame
      /// @param[out] pc_cm Ejectile 3-momentum magnitude (MeV) in the CM frame
      /// @param[out] Ed_cm Residue total energy (MeV) in the CM frame
      void two_two_scatter( double KEa, double& s, double& Ec_cm,
        double& pc_cm, double& Ed_cm ) const;

      /// @brief Helper function that makes an event object.
      /// @details This function should be called by
      /// marley::Reaction::create_event() after CM frame scattering angles
      /// have been sampled for the ejectile. For reactions where the residue
      /// may be left in an excited state, the excitation energy should be
      /// recorded by supplying it as the final argument.
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param pc_cm Ejectile 3-momentum magnitude (MeV) in the CM frame
      /// @param cos_theta_c_cm Cosine of ejectile's CM frame polar angle
      /// @param phi_c_cm Ejectile's CM frame azimuthal angle (radians)
      /// @param Ec_cm Ejectile total energy (MeV) in the CM frame
      /// @param Ed_cm Residue total energy (MeV) in the CM frame
      /// @param E_level Residue excitation energy (MeV)
      /// @param twoJ Two times the residue spin
      /// @param P Intrinsic parity of the residue
      virtual std::shared_ptr< HepMC3::GenEvent > make_event_object(
        double KEa, double pc_cm, double cos_theta_c_cm, double phi_c_cm,
        double Ec_cm, double Ed_cm, double E_level, int twoJ,
        const marley::Parity& P ) const;

      /// Returns a vector of PDG codes for projectiles that participate
      /// in a particular ProcessType
      static const std::vector<int>& get_projectiles( ProcessType proc_type );
  };

}

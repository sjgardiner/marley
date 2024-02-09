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
#include <functional>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "marley/CoulombCorrector.hh"
#include "marley/DecayScheme.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/MatrixElement.hh"
#include "marley/NuclearReaction.hh"
#include "marley/StructureDatabase.hh"

namespace marley {

  class Generator;

  /// @brief A neutrino-nucleus reaction whose cross section is calculated
  /// according to the allowed approximation
  class AllowedNuclearReaction : public NuclearReaction {

    public:

      /// @param pt Type of scattering process represented by this Reaction
      /// @param pdg_a Projectile PDG code
      /// @param pdg_b Target PDG code
      /// @param pdg_c Ejectile PDG code
      /// @param pdg_d Residue PDG code
      /// @param q_d Charge of the residue after the prompt 2->2 scatter
      /// represented by this AllowedNuclearReaction object
      /// @param mat_els A vector of MatrixElement objects that should
      /// be used to compute cross sections for this AllowedNuclearReaction
      /// @param mode Indicates the method to use when computing Coulomb
      /// corrections for the reaction cross section
      AllowedNuclearReaction( ProcessType pt, int pdg_a, int pdg_b, int pdg_c,
        int pdg_d, int q_d,
        const std::shared_ptr<std::vector<marley::MatrixElement> >& mat_els,
        CoulombCorrector::CoulombMode mode );

      virtual std::shared_ptr< HepMC3::GenEvent > create_event(
        int particle_id_a, double KEa, marley::Generator& gen ) const override;

      /// @brief Sets the DecayScheme object to use for sampling excited levels
      /// in the residue
      void set_decay_scheme( marley::DecayScheme* scheme );

      /// @brief Total reaction cross section (MeV<sup> -2</sup>), including
      /// all kinematically-allowed final nuclear levels
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @return Reaction total cross section (MeV<sup> -2</sup>)
      /// @note This function returns 0. if pdg_a != pdg_a_.
      virtual double total_xs( int pdg_a, double KEa ) const override;

      /// @brief Differential cross section
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$
      /// (MeV<sup> -2</sup>) evaluated in the center-of-momentum frame
      /// for a transition to a particular final nuclear level
      /// @param mat_el MatrixElement object describing the transition to the
      /// final nuclear level
      /// @param KEa Lab-frame projectile kinetic energy (MeV)
      /// @param cos_theta_c_cm Ejectile scattering cosine as measured
      /// in the CM frame
      double diff_xs( const marley::MatrixElement& mat_el, double KEa,
        double cos_theta_c_cm ) const;

      /// @brief Differential cross section
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$
      /// (MeV<sup> -2</sup>) summed over all kinematically-allowed final
      /// nuclear levels
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param cos_theta_c_cm CM frame scattering cosine of the ejectile
      /// @note This function returns 0. if pdg_a != pdg_a_.
      virtual double diff_xs( int pdg_a, double KEa,
        double cos_theta_c_cm ) const;

      /// @brief Total cross section (MeV<sup> -2</sup>) for a given final
      /// nuclear level
      /// @param me MatrixElement object describing the transition to the
      /// final nuclear level
      /// @param KEa Lab-frame projectile kinetic energy (MeV)
      /// @param[out] beta_c_cm After this function is called, this variable
      /// will contain the (dimensionless) speed of the ejectile as measured
      /// in the center of momentum frame
      /// @param check_max_E_level Whether to check the reaction threshold
      /// via a call to max_level_energy(). If this argument is set to false,
      /// then the check will be skipped.
      virtual double total_xs( const marley::MatrixElement& me, double KEa,
        double& beta_c_cm, bool check_max_E_level = true ) const;

      /// Allows access to the owned vector of MatrixElement objects
      inline const std::vector<marley::MatrixElement>& matrix_elements() const
        { return *matrix_elements_; }

    protected:

      virtual void set_description() override;

      /// Helper function used by AllowedNuclearReaction::create_event()
      virtual std::shared_ptr< HepMC3::GenEvent > make_event_object(
        double KEa, double pc_cm, double cos_theta_c_cm, double phi_c_cm,
        double Ec_cm, double Ed_cm, double E_level, int twoJ,
        const marley::Parity& P ) const override;

      /// @brief Samples a polar angle cosine for the ejectile using
      /// the relevant portion of the reaction nuclear matrix element
      /// @param m_type Integer representing the type of transition
      /// (0 = Fermi, 1 = Gamow-Teller)
      /// @param beta_c_cm <a
      /// href="http://scienceworld.wolfram.com/physics/RelativisticBeta.html">
      /// Dimensionless speed</a> of the ejectile in the CM frame
      /// @param gen Reference to the Generator to use for random sampling
      double sample_cos_theta_c_cm( const marley::MatrixElement& matrix_el,
        double beta_c_cm, marley::Generator& gen ) const;

      /// Helper function for total_xs and summed_diff_xs()
      /// @param pdg_a PDG code for the projectile
      /// @param KEa Lab-frame projectile kinetic energy (MeV)
      /// @param cos_theta_c_cm CM frame scattering cosine for the ejectile
      /// @param level_xsecs If this pointer is not nullptr, then
      /// the std::vector<double> that it points to will be cleared
      /// and loaded with the partial cross sections to each individual nuclear
      /// level included in the sum. This feature is helpful when computing
      /// weights for sampling a nuclear transition in create_event().
      /// @param differential Whether the total cross section (false)
      /// or the differential cross section (true,
      /// @f$d\sigma/d\cos\theta_{c}^{\mathrm{CM}}@f$) should
      /// be summed over all kinematically accessible nuclear levels
      /// by this function
      /// @return The requested cross section (MeV<sup> -2</sup>)
      double summed_xs_helper( int pdg_a, double KEa, double cos_theta_c_cm,
        std::vector<double>* level_xsecs, bool differential ) const;

      /// @brief Matrix elements representing all of the possible nuclear
      /// transitions that may be caused by this reaction
      std::shared_ptr< std::vector<marley::MatrixElement> > matrix_elements_;

      /// @brief Object that handles calculations of Coulomb correction factors
      CoulombCorrector coulomb_corrector_;
  };

}

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
#include <string>

#include "MassTable.hh"
#include "Reaction.hh"

namespace marley {

  class Generator;

  // Represents a reaction in which an incident (anti)neutrino of any flavor
  // elastically scatters off of an electron originally bound to an atom with
  // atomic number Z_atom_. The cross section computed for this reaction is
  // computed per atom (scaled up from the single-electron cross section by the
  // atomic number Z_atom_, which is the same as the number of electron targets
  // bound to the atom). The small binding energy of the electrons is currently
  // neglected.
  class ElectronReaction : public marley::Reaction {

    public:

      ElectronReaction(int pdg_a, int target_atom_pdg);

      inline virtual marley::TargetAtom atomic_target() const override final
        { return atom_; }

      // Total reaction cross section (in MeV^(-2)) for an incident
      // projectile with lab-frame kinetic energy Ea
      virtual double total_xs(int pdg_a, double KEa) const override;

      // Differential cross section (MeV^(-2)) in the CM frame
      virtual double diff_xs(int pdg_a, double KEa, double cos_theta_c_cm)
        const;

      // Creates an event object for this reaction using the generator gen
      virtual std::shared_ptr< HepMC3::GenEvent > create_event(
        int particle_id_a, double KEa, marley::Generator& gen) const override;

      inline virtual double threshold_kinetic_energy() const override
        { return KEa_threshold_; }

      inline double g1() const { return g1_; }
      inline double g2() const { return g2_; }

    private:

      // Atomic target involved in this reaction
      marley::TargetAtom atom_;

      // Coupling constants to use for cross section calculations (updated
      // each time based on the projectile's particle ID)
      double g1_, g2_;

      // Threshold kinetic energy of the projectile
      double KEa_threshold_;

      /// @brief Helper function for the constructor.
      /// @details Sets the g1_ and g2_ member variables to the appropriate
      /// values based on the projectile PDG code (pdg_a_). Complains if the
      /// projectile PDG code isn't recognized.
      void set_coupling_constants();

      // Sample a CM frame scattering cosine for the ejectile. The first
      // argument is Mandelstam s.  This function assumes that the coupling
      // constants g1 and g2 have already been updated, so make sure to use it
      // only after calling determine_coupling_constants.
      double sample_cos_theta_c_cm(double s, marley::Generator& gen);
  };

}

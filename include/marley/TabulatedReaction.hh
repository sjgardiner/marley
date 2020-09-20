/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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
#include <memory>

// MARLEY includes
#include "marley/Reaction.hh"
#include "marley/TabulatedXSec.hh"

namespace marley {

  /// @brief Generates inclusive scattering events and computes cross
  /// sections using tabulated nuclear responses
  class TabulatedReaction : public Reaction {

    public:

      TabulatedReaction( Reaction::ProcessType pt, int pdg_a,
        int pdg_b, int pdg_c, int pdg_d, int q_d,
        const std::shared_ptr<TabulatedXSec>& txsec );

      inline virtual TargetAtom atomic_target() const
        { return TargetAtom( pdg_b_ ); }

      virtual Event create_event( int pdg_a, double KEa,
        Generator& gen ) const override;

      // TODO: write this correctly
      inline virtual double diff_xs( int pdg_a, double KEa,
        double cos_theta_c_cm ) const override { return 0.; }

      virtual double threshold_kinetic_energy() const override;
      virtual double total_xs( int pdg_a, double KEa ) const override;

    protected:

      /// @brief Lab-frame projectile kinetic energy threshold for a transition
      /// to the residue ground state
      double KEa_threshold_;

      /// @brief Ground-state mass of the residue
      double md_gs_;

      /// @brief Net charge of the residue
      int q_d_;

      /// @brief Helper object that handles cross section calculations
      std::shared_ptr< TabulatedXSec > xsec_;
  };

}

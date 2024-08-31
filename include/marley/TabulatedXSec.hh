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

// Standard library includes
#include <map>

// MARLEY includes
#include "marley/CoulombCorrector.hh"
#include "marley/ChebyshevInterpolatingFunction.hh"
#include "marley/LeptonFactors.hh"
#include "marley/Parity.hh"
#include "marley/Reaction.hh"
#include "marley/ResponseTable.hh"

using ProcType = marley::Reaction::ProcessType;

namespace marley {

  /// @brief Computes inclusive lepton-nucleus cross sections using
  /// tabulated nuclear response functions
  class TabulatedXSec {

    public:

      explicit TabulatedXSec( int target_pdg,
        marley::Reaction::ProcessType p_type,
        marley::CoulombCorrector::CoulombMode mode, double Delta );

      void add_table( const std::string& file_name );

      /// @brief Get the shift used to compute the effective energy transfer
      double Delta() const;

      /// @brief Returns true if this cross section represents a CC process
      /// or false otherwise
      inline bool is_cc() const  {
        bool is_cc = ( proc_type_ == ProcType::NeutrinoCC_Discrete || proc_type_ == ProcType::AntiNeutrinoCC_Discrete
         || proc_type_ == ProcType::NeutrinoCC_Continuum || proc_type_ == ProcType::AntiNeutrinoCC_Continuum );
        return is_cc;
      }

      /// @brief Simple struct representing a given multipole order, e.g., 2+
      struct MultipoleLabel {
        MultipoleLabel( unsigned J, marley::Parity Pi ) : J_( J ), Pi_( Pi ) {}
        unsigned J_;
        marley::Parity Pi_;

        /// Define the order of MultipoleLabel objects so that they
        /// be used as keys for a std::map
        inline bool operator<( const MultipoleLabel& ml ) const {
          if ( this->J_ != ml.J_ ) {
            return ( this->J_ < ml.J_ );
          }
          int p1 = static_cast<int>( this->Pi_ );
          int p2 = static_cast<int>( ml.Pi_ );
          return p1 < p2;
        }
      };

      double diff_xsec( int pdg_a, double KE_a, double omega, double cos_theta,
        const MultipoleLabel& ml ) const;

      void add_table( unsigned J, marley::Parity Pi,
        const std::string& file_name );

      inline const ResponseTable& get_table( const MultipoleLabel& ml ) const
        { return responses_.at( ml ); }

      inline const std::map< MultipoleLabel, ResponseTable >& get_table_map()
        const { return responses_; }

      double integral( int pdg_a, double KEa, const MultipoleLabel& ml,
        double& diff_max ) const;

      double integral( int pdg_a, double KEa ) const;

      void optimize( int pdg_a, double max_KEa );

    protected:

      // Helper function for integral that does the actual integration
      double compute_integral( int pdg_a, double KEa, const MultipoleLabel& ml,
        double& diff_max ) const;

      struct OptimizationMapKey {
        OptimizationMapKey( int pdg_a, MultipoleLabel ml )
          : pdg_a_( pdg_a ), ml_( ml ) {}

        /// Define the order of OptimizationMapKey objects so that they
        /// be used as keys for a std::map
        inline bool operator<( const OptimizationMapKey& omk ) const {
          if ( this->pdg_a_ != omk.pdg_a_ ) {
            return ( this->pdg_a_ < omk.pdg_a_ );
          }
          return this->ml_ < omk.ml_;
        }

        int pdg_a_;
        MultipoleLabel ml_;
      };

      struct OptimizationMapValue {
        OptimizationMapValue( ChebyshevInterpolatingFunction tot_xsec,
          ChebyshevInterpolatingFunction max_diff_xsec )
          : tot_xsec_( tot_xsec ), max_diff_xsec_( max_diff_xsec ) {}

        ChebyshevInterpolatingFunction tot_xsec_;
        ChebyshevInterpolatingFunction max_diff_xsec_;
      };

      /// Nuclide represented by this table of nuclear responses
      TargetAtom ta_;

      /// Kind of process for which the cross section will be computed
      marley::Reaction::ProcessType proc_type_;

      /// Method to use for computing Coulomb corrections to the cross section
      marley::CoulombCorrector::CoulombMode coulomb_mode_;

      /// PDG code of the nuclear residue produced by the reaction described by
      /// this cross section
      int pdg_d_;

      /// Tables of nuclear responses organized by multipole
      std::map< MultipoleLabel, ResponseTable > responses_;

      /// Storage for interpolating functions used to optimize calculations
      /// requiring numerical integration
      std::map< OptimizationMapKey, OptimizationMapValue > optimization_map_;

      /// Effective shift in the energy transfer that accounts for the mass
      /// difference between the ground state of the initial nucleus and the
      /// isobaric analog state in the daughter nucleus
      double Delta_;
  };

}

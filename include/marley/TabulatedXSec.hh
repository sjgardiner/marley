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
#include <map>

// MARLEY includes
#include "marley/Parity.hh"
#include "marley/Reaction.hh"
#include "marley/ResponseTable.hh"

namespace marley {

  /// @brief Computes inclusive lepton-nucleus cross sections using
  /// tabulated nuclear response functions
  class TabulatedXSec {

    public:

      explicit TabulatedXSec( int target_pdg,
        marley::Reaction::ProcessType p_type );

      void add_table( const std::string& file_name );

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

      class LeptonFactors {

        public:

          LeptonFactors() {}

          LeptonFactors( double vCC, double vLL, double vCL, double vT,
            double vTprime ) : vCC_( vCC ), vLL_( vLL ), vCL_( vCL ),
            vT_( vT ), vTprime_( vTprime ) {}

          inline double vCC() const { return vCC_; }
          inline double vLL() const { return vLL_; }
          inline double vCL() const { return vCL_; }
          inline double vT() const { return vT_; }
          inline double vTprime() const { return vTprime_; }

          inline void set_vCC(double val) { vCC_ = val; }
          inline void set_vLL(double val) { vLL_ = val; }
          inline void set_vCL(double val) { vCL_ = val; }
          inline void set_vT(double val) { vT_ = val; }
          inline void set_vTprime(double val) { vTprime_ = val; }

          // Scalar product with the corresponding nuclear responses
          inline double operator*(const ResponseTable::NuclearResponses& nr) {
            double product = 0.;
            product += this->vCC_ * nr.RCC();
            product += this->vLL_ * nr.RLL();
            product += this->vCL_ * nr.RCL();
            product += this->vT_ * nr.RT();
            product += this->vTprime_ * nr.RTprime();
            return product;
          }

        protected:
          double vCC_ = 0.;
          double vLL_ = 0.;
          double vCL_ = 0.;
          double vT_ = 0.;
          double vTprime_ = 0.;
      };

      double diff_xsec( int pdg_a, double KE_a, double omega, double cos_theta,
        const MultipoleLabel& ml );

      void add_table( unsigned J, marley::Parity Pi,
        const std::string& file_name );

      inline const ResponseTable& get_table( const MultipoleLabel& ml ) const
        { return responses_.at( ml ); }

      inline const std::map< MultipoleLabel, ResponseTable >& get_table_map()
        const { return responses_; }

      double integral( int pdg_a, double KEa, const MultipoleLabel& ml,
        double& diff_max );

      double integral( int pdg_a, double KEa );

    protected:

      /// Nuclide represented by this table of nuclear responses
      TargetAtom ta_;

      /// Kind of process for which the cross section will be computed
      marley::Reaction::ProcessType proc_type_;

      /// Tables of nuclear responses organized by multipole
      std::map< MultipoleLabel, ResponseTable > responses_;
  };

}

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

// MARLEY includes
#include "marley/NuclearResponses.hh"

namespace marley {

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

      inline void set_vCC( double val ) { vCC_ = val; }
      inline void set_vLL( double val ) { vLL_ = val; }
      inline void set_vCL( double val ) { vCL_ = val; }
      inline void set_vT( double val ) { vT_ = val; }
      inline void set_vTprime( double val ) { vTprime_ = val; }

      // Scalar product with the corresponding nuclear responses
      inline double operator*( const NuclearResponses& nr ) {
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

}

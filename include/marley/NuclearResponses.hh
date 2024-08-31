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

// Standard library includes
#include <fstream>
#include <functional>

namespace marley {

  class NuclearResponses {

    public:

      NuclearResponses() : rCC_( 0. ), rLL_( 0. ), rCL_( 0. ), rTvv_( 0. ),
        rTaa_( 0. ), rTprime_( 0. ) {}

      NuclearResponses( double rCC, double rLL, double rCL, double rTvv,
        double rTaa, double rTprime ) : rCC_( rCC ), rLL_( rLL ),
        rCL_( rCL ), rTvv_( rTvv ), rTaa_( rTaa ), rTprime_( rTprime ) {}

      inline double RCC() const { return rCC_; }
      inline double RLL() const { return rLL_; }
      inline double RCL() const { return rCL_; }
      inline double RT() const { return rTvv_ + rTaa_; }
      inline double RTvv() const { return rTvv_; }
      inline double RTaa() const { return rTaa_; }
      inline double RTprime() const { return rTprime_; }

      inline void set_RCC(double val) { rCC_ = val; }
      inline void set_RLL(double val) { rLL_ = val; }
      inline void set_RCL(double val) { rCL_ = val; }
      inline void set_RTvv(double val) { rTvv_ = val; }
      inline void set_RTaa(double val) { rTaa_ = val; }
      inline void set_RTprime(double val) { rTprime_ = val; }

      // Scalar operations
      NuclearResponses operator*(double d) const
        { return apply_to_members(d, [](double a, double b)
          -> double { return a * b; }); }

      NuclearResponses operator/(double d) const
        { return apply_to_members(d, [](double a, double b)
          -> double { return a / b; }); }

      // NuclearResponses operations
      NuclearResponses operator+(const NuclearResponses& other) const
        { return apply_to_members(other, [](double a, double b)
          -> double { return a + b; }); }

      NuclearResponses operator-(const NuclearResponses& other) const
        { return apply_to_members(other, [](double a, double b)
          -> double { return a - b; }); }

      NuclearResponses operator*(const NuclearResponses& other) const
        { return apply_to_members(other, [](double a, double b)
          -> double { return a * b; }); }

      NuclearResponses operator/(const NuclearResponses& other) const
        { return apply_to_members(other, [](double a, double b)
          -> double { return a / b; }); }

      bool operator==(const NuclearResponses& other) const {
        bool are_equal = this->rCC_ == other.rCC_;
        if ( are_equal ) are_equal = this->rLL_ == other.rLL_;
        if ( are_equal ) are_equal = this->rCL_ == other.rCL_;
        if ( are_equal ) are_equal = this->rTvv_ == other.rTvv_;
        if ( are_equal ) are_equal = this->rTaa_ == other.rTaa_;
        if ( are_equal ) are_equal = this->rTprime_ == other.rTprime_;
        return are_equal;
      }

      bool operator!=(const NuclearResponses& other) const {
        return !operator==( other );
      }

      inline void print(std::ostream& out) const {
        out << rCC_ << ' ' << rLL_ << ' ' << rCL_ << ' '
          << rTvv_ << ' ' << rTaa_ << ' ' << rTprime_;
      }

  protected:

    NuclearResponses apply_to_members( double d,
      const std::function<double(double,double)>& my_function ) const
    {
      // Create a new NuclearResponses object by applying the function to
      // each pair of elements
      NuclearResponses result;

      result.rCC_ = my_function( this->rCC_, d );
      result.rLL_ = my_function( this->rLL_, d );
      result.rCL_ = my_function( this->rCL_, d );
      result.rTvv_ = my_function( this->rTvv_, d );
      result.rTaa_ = my_function( this->rTaa_, d );
      result.rTprime_ = my_function( this->rTprime_, d );

      return result;
    }

    NuclearResponses apply_to_members( const NuclearResponses& other,
      const std::function<double(double,double)>& my_function ) const
    {
      // Create a new NuclearResponses object by applying the function to
      // each pair of elements
      NuclearResponses result;

      result.rCC_ = my_function( this->rCC_, other.rCC_ );
      result.rLL_ = my_function( this->rLL_, other.rLL_ );
      result.rCL_ = my_function( this->rCL_, other.rCL_ );
      result.rTvv_ = my_function( this->rTvv_, other.rTvv_ );
      result.rTaa_ = my_function( this->rTaa_, other.rTaa_ );
      result.rTprime_ = my_function( this->rTprime_, other.rTprime_ );

      return result;
    }

    double rCC_;
    double rLL_;
    double rCL_;
    double rTvv_;
    double rTaa_;
    double rTprime_;

  };

}

inline std::ostream& operator<<( std::ostream& out,
  const marley::NuclearResponses& np )
{
  np.print( out );
  return out;
}

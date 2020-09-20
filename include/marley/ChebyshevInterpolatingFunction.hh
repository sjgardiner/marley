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
#include <vector>

#include "marley/marley_utils.hh"

namespace marley {

  // Default Chebyshev grid size to use (when not using adaptive sizing)
  constexpr size_t DEFAULT_N_CHEBYSHEV = 64u;

  /// @brief Approximate representation of a 1D continuous function
  class ChebyshevInterpolatingFunction {

    public:

      // N = 0 triggers adaptive grid sizing
      ChebyshevInterpolatingFunction(const std::function<double(double)>& func,
        double x_min, double x_max, size_t N = 0);

      /// @brief Approximates the represented function using the barycentric
      /// formula
      inline double evaluate(double x) const {
        double numer = 0.;
        double denom = 0.;
        double w_j = -1.;
        for ( size_t j = 0; j <= N_; ++j ) {
          w_j = -w_j; // w_j = (-1)^(j)
          double x_j = Xs_.at( j );
          double f_j = Fs_.at( j );
          // If the requested x value is exactly equal to one
          // of the grid points where we previously evaluated the
          // function, then just return that
          if ( x_j == x ) return f_j;

          double temp = w_j / ( x - x_j );
          if ( j == 0 || j == N_ ) temp /= 2.;
          denom += temp;
          numer += temp * f_j;
        }

        double px = numer / denom;
        return px;
      }

      // @brief Returns the integral of this function on the interval [x_min_,
      // x_max_]
      inline double integral() const { return integral_; };

      inline const std::vector<double>& chebyshev_coeffs() const
        { return chebyshev_coeffs_; }

      inline const std::vector<double>& Fs() const
        { return Fs_; }

      inline const std::vector<double>& Xs() const
        { return Xs_; }

      inline int N() const { return N_; }

      ChebyshevInterpolatingFunction cdf() const;

      inline double x_min() const { return x_min_; }
      inline double x_max() const { return x_max_; }

    protected:

      /// @brief Default constructor used by cdf()
      inline ChebyshevInterpolatingFunction() {};

      /// @brief Maximum allowed value of the grid size parameter
      static constexpr size_t N_MAX_ = 65536u; // 16 adaptive iterations

      /// @brief Grid size parameter (N + 1 total points)
      size_t N_;

      /// @brief Lower edge of the grid
      double x_min_;

      /// @brief Upper edge of the grid
      double x_max_;

      /// @brief Approximate integral of the function on [x_min_, x_max_]
      double integral_;

      /// @brief Chebyshev points at which the function was evaluated
      std::vector<double> Xs_;

      /// @brief Function values at the grid points
      std::vector<double> Fs_;

      /// @brief Coefficients of the Chebyshev expansion of this function
      std::vector<double> chebyshev_coeffs_;

      // Helper arrays for discrete cosine transforms calculated using FFTPACK4
      std::vector<double> wsave_;
      std::vector<int> ifac_;

      /// @brief For a given N, returns the x position of the
      /// jth Chebyshev point (of the second kind)
      /// @todo If needed for speed, consider caching the std::cos evaluations
      /// used here
      inline double chebyshev_point(size_t j) const {
        // jth Chebyshev point (of the second kind) on [-1, 1]
        double x = std::cos( static_cast<double>(marley_utils::pi * j) / N_ );
        // Affine transformation to [x_min_, x_max_] (primed basis)
        double x_prime = (( x_max_ - x_min_ )*x + ( x_max_ + x_min_ )) / 2.;
        return x_prime;
      }

      void compute_integral();
  };

}

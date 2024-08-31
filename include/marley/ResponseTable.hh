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
#include <fstream>
#include <functional>
#include <vector>

// MARLEY includes
#include "marley/NuclearResponses.hh"

namespace marley {

  /// @brief A table of nuclear responses for computing inclusive
  /// cross sections
  class ResponseTable {
    public:

      ResponseTable() {}

      ResponseTable( std::shared_ptr< std::vector<double> >& wvec,
        std::shared_ptr< std::vector<double> >& qvec,
        std::shared_ptr< std::vector<NuclearResponses> >& resp )
        : wvec_( wvec ), qvec_( qvec ), responses_( resp ) {}

      /// Retrieve the minimum &omega; value
      inline double w_min() const { return wvec_->front(); }

      /// Retrieve the maximum &omega; value
      inline double w_max() const { return wvec_->back(); }

      /// Retrieve the minimum q value
      inline double q_min() const { return qvec_->front(); }

      /// Retrieve the maximum q value
      inline double q_max() const { return qvec_->back(); }

      /// Access the energy transfer grid points
      inline const std::vector<double>& w_grid() const { return *wvec_; }

      /// Access the 3-momentum transfer grid points
      inline const std::vector<double>& q_grid() const { return *qvec_; }

      /// Access the nuclear responses
      inline const std::vector<NuclearResponses>& get_responses() const
        { return *responses_; }

      /// Calculates the index in the vector of responses that
      /// corresponds to a given set of &omega; and q indices
      /// \param[in] iw Index of the desired grid point on the &omega; axis
      /// \param[in] iq Index of the desired grid point on the q axis
      inline size_t response_index( size_t iw, size_t iq ) const {
        size_t num_q_points = qvec_->size();
        size_t r_idx = ( num_q_points * iw ) + iq;
        return r_idx;
      }

      /// Use simple bilinear interpolation to compute a set of nuclear
      /// responses
      NuclearResponses interpolate(double w, double q) const;

    protected:

      /// Grid points along the w-axis
      std::shared_ptr< std::vector<double> > wvec_;

      /// Grid points along the q-axis
      std::shared_ptr< std::vector<double> > qvec_;

      /// Sets of nuclear responses evaluated at each grid point
      std::shared_ptr< std::vector<NuclearResponses> > responses_;

      /// Determines the indices for the two gridpoints surrounding a requested
      /// one-dimensional coordinate. If the coordinate is outside of the grid,
      /// then this function returns the two closest grid points.
      /// \param[in] vec A vector of grid point coordinates
      /// \param[in] val The requested x or y coordinate
      /// \param[out] lower_index The index of the closest grid point less than
      /// or equal to the requested value, or the lower of the two nearest grid
      /// points if the value falls outside of the grid
      /// \param[out] upper_index The index of the closest grid point greater
      /// than the requested value, or the higher of the two nearest grid
      /// points if the value falls outside of the grid
      /// \return true if the requested value is within the grid, or false
      /// otherwise
      bool get_neighbor_indices( const std::vector<double>& vec,
        double val, size_t& lower_index, size_t& upper_index ) const;
  };

}

inline std::ostream& operator<<(std::ostream& out,
  const marley::ResponseTable& rt)
{
  for ( size_t iw = 0u; iw < rt.w_grid().size(); ++iw ) {
    for ( size_t iq = 0u; iq < rt.q_grid().size(); ++iq ) {
      out << rt.w_grid().at( iw ) << "  " << rt.q_grid().at( iq ) << "  ";
      size_t ir = rt.response_index( iw, iq );
      out << rt.get_responses().at( ir ) << '\n';
    }
  }
  return out;
}

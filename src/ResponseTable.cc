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

#include "marley/Error.hh"
#include "marley/ResponseTable.hh"

marley::ResponseTable::NuclearResponses
  marley::ResponseTable::interpolate(double w, double q) const
{
  // Return a dummy value if the point lies outside the grid boundaries
  if ( w < this->w_min() || w > this->w_max()
    || q < this->q_min() || q > this->q_max() )
  {
    throw marley::Error( "Outside of grid bounds in marley::ResponseTable::"
      "interpolate()" );
    return marley::ResponseTable::NuclearResponses();
  }

  // Find the indices of the grid points on either side of the desired w and
  // q values.
  size_t iw_low, iw_hi, iq_low, iq_high;
  this->get_neighbor_indices( *wvec_, w, iw_low, iw_hi );
  this->get_neighbor_indices( *qvec_, q, iq_low, iq_high );

  // Get the w and q values corresponding to the nearest-neighbor grid
  // points found above
  double w1 = wvec_->at( iw_low );
  double w2 = wvec_->at( iw_hi );
  double q1 = qvec_->at( iq_low );
  double q2 = qvec_->at( iq_high );

  // Retrieve the nuclear responses at each of the four grid points of
  // interest
  const auto& r11 = responses_->at( this->response_index(iw_low, iq_low) );
  const auto& r21 = responses_->at( this->response_index(iw_hi, iq_low) );
  const auto& r12 = responses_->at( this->response_index(iw_low, iq_high) );
  const auto& r22 = responses_->at( this->response_index(iw_hi, iq_high) );

  // Perform the interpolation (first q, then w)
  NuclearResponses r1 = r11 * (q2 - q)/(q2 - q1) + r12 * (q - q1)/(q2 - q1);
  NuclearResponses r2 = r21 * (q2 - q)/(q2 - q1) + r22 * (q - q1)/(q2 - q1);
  NuclearResponses r  = r1  * (w2 - w)/(w2 - w1) + r2  * (w - w1)/(w2 - w1);

  // Return the result
  return r;
}

bool marley::ResponseTable::get_neighbor_indices(
  const std::vector<double>& vec, double val, size_t& lower_index,
  size_t& upper_index ) const
{
  if ( vec.size() < 2u ) {
    throw marley::Error( "Vector with fewer than two entries passed to"
      " marley::ResponseTable::get_neighbor_indices()" );
  }

  bool within = true;

  auto begin = vec.cbegin();
  auto end = vec.cend();

  // std::lower_bound returns an iterator to the first element of the
  // container which is not less than the supplied value
  auto not_less_point = std::lower_bound( begin, end, val );

  decltype( begin ) lower_point;

  // Check whether the requested grid point is within the grid limits
  if ( not_less_point == begin ) {
    lower_point = begin;
    // first element of vec > val
    if ( *begin != val ) within = false;
  }
  else if ( not_less_point == end ) {
    // last element of vec < val
    within = false;
    lower_point = end - 2;
  }
  else {
    // x is within the grid limits
    lower_point = not_less_point - 1;
  }

  lower_index = std::distance( begin, lower_point );
  upper_index = lower_index + 1;

  return within;
}

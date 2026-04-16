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

// Standard library includes
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"
#include "marley/DecayScheme.hh"
#include "marley/Generator.hh"
#include "marley/marley_kinematics.hh"
#include "marley/Logger.hh"
#include "marley/HauserFeshbachDecay.hh"

// Returns a pointer to the level owned by this decay scheme object
// that has the closest excitation energy to E_level (E_level
// has units of MeV).
marley::Level* marley::DecayScheme::get_pointer_to_closest_level(
  double E_level)
{
  // If this decay scheme doesn't own any levels, return nullptr immediately
  size_t num_levels = levels_.size();
  if (num_levels == 0) return nullptr;

  // Search for the level whose energy is closest to the given value of E_level
  size_t e_index = level_lower_bound_index(E_level);

  if (e_index == num_levels) {
    // The given energy is greater than every level energy in our decay scheme.
    // We will therefore assume that the desired level is the highest level.
    // Its index is given by one less than the number of elements in the sorted
    // vector, so subtract one from our previous result.
    --e_index;
  }
  else if (e_index > 0) {
    // If the calculated index does not correspond to the first element, we
    // still need to check which of the two levels found (one on each side) is
    // really the closest. Do so and reassign the index if needed.
    if (std::abs(E_level - levels_.at(e_index)->energy())
      > std::abs(E_level - levels_.at(e_index - 1)->energy()))
    {
      --e_index;
    }
  }

  // Return a pointer to the selected level object
  return levels_.at(e_index).get();
}

int marley::DecayScheme::pdg() const {
  return marley_utils::get_nucleus_pid( Z_, A_ );
}

void marley::DecayScheme::do_cascade( marley::Level& initial_level,
  HepMC3::GenEvent& event, marley::Generator& gen,
  std::shared_ptr< HepMC3::GenParticle >& residue )
{
  int qIon = marley_hepmc3::get_particle_charge( *residue );

  MARLEY_LOG_DEBUG() << "Beginning gamma cascade at level with energy "
    << initial_level.energy() << " MeV";

  marley::Level* p_current_level = &initial_level;

  const marley::MassTable& mt = marley::MassTable::Instance();

  // Initialize some variables used in the gamma cascade loop below
  bool cascade_finished = false;
  double gamma_branching_ratio = 0.;
  double level_total_width = 0.;

  while ( !cascade_finished ) {

    // Randomly select a gamma to produce while storing its branching ratio
    const marley::Gamma* p_gamma = p_current_level->sample_gamma( gen,
      &gamma_branching_ratio );

    if ( !p_gamma ) {
      MARLEY_LOG_DEBUG() << "  this level does not have any gammas";
      cascade_finished = true;
    }
    else {

      // Get the total decay width (MeV) for the initial level
      level_total_width = marley_utils::hbar * marley_utils::log_2
        / p_current_level->half_life();

      // Update the current level now that the gamma has been emitted
      p_current_level = p_gamma->end_level();
      if ( !p_current_level ) {
        throw marley::Error( "This gamma does not have an end level."
          " Cannot continue cascade." );
      }
      MARLEY_LOG_DEBUG() << std::setprecision(15) << std::scientific
        << "  emitted gamma with energy "
        << p_gamma->energy() << " MeV. New level has energy "
        << p_current_level->energy() << " MeV.";
      MARLEY_LOG_DEBUG() << "gamma energy = " << p_gamma->energy();

      // Get the excitation energy of the end level. This will be added to
      // the ground state mass of the nucleus to determine its
      // post-gamma-emission mass.
      double Exf = p_current_level->energy();

      // Create new particle objects to represent the emitted gamma and
      // recoiling nucleus
      auto gamma = marley_hepmc3::make_particle( marley_utils::PHOTON,
        marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, 0.0 );

      int pdg = marley_utils::get_nucleus_pid( Z_, A_ );

      auto nucleus = marley_hepmc3::make_particle( pdg,
        marley_hepmc3::NUHEPMC_INTERMEDIATE_RESIDUE_STATUS,
        mt.get_atomic_mass(pdg) + Exf
        - qIon*mt.get_particle_mass(marley_utils::ELECTRON) );

      // Create a new binary decay vertex
      auto decay_vtx = std::make_shared< HepMC3::GenVertex >();
      decay_vtx->set_status( marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX );

      decay_vtx->add_particle_in( residue );
      decay_vtx->add_particle_out( gamma );
      decay_vtx->add_particle_out( nucleus );

      // Sample a decay time (MeV^{-1}) for emission of the chosen gamma-ray
      // and store this timing information in the new binary decay vertex
      double gamma_partial_width = gamma_branching_ratio * level_total_width;
      marley_hepmc3::store_decay_time( gamma_partial_width, gen, decay_vtx,
        residue );

      // Add the decay vertex to the event record
      event.add_vertex( decay_vtx );

      // We can set the charge attribute now that the daughter nucleus
      // belongs to the decay vertex (and thus the parent event)
      marley_hepmc3::set_particle_charge( *nucleus, qIon );

      // Sample a direction assuming that the gammas are emitted isotropically
      // in the nucleus's rest frame.
      // sample from [-1, 1]
      double gamma_cos_theta = gen.uniform_random_double( -1.0, 1.0, true );
      // sample from [0, 2*pi)
      double gamma_phi = gen.uniform_random_double( 0., 2.*marley_utils::pi,
        false );

      // Determine the final energies and momenta for the recoiling nucleus and
      // emitted gamma ray. Store them in the final state particle objects.
      marley_kinematics::two_body_decay( residue, gamma, nucleus,
        gamma_cos_theta, gamma_phi );

      // Update the residue for this event to take into account changes from
      // gamma ray emission
      residue.swap( nucleus );
    }
  }

  MARLEY_LOG_DEBUG() << "Finished gamma cascade at level with energy "
    << p_current_level->energy();

  residue->set_status( marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS );
}

marley::DecayScheme::DecayScheme( int Z, int A ) : Z_( Z ), A_( A )
{
}

marley::DecayScheme::DecayScheme( int Z, int A, const std::string& filename,
  marley::DecayScheme::FileFormat ff ) : Z_( Z ), A_( A )
{
  this->parse( filename, ff );
}

void marley::DecayScheme::parse_talys( const std::string& filename ) {
  // First line in a TALYS level dataset has fortran
  // format (2i4, 2i5, 56x, i4, a2)
  // General regex for this line:
  // std::regex nuclide_line("[0-9 ]{18} {56}[0-9 ]{4}.{2}");
  std::string nuc_id = marley_utils::nuc_id( Z_, A_ );

  // Make the last character of the nuc_id lowercase to follow the TALYS
  // convention
  nuc_id.back() = tolower( nuc_id.back() );

  const std::regex nuclide_line( "[0-9 ]{18} {57}" + nuc_id );

  // Open the TALYs level data file for parsing
  std::ifstream file_in( filename );

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if ( !file_in.good() ) throw marley::Error( "Could not read from the"
    " TALYS data file " + filename );

  // String to store the current line of the TALYS file during parsing
  std::string line;
  bool found_decay_scheme = false;

  while ( std::getline(file_in, line) ) {
    if ( std::regex_match(line, nuclide_line) ) {
      found_decay_scheme = true;
      break;
    }
  }

  if ( !found_decay_scheme ) throw marley::Error( "Gamma decay scheme data"
    " (adopted levels, gammas) for " + marley_utils::nucid_to_symbol( nuc_id )
    + " could not be found in the TALYS data file " + filename );

  MARLEY_LOG_DEBUG() << "Gamma decay scheme data for " + nuc_id
    << " found. Using TALYS dataset ";
  MARLEY_LOG_DEBUG() << line;

  // Dummy integer and number of excited levels for this nuclide
  int dummy, num_excited_levels;

  // Read in the number of excited levels from the first line of data
  std::istringstream iss( line );
  iss >> dummy >> dummy >> dummy >> num_excited_levels;

  for ( int l_idx = 0; l_idx <= num_excited_levels; ++l_idx ) {

    // Get the next line of the file. This will be a discrete level record
    std::getline( file_in, line );

    // Load the new line into our istringstream object for parsing. Reset
    // the stream so that we start parsing from the beginning of the string.
    iss.str( line );
    iss.clear();

    // Read in this level's index, energy, spin, parity, and
    // number of gamma transitions
    int level_num, pi, num_gammas;
    double level_energy, spin, half_life;
    iss >> level_num >> level_energy >> spin >> pi >> num_gammas >> half_life;

    // Compute two times the spin so that we can represent half-integer
    // nuclear level spins as integers
    int twoJ = std::round( 2 * spin );

    // Create a parity object to use when constructing the level
    marley::Parity parity = marley::Parity( pi );

    // Construct a new level object and add it to the decay scheme. Get
    // a pointer to the newly-added level
    marley::Level& current_level = add_level( marley::Level(level_energy,
      twoJ, parity, half_life) );

    for ( int g_idx = 0; g_idx < num_gammas; ++g_idx ) {

      // Get the next line of the file. This will be a gamma record
      std::getline( file_in, line );

      // Load the new line into our istringstream object for parsing. Reset
      // the stream so that we start parsing from the beginning of the string.
      iss.str( line );
      iss.clear();

      // Read in the index of the final level and branching ratio
      // for this gamma transition
      int gamma_final_level_num;
      double br;
      iss >> gamma_final_level_num >> br;

      // Process this gamma if it has a nonvanishing branching ratio
      if ( br > 0. ) {

        marley::Level* final_level = levels_.at( gamma_final_level_num ).get();

        // Compute the gamma ray's energy in MeV by subtracting the energy
        // of the final level from the energy of the initial level
        double gamma_energy = level_energy - final_level->energy();

        // Create the new Gamma object for the current level
        current_level.add_gamma( gamma_energy, br, final_level );
      }
    }
  }

  file_in.close();
}

void marley::DecayScheme::print_report( std::ostream& ostr ) const {
  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for ( const auto& lev : levels_ ) {
    int twoj = lev->twoJ();
    std::string spin = std::to_string( twoj / 2 );
    // If 2*J is odd, then the level has half-integer spin
    if ( twoj % 2 ) spin += "/2";
    marley::Parity parity = lev->parity();

    ostr << "Level at " << lev->energy() << " MeV has spin-parity "
      << spin << parity << " and half-life " << lev->half_life() << " s\n";

    std::vector< marley::Gamma >& gammas = lev->gammas();

    // Cycle through each of the gammas owned by the current level
    // (according to the ENSDF specification, these will already be
    // sorted in order of increasing energy)
    for ( const auto& g : gammas ) {
      ostr << "  has a gamma with energy " << g.energy() << " MeV";
      ostr << " (transition to level at "
        << g.end_level()->energy() << " MeV)" << '\n';
      ostr << "    and relative intensity " << g.relative_intensity() << '\n';
    }
  }
}

void marley::DecayScheme::print_latex_table( std::ostream& ostr ) {

  std::string nuc_id = marley_utils::nuc_id( Z_, A_ );

  std::string caption_beginning =
    std::string("{\\textbf{Levels") +
    " and $\\boldsymbol{\\gamma}$ transitions \n for " +
    "\\isotope[\\boldsymbol{" +
    marley_utils::trim_copy(nuc_id.substr(0,3)) +
    "}]{\\textbf{" + nuc_id.substr(3,1) +
    marley_utils::trim_copy(
      marley_utils::to_lowercase(nuc_id.substr(4,1))) +
    "}} \n";// from file " + filename + " ";

  ostr << marley_utils::latex_table_1;

  ostr << caption_beginning + "}}\\\\\n";

  ostr << marley_utils::latex_table_2;

  ostr << caption_beginning + " -- \\textit{continued}}} \\\\\n";

  ostr << marley_utils::latex_table_3;

  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for (const auto& lev : levels_ ) {

    std::string sp = lev->spin_parity_string();

    ostr << lev->energy() << " & " << sp  << " & ";

    const auto& gammas = lev->gammas();

    // If there aren't any gammas for this level, finish writing
    // the current row of the table. Add extra space between this
    // level and the next one.
    if ( gammas.empty() ) {
      ostr << " &  &";
      // If this is the last row of the table, don't add extra space.
      if (lev == levels_.back()) ostr << '\n';
      else ostr << " \\\\ \\addlinespace[\\ExtraRowSpace]\n";
    }

    // Cycle through each of the gammas owned by the current level
    for ( const auto& g : gammas ) {
      // If this is not the first gamma, add empty columns
      // for the level energy and spin-parity
      if ( &g != &gammas.front() ) ostr << " & & ";
      // Output information about the current gamma
      ostr << g.energy() << " & " << g.relative_intensity()
        << " & " << g.end_level()->energy();
      // Add vertical space after the final gamma row. Also prevent page breaks
      // in the middle of a list of gammas by outputting a star at the end of
      // each row except the final gamma row.
      if ( &g == &gammas.back() ) {
	// Don't add the extra row space for the very last row in the table
	if ( lev == levels_.back() ) ostr << '\n';
        else ostr << " \\\\ \\addlinespace[\\ExtraRowSpace]" << '\n';
      }
      else ostr << " \\\\*" << '\n';
    }
  }
  ostr << marley_utils::latex_table_4 << '\n';
}

// Finds the index for the first level with excitation energy not less than Ex
size_t marley::DecayScheme::level_lower_bound_index( double Ex ) {
  const auto E_begin = marley::Level::make_energy_iterator( levels_.cbegin() );
  const auto E_end = marley::Level::make_energy_iterator( levels_.cend() );

  const auto closest_E_iter = std::lower_bound( E_begin, E_end, Ex );
  return std::distance( E_begin, closest_E_iter );
}

// Adds a new level to the decay scheme and returns a reference to it
marley::Level& marley::DecayScheme::add_level( const marley::Level& level )
{
  // Compute the numerical index for where we will insert the new level
  size_t index = level_lower_bound_index( level.energy() );

  // Insert the new level into the decay scheme
  levels_.insert( levels_.begin() + index,
    std::make_unique< marley::Level >(level) );

  // Return a reference to the newly-added level
  return *levels_.at( index );
}

void marley::DecayScheme::print( std::ostream& out ) const {

  size_t num_levels = levels_.size();

  out << Z_ << ' ' << A_ << ' ' << num_levels << '\n';

  for ( const auto& lev : levels_ ) {
    out << "  " << lev->energy() << ' ' << lev->twoJ() << ' '
      << lev->parity() << ' ' << lev->gammas().size() << '\n';
    for ( const auto& g : lev->gammas() ) {
      out << "    " << g.energy() << ' ' << g.relative_intensity();

      const auto cit = std::find_if( levels_.cbegin(), levels_.cend(),
        [&g]( const std::unique_ptr< marley::Level >& l )
        -> bool { return l.get() == g.end_level(); } );

      int level_f_idx = -1;
      if ( cit != levels_.cend() ) {
        level_f_idx = std::distance(levels_.cbegin(), cit);
      }
      out << " " << level_f_idx << '\n';
    }
  }
}

void marley::DecayScheme::read_from_stream( std::istream& in ) {

  levels_.clear();

  int num_levels;
  in >> Z_ >> A_ >> num_levels;

  // If we had trouble parsing the decay scheme header, then
  // just return the stream without doing anything else.
  if ( !in ) return;

  double energy, ri, half_life;
  int two_j, num_gammas, level_f_idx;
  marley::Parity pi;

  for ( int i = 0; i < num_levels; ++i ) {
    in >> energy >> two_j >> pi >> num_gammas >> half_life;

    marley::Level& l = add_level( marley::Level(energy, two_j, pi, half_life) );
    for ( int j = 0; j < num_gammas; ++j ) {
      in >> energy >> ri >> level_f_idx;
      l.add_gamma( energy, ri, levels_.at(level_f_idx).get() );
    }
  }
}

void marley::DecayScheme::parse( const std::string& filename,
  marley::DecayScheme::FileFormat ff )
{
  // Parse the data file using the appropriate format
  switch (ff) {

    case FileFormat::native:
      this->parse_native( filename );
      break;

    case FileFormat::talys:
      this->parse_talys( filename );
      break;

    // Add more data file formats as needed

    default:
      throw marley::Error( "Unsupported file format passed to"
        " marley::DecayScheme constructor." );
  }
}

void marley::DecayScheme::parse_native( const std::string& filename ) {

  // Open the level data file for parsing
  std::ifstream file_in( filename );

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if ( !file_in.good() ) throw marley::Error( "Could not read from the"
    " data file " + filename );

  read_from_stream( file_in );

  file_in.close();
}

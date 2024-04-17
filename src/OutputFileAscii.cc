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

// Standard library includes
#include <limits>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

// MARLEY includes
#include "marley/OutputFileAscii.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"

namespace {

  constexpr char DUMMY_CHAR = 'a';

  // Moves to the end of an ASCII-format file containing HepMC3 events and
  // opened as the input std::fstream. Backs up until the start of the last
  // HepMC3 event stored in the file. The stream is left in a state that is
  // ready for reading in the last event using, e.g., HepMC3::ReaderAscii.
  void seek_to_last_genevent( std::fstream& stream ) {
    char c = DUMMY_CHAR;
    char old_c = DUMMY_CHAR;
    stream.seekg( 0, std::ios::end );
    std::streampos size = stream.tellg();
    for ( int i = 1; i <= size; ++i ) {
      stream.seekg( -i, std::ios::end );
      old_c = c;
      stream.get( c );
      if ( c == '\n' && old_c == 'E' ) {
        stream.seekg( -i + 1, std::ios::end );
        stream.clear();
        break;
      }
    }
  }

}

marley::OutputFileAscii::OutputFileAscii( const marley::JSON& output_config )
  : marley::OutputFile( output_config )
{
  format_ = Format::ASCII;
  
  this->open();
  reader_ = std::make_shared< HepMC3::ReaderAscii >( stream_ );
  writer_ = std::make_shared< HepMC3::WriterAscii >( stream_ );

  // Ensure that all floating-point values are output with full precision
  writer_->set_precision( std::numeric_limits<double>::max_digits10 );
}

void marley::OutputFileAscii::open() {

  bool file_exists = check_if_file_exists( name_ );

  auto open_mode_flag = std::ios::in | std::ios::out | std::ios::trunc;

  if ( mode_ == Mode::OVERWRITE && file_exists && !force_ ) {
    bool overwrite = marley_utils::prompt_yes_no( "Overwrite file " + name_ );
    if ( !overwrite ) {
      MARLEY_LOG_INFO() << "Cancelling overwrite of output file \""
        << name_ << '\"';
      open_mode_flag = std::ios::in | std::ios::out;
      mode_ = Mode::RESUME;
    }
  }

  if ( mode_ == Mode::RESUME ) {
    if ( !file_exists ) throw marley::Error( "Cannot resume run. Could"
      " not open the file \"" + name_ + '\"' );
    else open_mode_flag = std::ios::in | std::ios::out;
  }
  else if ( mode_ != Mode::OVERWRITE ) {
    throw marley::Error( "Unrecognized file mode encountered in"
      " OutputFileAscii::open()" );
  }

  stream_.open( name_, open_mode_flag );

}

bool marley::OutputFileAscii::resume( std::unique_ptr<marley::Generator>& gen,
  long& num_previous_events )
{
  if ( mode_ != Mode::RESUME ) {
    throw marley::Error( "Cannot call OutputFileAscii::resume() for an output"
      " mode other than \"resume\"" );
    return false;
  }

  MARLEY_LOG_INFO() << "Continuing previous run from the file " << name_;

  // Create a temporary event to use for storage while parsing the file
  auto evt = std::make_shared< HepMC3::GenEvent >();

  // Read back the first event from the file so that we can retrieve the run
  // information
  stream_.seekg( 0, std::ios::beg );
  bool read_ok = reader_->read_event( *evt );

  auto run_info = evt->run_info();

  if ( !read_ok || !run_info ) {
    throw marley::Error( "Failed to retrieve run information from the file "
      + name_ );
    return false;
  }

  // Extract the prior generator JSON configuration from the RunInfo
  auto config_str = run_info->attribute< HepMC3::StringAttribute >(
    "MARLEY.JSONconfig" );

  if ( !config_str ) {
    throw marley::Error( "Failed to retrieve generator configuration from the"
      " file " + name_ );
    return false;
  }

  auto seed_str = run_info->attribute< HepMC3::StringAttribute >(
    "MARLEY.RNGseed" );

  if ( !seed_str ) {
    throw marley::Error( "Failed to load previous random number"
      " generator seed from the file " + name_ );
    return false;
  }

  // Parse the prior generator configuration into a JSON object
  auto json_config = marley::JSON::load( config_str->value() );

  // Move to the beginning of the last HepMC3 event in the file
  seek_to_last_genevent( stream_ );

  // Parse the last HepMC3 event so that we can retrieve the generator state
  read_ok = reader_->read_event( *evt );

  auto state_str = evt->attribute< HepMC3::StringAttribute >(
    "MARLEY.GeneratorState" );

  if ( !read_ok || !state_str ) {
    throw marley::Error( "Failed to retrieve generator state from the file "
      + name_ );
    return false;
  }

  // We assume that the events in the file were originally written out using
  // the convention in the marley executable: the event number of the last
  // event in the file is equal to the total number of events.
  num_previous_events = evt->event_number();

  // We're done retrieving the information. Restore the Generator to its
  // previous state
  gen = this->restore_generator( json_config );
  gen->seed_using_state_string( state_str->value() );

  MARLEY_LOG_INFO() << "The previous run was initialized using"
    << " the random number generator seed " << seed_str->value();

  // TODO: maybe use seekg() to back up to remove HepMC3 footer after the last
  // event

  return true;
}

int_fast64_t marley::OutputFileAscii::bytes_written() {
  // If the stream is open, then update the byte count. Otherwise, just
  // use the saved value.
  if ( stream_.is_open() ) {
    stream_.flush();
    byte_count_ = static_cast<int_fast64_t>( stream_.tellp() );
  }
  return byte_count_;
}

void marley::OutputFileAscii::write_event( HepMC3::GenEvent* event ) {

  if ( !event ) throw marley::Error( "Null pointer passed to"
    " OutputFileAscii::write_event()" );

  writer_->write_event( *event );
}

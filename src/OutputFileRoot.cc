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
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Data/GenRunInfoData.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/OutputFileRoot.hh"

marley::OutputFileRoot::OutputFileRoot( const marley::JSON& output_config )
  : marley::OutputFile( output_config )
{
  this->open();
}

marley::OutputFileRoot::~OutputFileRoot() {
  out_tfile_->cd();
  out_tree_->Write();
}

void marley::OutputFileRoot::open() {

  bool file_exists = check_if_file_exists( name_ );
  std::string open_mode_str( "update" );
  if ( mode_ == Mode::OVERWRITE ) open_mode_str = "recreate";

  // Open the ROOT file with the appropriate mode
  out_tfile_ = std::make_unique< TFile >( name_.c_str(),
    open_mode_str.c_str() );

  // Set up storage for event data to be written to the output
  temp_event_data_ptr_ = new HepMC3::GenEventData;
  temp_event_data_.reset( temp_event_data_ptr_ );

  // Check if the event TTree already exists in the file. If it does,
  // use it. Otherwise, create a new one and set up the "event" branch.
  out_tfile_->GetObject( "MARLEY_event_tree", out_tree_ );
  if ( !out_tree_ ) {
    out_tfile_->cd();
    out_tree_ = new TTree( "MARLEY_event_tree",
      "HepMC3-format MARLEY events" );
    out_tree_->Branch( "event", temp_event_data_ptr_ );
  }

  // Set the "event" branch address for the event data
  out_tree_->SetBranchAddress( "event", &temp_event_data_ptr_ );
}

bool marley::OutputFileRoot::resume(
  std::unique_ptr< marley::Generator >& gen, long& num_previous_events )
{
  if ( mode_ != Mode::RESUME ) {
    throw marley::Error( "Cannot call OutputFileRoot::resume() for an output"
      " mode other than \"resume\"" );
    return false;
  }

  // Retrieve the run information from the file
  HepMC3::GenRunInfoData* temp_run_info_data_ptr = nullptr;
  out_tfile_->GetObject( "MARLEY_run_info", temp_run_info_data_ptr );
  temp_run_info_data_.reset( temp_run_info_data_ptr );

  if ( !temp_run_info_data_ ) {
    throw marley::Error( "Failed to retrieve run information from the ROOT"
      " file \"" + name_ + "\"" );
    return false;
  }

  run_info_ = std::make_shared< HepMC3::GenRunInfo >();
  run_info_->read_data( *temp_run_info_data_ );

  // Get the generator JSON configuration
  auto config_attr = run_info_->attribute< HepMC3::StringAttribute >(
    "MARLEY.JSONconfig" );

  if ( !config_attr ) {
    throw marley::Error( "Failed to retrieve generator configuration from the"
      " ROOT file \"" + name_ + "\"" );
    return false;
  }

  std::string config_str = config_attr->value();

  // Initialize the generator using the JSON configuration
  marley::JSON config = marley::JSON::load( config_str );
  gen = this->restore_generator( config );

  // Retrieve the last event from the TTree
  num_previous_events = out_tree_->GetEntries();

  this->clear_event_data();
  out_tree_->GetEntry( num_previous_events - 1 );

  auto ev = std::make_shared< HepMC3::GenEvent >();
  ev->read_data( *temp_event_data_ );

  // Retrieve the generator state string from this event
  auto state_attr = ev->attribute< HepMC3::StringAttribute >(
    "MARLEY.GeneratorState" );

  if ( !state_attr ) {
    throw marley::Error( "Failed to retrieve generator state from the ROOT"
      " file \"" + name_ + "\"" );
    return false;
  }

  std::string state_str = state_attr->value();

  // Restore the random number generator state of the generator
  gen->seed_using_state_string( state_str );

  // Use the retrieved run information for making new events
  gen->set_run_info( run_info_ );

  return true;
}

void marley::OutputFileRoot::write_event( HepMC3::GenEvent* event ) {

  if ( !event ) throw marley::Error( "Null pointer passed to"
    " OutputFileRoot::write_event()" );

  // If the run information hasn't been set, then initialize it and save it
  auto ev_run_info = event->run_info();
  if ( !run_info_ && ev_run_info ) {
    run_info_ = ev_run_info;

    temp_run_info_data_ = std::make_unique< HepMC3::GenRunInfoData >();
    run_info_->write_data( *temp_run_info_data_ );

    out_tfile_->cd();
    out_tfile_->WriteObject( temp_run_info_data_.get(),
      "MARLEY_run_info", "WriteDelete" );
  }
  else if ( run_info_ != ev_run_info ) {
    throw marley::Error( "Unexpected change in MARLEY run information" );
  }

  // Clear out any pre-existing contents from the temporary event storage
  this->clear_event_data();
  event->write_data( *temp_event_data_ );
  out_tree_->Fill();
}

void marley::OutputFileRoot::clear_event_data() {
  temp_event_data_->particles.clear();
  temp_event_data_->vertices.clear();
  temp_event_data_->links1.clear();
  temp_event_data_->links2.clear();
  temp_event_data_->attribute_id.clear();
  temp_event_data_->attribute_name.clear();
  temp_event_data_->attribute_string.clear();
}

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
#include <string>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenParticle.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
// #include "TDirectory.h"
// #include "TVersionCheck.h"
// #include "TClass.h"
// #include "TObject.h"
// #include "TStorage.h"
// #include "TDataType.h"

// MARLEY includes
#include "marley/OutputFilePlainRoot.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"


#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#endif

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

marley::OutputFilePlainRoot::OutputFilePlainRoot( const marley::JSON& output_config )
  : marley::OutputFile( output_config )
{
  // We still need to initialise the HepMC3 objects, as our strategy is to produce
  // the plain ROOT file by applying marsum logic to the HepMC3 format.
  // this->open();
  reader_ = std::make_shared< HepMC3::ReaderAscii >( stream_ );
  writer_ = std::make_shared< HepMC3::WriterAscii >( stream_ );

  // Ensure that all floating-point values are output with full precision
  writer_->set_precision( std::numeric_limits<double>::max_digits10 );

  // Create the TFile and TTree pointers, and define the branches
  /// @todo Maybe you should move this to the header? Idk
  out_tfile_ = new TFile( name_.c_str(), "recreate" );
  out_tree_ = new TTree( "mst", "MARLEY summary tree" );

  // projectile branches
  out_tree_->Branch( "pdgv", &pdgv_, "pdgv/I" );
  out_tree_->Branch( "Ev", &Ev_, "Ev/D" );
  out_tree_->Branch( "KEv", &KEv_, "KEv/D" );
  out_tree_->Branch( "pxv", &pxv_, "pxv/D" );
  out_tree_->Branch( "pyv", &pyv_, "pyv/D" );
  out_tree_->Branch( "pzv", &pzv_, "pzv/D" );

  // target branches
  out_tree_->Branch( "pdgt", &pdgt_, "pdgt/I" );
  out_tree_->Branch( "Mt", &Mt_, "Mt/D" );

  // ejectile branches
  out_tree_->Branch( "pdgl", &pdgl_, "pdgl/I" );
  out_tree_->Branch( "El", &El_, "El/D" );
  out_tree_->Branch( "KEl", &KEl_, "KEl/D" );
  out_tree_->Branch( "pxl", &pxl_, "pxl/D" );
  out_tree_->Branch( "pyl", &pyl_, "pyl/D" );
  out_tree_->Branch( "pzl", &pzl_, "pzl/D" );

  // residue branches
  out_tree_->Branch( "pdgr", &pdgr_, "pdgr/I" );
  out_tree_->Branch( "Er", &Er_, "Er/D" );
  out_tree_->Branch( "KEr", &KEr_, "KEr/D" );
  out_tree_->Branch( "pxr", &pxr_, "pxr/D" );
  out_tree_->Branch( "pyr", &pyr_, "pyr/D" );
  out_tree_->Branch( "pzr", &pzr_, "pzr/D" );

  // Nuclear excitation energy branch
  out_tree_->Branch( "Ex", &Ex_, "Ex/D" );

  // Spin and parity branches
  out_tree_->Branch( "twoJ", &twoJ_, "twoJ/I" );
  out_tree_->Branch( "parity", &par_, "parity/I" );

  // De-excitation products (final-state particles other than the
  // ejectile and ground-state residue)
  out_tree_->Branch( "np", &np_, "np/I" );
  out_tree_->Branch( "pdgp", PDGs_.data(), "pdgp[np]/I" );
  out_tree_->Branch( "Ep",  Es_.data(), "Ep[np]/D" );
  out_tree_->Branch( "KEp", KEs_.data(), "KEp[np]/D" );
  out_tree_->Branch( "pxp", pXs_.data(), "pxp[np]/D" );
  out_tree_->Branch( "pyp", pYs_.data(), "pyp[np]/D" );
  out_tree_->Branch( "pzp", pZs_.data(), "pzp[np]/D" );

  // Flux-averaged total cross section
  out_tree_->Branch( "xsec", &flux_avg_tot_xsec_, "xsec/D" );

}

/// @todo check if this function is neccessary and if it is, implement it
void marley::OutputFilePlainRoot::open() {

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
      " OutputFilePlainRoot::open()" );
  }

  stream_.open( name_, open_mode_flag );

}

/// @todo Implement this function
bool marley::OutputFilePlainRoot::resume( std::unique_ptr<marley::Generator>& gen,
  long& num_previous_events )
{
  return false;
}

/// @todo Implement this function  
int_fast64_t marley::OutputFilePlainRoot::bytes_written() {
  // If the stream is open, then update the byte count. Otherwise, just
  // use the saved value.
  return 0;
}

void marley::OutputFilePlainRoot::write_event( HepMC3::GenEvent* ev ) {

  if ( !ev ) throw marley::Error( "Null pointer passed to"
    " OutputFilePlainRoot::write_event()" );

  //writer_->write_event( *event );

  // Instead of writing the event directly into ASCII, we "convert" it and write
  // it in plain ROOT format.
  PDGs_.clear();
  Es_.clear();
  KEs_.clear();
  pXs_.clear();
  pYs_.clear();
  pZs_.clear();

  auto projectile = marley_hepmc3::get_projectile( *ev );
  pdgv_ = projectile->pid();

  double mv = projectile->generated_mass();
  const HepMC3::FourVector& p4v = projectile->momentum();

  Ev_ = p4v.e();
  KEv_ = std::max( 0., Ev_ - mv );
  pxv_ = p4v.px();
  pyv_ = p4v.py();
  pzv_ = p4v.pz();

  auto target = marley_hepmc3::get_target( *ev );
  pdgt_ = target->pid();
  Mt_ = target->generated_mass();

  auto ejectile = marley_hepmc3::get_ejectile( *ev );
  pdgl_ = ejectile->pid();

  // Object ID number for the ejectile in the event (*not* the same as the
  // PDG code)
  int ej_id = ejectile->id();

  double ml = ejectile->generated_mass();
  const HepMC3::FourVector& p4l = ejectile->momentum();

  El_ = p4l.e();
  KEl_ = std::max( 0., El_ - ml );
  pxl_ = p4l.px();
  pyl_ = p4l.py();
  pzl_ = p4l.pz();

  auto residue = marley_hepmc3::get_residue( *ev );
  pdgr_ = residue->pid();

  double mr = residue->generated_mass();
  const HepMC3::FourVector& p4r = residue->momentum();

  Er_ = p4r.e();
  KEr_ = std::max( 0., Er_ - mr );
  pxr_ = p4r.px();
  pyr_ = p4r.py();
  pzr_ = p4r.pz();

  // TODO: add error handling here for missing attributes
  auto Ex_attr = residue->attribute< HepMC3::DoubleAttribute >( "Ex" );
  Ex_ = Ex_attr->value();

  auto twoJ_attr = residue->attribute< HepMC3::IntAttribute >( "twoJ" );
  twoJ_ = twoJ_attr->value();

  auto parity_attr = residue->attribute< HepMC3::IntAttribute >( "parity" );
  par_ = parity_attr->value();

  /// @todo Implement this
  flux_avg_tot_xsec_ = 0.; 

  np_ = 0;
  const auto& particles = ev->particles();
  for ( const auto& p : particles ) {
    // Only save information about final-state particles in the array of
    // nuclear de-excitation products
    if ( p->status() != marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS ) {
      continue;
    }

    // Skip the ejectile since its information is already saved elsewhere in
    // the output branches. We use the unique object ID number in the event
    // (not the same as the PDG code) to check for this
    if ( p->id() == ej_id ) continue;

    ++np_;

    PDGs_.push_back( p->pid() );

    double mp = p->generated_mass();
    const HepMC3::FourVector& p4p = p->momentum();

    double Ep = p4p.e();
    Es_.push_back( Ep );

    double KEp = std::max( 0., Ep - mp );
    KEs_.push_back( KEp );

    pXs_.push_back( p4p.px() );
    pYs_.push_back( p4p.py() );
    pZs_.push_back( p4p.pz() );
  }

  // Update the branch addresses (manipulating the vectors may have
  // invalidated them)
  out_tree_->SetBranchAddress( "pdgp", PDGs_.data() );
  out_tree_->SetBranchAddress( "Ep",  Es_.data() );
  out_tree_->SetBranchAddress( "KEp", KEs_.data() );
  out_tree_->SetBranchAddress( "pxp", pXs_.data() );
  out_tree_->SetBranchAddress( "pyp", pYs_.data() );
  out_tree_->SetBranchAddress( "pzp", pZs_.data() );

  out_tree_->Fill();
}
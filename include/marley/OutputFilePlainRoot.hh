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

// Standard library includes
#include <memory>

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// MARLEY includes
#include "marley/OutputFile.hh"

namespace HepMC3 {
  class ReaderAscii;
  class WriterAscii;
}

namespace marley {

  /// @brief Class that represents an output file in in the "plain" ROOT format
  class OutputFilePlainRoot : public OutputFile {

    public:

      OutputFilePlainRoot( const JSON& output_config );

      virtual ~OutputFilePlainRoot();

      virtual bool resume( std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events ) override;

      virtual int_fast64_t bytes_written() override;

      virtual void write_event( HepMC3::GenEvent* ev ) override;

    protected:

      /// Basic checks on the file before initializing
      virtual void open();

      /// Storage for the number of bytes written to disk
      int_fast64_t byte_count_ = 0;

      // TFile and TTree pointer objects for ROOT output
      TFile* out_tfile_;
      TTree* out_tree_;

      // Storage for output TTree branch variables
      double Ex_; // nuclear excitation energy
      int twoJ_; // two times the residue spin immediately after the two-two
                // scattering reaction
      int par_; // integer representation of the intrinsic parity of the
              // residue immediately following the two-two reaction
      double flux_avg_tot_xsec_; // flux-averaged total cross section
      double Ev_, KEv_, pxv_, pyv_, pzv_; // projectile
      double Mt_; // target mass
      double El_, KEl_, pxl_, pyl_, pzl_; // ejectile
      double Er_, KEr_, pxr_, pyr_, pzr_; // residue (after de-excitations)
      int pdgv_, pdgt_, pdgl_, pdgr_; // PDG codes
      int np_; // number of de-excitation products (final-state particles other
              // than the ejectile and residue)

      // Inpormation about each of the other final-state particles
      std::vector<int> PDGs_;
      std::vector<double> Es_, KEs_, pXs_, pYs_, pZs_;
  };

}

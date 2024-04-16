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

// standard library includes
#include <fstream>
#include <memory>
#include <string>

namespace marley {

  class Generator;
  class JSON;

  /// @brief Abstract base class for objects that deliver output to a file
  /// opened by the marley command-line executable
  class OutputFile {

    protected:

      // The input JSON object typically appears under the "output" key in the
      // "executable_settings" section of a MARLEY job configuration file
      OutputFile( const JSON& output_config );

    public:

      // Factory method that constructs an appropriate derived object given the
      // input JSON configuration
      static std::shared_ptr< OutputFile > make_OutputFile(
        const JSON& output_config );

      virtual ~OutputFile() = default;

      const std::string& name() const { return name_; }

      /// @brief Load a marley::Generator object whose configuration and state
      /// were saved to the metadata in the output file
      /// @param[out] gen A std::unique_ptr that will be filled with a restored
      /// Generator object constructed using the saved metadata
      /// @param[out] num_previous_events The number of events previously saved
      /// to the output file
      /// @return true if the Generator is successfully restored using the
      /// file's saved metadata, or false otherwise.
      virtual bool resume( std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events ) = 0;

      /// @brief The number of bytes that have been written during the
      /// current MARLEY session to this file
      virtual int_fast64_t bytes_written() = 0;

      /// @brief Write a new HepMC3::GenEvent to this output file
      /// @details If a nullptr is passed to this function, then
      /// a marley::Error will be thrown
      virtual void write_event( HepMC3::GenEvent* event ) = 0;

      bool mode_is_resume() const { return mode_ == Mode::RESUME; }

      // Add more formats to the enum class as needed. This should include
      // every event format that MARLEY knows how to write. The "ASCII" format
      // is MARLEY's native format for textual input and output of
      // HepMC3::GenEvent objects.
      enum class Format { ASCII, ROOT, PLAIN_ROOT };

    protected:

      /// @brief Checks that a given file exists (and is readable)
      inline bool check_if_file_exists( const std::string& filename ) {
        std::ifstream test_stream( filename );
        return test_stream.good();
      }

      /// @brief Helper function for resume() that instantiates a
      /// marley::Generator object given the previous JSON configuration object
      std::unique_ptr< marley::Generator > restore_generator(
        const marley::JSON& config );

      // Modes to use when writing output to files that are not initially empty
      // OVERWRITE = removes any previous contents of the file, then writes new
      //   events in the requested format.
      // RESUME = uses saved metadata in a file to restore a previous MARLEY
      //   configuration (including the random number generator state), then
      //   appends new events after those currently saved in the file.
      enum class Mode { OVERWRITE, RESUME };

      std::string name_; ///< Name of the file to receive output
      Format format_; ///< Format to use when writing events to the file
      Mode mode_;

      // Whether to force the current output mode without prompting the user.
      // Currently only used when setting the OVERWRITE output mode.
      bool force_;
  };

}

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

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/OutputFile.hh"
#include "marley/OutputFileAscii.hh"
#include "marley/OutputFilePlainRoot.hh"

#ifdef USE_ROOT
  #include "marley/OutputFileRoot.hh"
#endif

marley::OutputFile::OutputFile( const marley::JSON& config ) {

  if ( !config.has_key("file") ) throw marley::Error( "Missing file"
    " name for an output file specification in the configuration file." );

  name_ = config.at("file").to_string();

  force_ = false; // default behavior is to prompt before overwriting
  if ( config.has_key("force") ) force_ = config.at( "force" ).to_bool();

  std::string mode_str( "overwrite" ); // default mode is "overwrite"
  if ( config.has_key("mode") ) mode_str = config.at( "mode" ).to_string();

  if ( mode_str == "overwrite" ) mode_ = Mode::OVERWRITE;
  else if ( mode_str == "resume" ) mode_ = Mode::RESUME;
  else throw marley::Error( "Invalid output mode \"" + mode_str
    + "\" given in an output file specification" );
}

std::unique_ptr< marley::Generator > marley::OutputFile::restore_generator(
  const marley::JSON& config )
{
  marley::JSONConfig jc( config );
  auto gen = std::make_unique< marley::Generator >( jc.create_generator() );
  return gen;
}

// Factory method that constructs an appropriate derived object given the
// input JSON configuration
std::shared_ptr< marley::OutputFile > marley::OutputFile::make_OutputFile(
  const JSON& output_config )
{
  if ( !output_config.has_key("format") ) throw marley::Error( "Missing format"
    " for an output file specification in the configuration file." );

  std::string format = output_config.at( "format" ).to_string();

  if ( format == "ascii" ) {
    return std::make_shared< marley::OutputFileAscii >( output_config );
  }
  #ifdef USE_ROOT
  else if ( format == "root" ) {
    return std::make_shared< marley::OutputFileRoot >( output_config );
  }
  else if ( format == "plain-root" ) {
    return std::make_shared< marley::OutputFilePlainRoot >( output_config );
  }
  #endif
  else throw marley::Error( "Invalid output file format \"" + format
    + "\" given in an output file specification" );
}

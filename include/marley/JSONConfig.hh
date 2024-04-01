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
#include <string>

// MARLEY includes
#include "marley/CoulombCorrector.hh"
#include "marley/FormFactor.hh"
#include "marley/Generator.hh"
#include "marley/InterpolationGrid.hh"
#include "marley/JSON.hh"

namespace marley {

  class JSONConfig {

    public:

      using InterpMethod = marley::InterpolationGrid<double>
        ::InterpolationMethod;

      explicit JSONConfig( const marley::JSON& object );
      explicit JSONConfig( const std::string& json_filename );

      marley::Generator create_generator() const;

      void prepare_direction( marley::Generator& gen ) const;
      void prepare_neutrino_source( marley::Generator& gen ) const;
      void prepare_reactions( marley::Generator& gen,
        marley::CoulombCorrector::CoulombMode coulomb_mode,
        marley::FormFactor::FFScalingMode ff_scaling_mode,
        bool superallowed ) const;
      void prepare_structure( marley::Generator& gen ) const;
      void prepare_target( marley::Generator& gen ) const;

      void update_logger_settings() const;

      InterpMethod get_interpolation_method( const std::string& rule ) const;
      int neutrino_pdg( const std::string& nu ) const;

      /// @brief Helper function used to define ROOT-based neutrino source
      /// types
      /// @detail This function is a no-op when MARLEY is built without
      /// ROOT support
      bool process_extra_source_types( const std::string& type,
        const marley::JSON& source_spec, int pdg_code,
        std::unique_ptr< marley::NeutrinoSource >& source ) const;

      inline const marley::JSON& get_json() const;
      inline void set_json( const marley::JSON& json );

      static void handle_json_error( const std::string& name,
        const marley::JSON& json );

    protected:

      /// @brief Helper function for loading strings from the JSON
      /// configuration
      std::string source_get( const char* name,
        const marley::JSON& source_spec, const char* description,
        const char* default_str ) const;

      /// @brief JSON object describing this configuration
      marley::JSON json_;
  };

  inline const marley::JSON& marley::JSONConfig::get_json() const
    { return json_; }

  inline void marley::JSONConfig::set_json( const marley::JSON& json )
    { json_ = json; }

}

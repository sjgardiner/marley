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

// standard library includes
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/JSONConfig.hh"
#include "marley/NeutrinoSource.hh"
#include "marley/NuclearReaction.hh"
#include "marley/Logger.hh"
#include "marley/StructureDatabase.hh"
#include "marley/FormFactor.hh"

#ifdef USE_ROOT
#include "marley/marley_root.hh"
#endif

using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;
using ProcType = marley::Reaction::ProcessType;
using CMode = marley::CoulombCorrector::CoulombMode;
using CRPADiscreteMode = marley::Generator::CRPADiscreteMode;
using FFScalingMode = marley::FormFactor::FFScalingMode;

// anonymous namespace for helper functions, etc.
namespace {

  void source_check_positive( double x, const char* description,
    const char* source_type )
  {
    if ( x <= 0. ) throw marley::Error( std::string("Non-positive ")
      + description + " value defined for a " + source_type
      + " neutrino source" );
  }

  void source_check_nonnegative(double x, const char* description,
    const char* source_type)
  {
    if ( x < 0. ) throw marley::Error( std::string("Negative ") + description
      + " value defined for a " + source_type + " neutrino source" );
  }

  double source_get_double( const char* name, const marley::JSON& source_spec,
    const char* description )
  {
    if ( !source_spec.has_key(name) ) throw marley::Error(
      std::string("Missing source.") + name + " key for " + description
      + " source" );
    bool ok;
    double result = source_spec.at( name ).to_double( ok );
    if ( !ok ) throw marley::Error( std::string("Invalid value given for"
      " source.") + name + " key for " + description + " source" );
    return result;
  }

  std::vector<double> get_vector( const char* name, const marley::JSON& spec,
    const char* description )
  {
    if ( !spec.has_key(name) ) throw marley::Error( std::string("The")
      + " specification for a " + description + " source should include"
      + " a " + name + " key." );

    const marley::JSON& vec = spec.at( name );

    if ( !vec.is_array() ) throw marley::Error( std::string("The")
      + " value given for the " + name + " key for a " + description
      + " source should be an array." );

    std::vector<double> result;

    auto elements = vec.array_range();
    if ( elements.begin() != elements.end() ) {
      bool ok;
      for ( const auto& el : elements ) {
        double dub = el.to_double( ok );
        if ( ok ) result.push_back( dub );
        else throw marley::Error( "Invalid array entry '"
          + el.dump_string() + "' given for the " + name + " key for a "
          + description + " source specification." );
      }
    }

    return result;
  }

}

marley::JSONConfig::JSONConfig( const marley::JSON& json ) : json_( json )
{
  update_logger_settings();
}

marley::JSONConfig::JSONConfig( const std::string& json_filename )
{
  // First update the Logger settings so we can have default logging
  // up and running when we parse the JSON configuration
  update_logger_settings();

  // Parse the config file
  json_ = marley::JSON::load_file( json_filename );

  // Now update the logger settings based on the file contents
  update_logger_settings();
}

int marley::JSONConfig::neutrino_pdg( const std::string& nu ) const {

  int pdg = 0;

  bool bad = false;

  // Matches integers
  static const std::regex rx_int = std::regex( "[-+]?[0-9]+" );
  if ( std::regex_match(nu, rx_int) ) {
    pdg = std::stoi( nu );
    if ( !marley::NeutrinoSource::pdg_is_allowed(pdg) ) bad = true;
  }
  else if ( !marley_utils::string_to_neutrino_pdg(nu, pdg) ) {
    bad = true;
  }

  if ( bad ) throw marley::Error( "Invalid neutrino type specification '"
    + nu + "' given for the MARLEY neutrino source." );

  return pdg;
}

marley::Generator marley::JSONConfig::create_generator() const
{
  uint_fast64_t seed;
  if ( json_.has_key("seed") ) {
    bool ok;
    seed = static_cast< uint_fast64_t >( json_.at("seed").to_long(ok) );
    if ( !ok ) handle_json_error( "seed", json_.at("seed") );
  }
  else seed = std::chrono::system_clock::now().time_since_epoch().count();

  // Start with a default-constructed Generator seeded with either a
  // user-supplied seed or the current number of seconds since the Unix epoch.
  marley::Generator gen( seed );

  // Set the method to use for computing Coulomb corrections in all reactions.
  // If the user gave an explicit setting for this, use that.
  // Otherwise, interpolate between the Fermi function and the modified
  // effective momentum approximation.
  CMode coulomb_mode = CMode::FERMI_AND_MEMA; // Default method
  if ( json_.has_key("coulomb_mode") ) {
    const auto& cmode = json_.at( "coulomb_mode" );
    if ( !cmode.is_string() ) throw marley::Error("Invalid Coulomb mode"
      " specification " + cmode.dump_string() );
    std::string my_mode = cmode.to_string();
    coulomb_mode = marley::CoulombCorrector
      ::coulomb_mode_from_string( my_mode );
  }

  // Set the way the form factors depend on Q^2 in all reactions.
  FFScalingMode ff_scaling_mode = FFScalingMode::FLAT; // Default mode
  if ( json_.has_key("ff_scaling_mode") ) {
    const auto& ffmode = json_.at( "ff_scaling_mode" );
    if ( !ffmode.is_string() ) throw marley::Error("Invalid form factor"
      " scaling mode specification " + ffmode.dump_string() );
    std::string my_mode = ffmode.to_string();
    ff_scaling_mode = marley::FormFactor
      ::ff_scaling_mode_from_string( my_mode );
  }
  // Inform the user about the set form factor scaling mode
  MARLEY_LOG_INFO() << "Configured form factor scaling mode: " << marley::FormFactor::string_from_ff_scaling_mode( ff_scaling_mode );


  // Turn off calls to Generator::normalize_E_pdf() until we
  // have set up all the needed pieces
  gen.dont_normalize_E_pdf_ = true;

  // Use the JSON settings to update the generator's parameters
  prepare_direction( gen );
  prepare_structure( gen );
  prepare_neutrino_source( gen );
  prepare_reactions( gen, coulomb_mode, ff_scaling_mode );
  prepare_target( gen );

  // If the user has disabled nuclear de-excitations, then set the
  // flag appropriately.
  if ( json_.has_key("do_deexcitations") ) {
    const auto& do_deex = json_.at( "do_deexcitations" );
    if ( do_deex.is_bool() ) {
      bool deexcite_or_not = do_deex.to_bool();
      gen.set_do_deexcitations( deexcite_or_not );
      if ( !deexcite_or_not ) {
        MARLEY_LOG_INFO() << "Nuclear de-excitations will not be simulated";
      }
    }
  }

  // Set the flag to specify how to deal with CRPA strength 
  // leaking below the unbound threshold
  CRPADiscreteMode crpa_discrete_mode = CRPADiscreteMode::IGNORE; // Ignore by default
  if ( json_.has_key("crpa_discrete_mode") ) {
    const auto& crpamode = json_.at( "crpa_discrete_mode" );
    if ( !crpamode.is_string() ) throw marley::Error("Invalid CRPA discrete mode"
      " specification " + crpamode.dump_string() );
    std::string my_mode = crpamode.to_string();
    crpa_discrete_mode = marley::Generator
      ::crpa_discrete_mode_from_string( my_mode );
    
    gen.set_crpa_discrete_mode( crpa_discrete_mode );
  }
  // Inform the user about the set CRPA discrete mode
  MARLEY_LOG_INFO() << "Configured CRPA discrete mode: " << marley::Generator::string_from_crpa_discrete_mode( crpa_discrete_mode );

  // Skip the rest of initialization if we've disabled all reactions.
  // This can be used to partially initialize the Generator in unusual
  // situations.
  if ( json_.has_key("reactions") ) {
    const auto& reactions = json_.at( "reactions" );
    if ( reactions.is_null() ) {
      MARLEY_LOG_INFO() << "Null reactions array detected."
        << " Initialization of reactions will be skipped.";
      return gen;
    }
  }

  // Iterate through all configured nuclear reactions. Set a flag indicating
  // whether at least one of them was a CC reaction. We will only bother to
  // print the logging message below if one such reaction was found.
  bool found_cc = false;
  for ( auto& react : gen.reactions_ ) {

    ProcType pt = react->process_type();
    if ( pt == ProcType::NeutrinoCC || pt == ProcType::AntiNeutrinoCC ) {
      found_cc = true;
    }
  }

  // If a CC reaction is configured, then print a logging message indicating
  // which Coulomb correction method is active. Otherwise, don't bother.
  if ( found_cc ) {
    std::string cmode_str = marley::CoulombCorrector
      ::string_from_coulomb_mode( coulomb_mode );
    MARLEY_LOG_INFO() << "Configured Coulomb correction method: " << cmode_str;
  }

  // Now that the reactions and source are both prepared, check that a neutrino
  // from the source can interact via at least one of the enabled reactions
  bool found_matching_pdg = false;
  int source_pdg = gen.get_source().get_pid();
  for ( const auto& react : gen.get_reactions() ) {
    if ( source_pdg == react->pdg_a() ) found_matching_pdg = true;
  }
  // If neutrinos from the source can never interact, then complain about it
  if ( !found_matching_pdg ) throw marley::Error( "The neutrino source"
    " produces " + marley_utils::get_particle_symbol(source_pdg)
    + ", which cannot participate in any of the configured reactions." );

  // Before returning the newly-created Generator object, print logging
  // messages describing the reactions that are active.
  MARLEY_LOG_INFO() << "Generator configuration complete. Active reactions:";
  for ( const auto& r : gen.get_reactions() ) {

    const marley::TargetAtom ta = r->atomic_target();
    double atom_frac = gen.get_target().atom_fraction( ta );

    if ( r->pdg_a() == source_pdg && atom_frac > 0. ) {

      std::string proc_type_str;
      if ( r->process_type() == ProcType::NeutrinoCC
        || r->process_type() == ProcType::AntiNeutrinoCC )
      {
        proc_type_str = "CC";
      }
      else if ( r->process_type() == ProcType::NC )
      {
        proc_type_str = "NC";
      }
      else if ( r->process_type() == ProcType::NuElectronElastic )
      {
        proc_type_str = "ES on " + ta.to_string();
      }
      else throw marley::Error( "Unrecognized process type encountered in"
        " marley::JSONConfig::prepare_reactions()" );

      // Show the threshold in red if it's above the maximum energy
      // produced by the source
      std::ostringstream temp_oss;
      double threshold_KE = r->threshold_kinetic_energy();
      bool no_flux = ( threshold_KE > gen.get_source().get_Emax() );
      if ( no_flux ) temp_oss << "\u001b[31m";
      temp_oss << threshold_KE << " MeV";
      if ( no_flux ) temp_oss << "\u001b[30m";
      temp_oss << ')';

      MARLEY_LOG_INFO() << "  " << proc_type_str << ": "
        << r->get_description() << " (KE @ threshold: "
        << temp_oss.str();
    }
  }

  // We've prepared the Generator and checked that at least one cross section
  // should be non-vanishing. Now actually integrate the cross section to check
  // that there is flux above threshold (and normalize the energy PDF at the
  // same time). An exception will be thrown if no neutrinos can interact.
  gen.dont_normalize_E_pdf_ = false;
  gen.normalize_E_pdf();

  // Save a copy of the JSON settings used to configure the generator
  gen.set_json_config( json_ );

  // Now we're all ready to go. Log the flux-averaged total cross section
  // value before returning the fully-configured generator.
  double avg_tot_xs = gen.flux_averaged_total_xs(); // MeV^(-2)
  MARLEY_LOG_INFO() << "Flux-averaged total cross section per atom: "
    << marley_utils::hbar_c2 * avg_tot_xs * marley_utils::fm2_to_minus40_cm2
    << " * 10^(-40) cm^2";

  return gen;
}

//------------------------------------------------------------------------------
void marley::JSONConfig::prepare_direction( marley::Generator& gen ) const {
  // Get the incident neutrino direction if the user has specified one
  if ( json_.has_key("direction") ) {

    const marley::JSON& direction = json_.at("direction");
    bool ok;

    // The usual use case is for the user to specify the components
    // of a direction 3-vector using a JSON object
    if ( direction.is_object() ) {

      std::array<double, 3> dir_vec = gen.neutrino_direction();

      if ( direction.has_key("x") ) {
        dir_vec.at(0) = direction.at( "x" ).to_double( ok );
        if ( !ok ) handle_json_error( "direction.x", direction.at("x") );
      }

      if ( direction.has_key("y") ) {
        dir_vec.at(1) = direction.at( "y" ).to_double( ok );
        if ( !ok ) handle_json_error( "direction.y", direction.at("y") );
      }

      if ( direction.has_key("z") ) {
        dir_vec.at(2) = direction.at( "z" ).to_double( ok );
        if ( !ok ) handle_json_error( "direction.z", direction.at("z") );
      }

      gen.set_neutrino_direction( dir_vec );
    }

    // The user may also request sampling of an isotropic projectile direction
    // for every event by associating the string value "isotropic" with the
    // direction key in the job configuration file
    else if ( direction.is_string() && direction.to_string() == "isotropic" ) {
      gen.get_rotator().set_randomize_directions( true );

      MARLEY_LOG_INFO() << "Projectile directions will be sampled"
        << " isotropically";
    }
    else {
      throw marley::Error( "Unrecognized value "
        + direction.dump_string() + " given for the job configuration file"
        " key \"direction\"" );
    }
  }
}

void marley::JSONConfig::prepare_reactions( marley::Generator& gen,
  CMode coulomb_mode, FFScalingMode ff_scaling_mode ) const
{
  const auto& fm = marley::FileManager::Instance();

  if ( json_.has_key("reactions") ) {

    const marley::JSON& rs = json_.at( "reactions" );

    // If the reactions key has a null value, skip trying
    // to load any reaction data. This can be used in unusual
    // situations when we don't actually want to simulate any
    // reactions.
    if ( rs.is_null() ) return;

    if ( rs.is_array() ) {

      auto reactions = rs.array_range();
      if ( reactions.begin() != reactions.end() ) {

        // Create a temporary vector to cache the (TargetAtom, ProcessType)
        // pairs for which reactions have already been loaded. Complain if
        // there is duplication.
        std::vector< std::pair<marley::TargetAtom, ProcType> >
          loaded_proc_types;

        for ( const auto& r : reactions ) {

          std::string filename = r.to_string();

          // Find the reaction data file using the MARLEY search path
          std::string full_file_name = fm.find_file( filename );
          if ( full_file_name.empty() ) {
            throw marley::Error( "Could not locate the reaction data file "
              + filename + ". Please check that the file name is spelled"
              " correctly and that the file is in a folder"
              " on the MARLEY search path." );
          }

          auto reacts = marley::Reaction::load_from_file(
            full_file_name, gen.get_structure_db(), coulomb_mode, ff_scaling_mode );

          if ( reacts.empty() ) throw marley::Error( "Failed to load"
            " any reactions from the file " + full_file_name + ". Please"
            " check that it is readable and conforms to the correct input"
            " format." );

          // All of the Reaction objects loaded from a single file will have
          // the same process type and atomic target, so just save this
          // information from the first one
          auto temp_atom = reacts.front()->atomic_target();
          auto temp_pt = reacts.front()->process_type();
          std::pair< marley::TargetAtom, ProcType >
            temp_pair( temp_atom, temp_pt );

          // If we have a duplicate, warn the user that we'll ignore it
          auto begin = loaded_proc_types.cbegin();
          auto end = loaded_proc_types.cend();
          if ( std::find(begin, end, temp_pair) != end ) {
            MARLEY_LOG_WARNING() << "Reaction settings for the "
              << marley::Reaction::proc_type_to_string( temp_pt )
              << " process on " << temp_atom << " were already loaded."
              << " To avoid duplication, those in " << full_file_name
              << " will be ignored.";
            continue;
          }
          // Otherwise, save the process type for later checks of this kind
          else {
            MARLEY_LOG_INFO() << "Loaded "
              << marley::Reaction::proc_type_to_string( temp_pt )
              << " reaction data for " << temp_atom << " from "
              << full_file_name;
            loaded_proc_types.push_back( temp_pair );
          }

          // Transfer ownership of the new reactions to the generator
          for ( auto& rct : reacts ) gen.add_reaction( std::move(rct) );
        }

        return;
      }
      else {
        throw marley::Error( "At least one reaction matrix data file must be"
          " specified using the \"reactions\" parameter" );
      }
    }

    handle_json_error( "reactions", rs );
  }

  throw marley::Error( "Missing \"reactions\" key in the MARLEY configuration"
    " file." );
}

void marley::JSONConfig::prepare_structure( marley::Generator& gen ) const
{
  auto& sdb = gen.get_structure_db();

  // Check for a custom configuration of nuclear optical model parameters.
  // If the user asked for one, use it instead of the default settings.
  const std::string om_key = "opt_mod";
  if ( json_.has_key(om_key) ) {
    const marley::JSON om_config = json_.at( om_key );
    MARLEY_LOG_INFO() << "Loading custom optical model configuration";
    MARLEY_LOG_DEBUG() << om_config.dump_string();

    sdb.load_optical_model_params( &om_config );
  }

  // If the user specified a non-default value of either the
  // maximum orbital angular momentum or the maximum multipolarity
  // to consider when simulating decays to the continuum, set
  // the appropriate member variable of the StructureDatabase
  // object owned by the Generator

  std::string flmax_key( "fragment_lmax" );
  if ( json_.has_key(flmax_key) ) {
    bool ok;
    const marley::JSON& flmax_json = json_.at( flmax_key );
    int f_lmax = flmax_json.to_long( ok );
    if ( !ok ) handle_json_error( flmax_key.c_str(), flmax_json );

    if ( f_lmax < 0 ) throw marley::Error( "Negative value of "
      + flmax_key + " = " + std::to_string(f_lmax) + " encountered in"
      " marley::JSONConfig::prepare_structure()" );

    sdb.set_fragment_l_max( f_lmax );

    MARLEY_LOG_INFO() << "Orbital angular momentum cutoff for fragment"
      << " differential decay widths set to l_max = " << f_lmax;
  }

  // TODO: reduce code duplication here
  std::string glmax_key( "gamma_lmax" );
  if ( json_.has_key(glmax_key) ) {
    bool ok;
    const marley::JSON& glmax_json = json_.at( glmax_key );
    int g_lmax = glmax_json.to_long( ok );
    if ( !ok ) handle_json_error( glmax_key.c_str(), glmax_json );

    if ( g_lmax < 1 ) throw marley::Error( "Nonpositive value of "
      + glmax_key + " = " + std::to_string(g_lmax) + " encountered in"
      " marley::JSONConfig::prepare_structure()" );

    sdb.set_gamma_l_max( g_lmax );

    MARLEY_LOG_INFO() << "Multipolarity cutoff for gamma-ray"
      << " differential decay widths set to l_max = " << g_lmax;
  }

}

//------------------------------------------------------------------------------
InterpMethod marley::JSONConfig::get_interpolation_method(
  const std::string& rule ) const
{
  // Try using the ENDF-style numerical codes first
  static const std::regex rx_nonneg_int( "[0-9]+" );

  if ( std::regex_match(rule, rx_nonneg_int) ) {
    int endf_interp_code = std::stoi( rule );
    if ( endf_interp_code == 1 ) return InterpMethod::Constant;
    else if ( endf_interp_code == 2 ) return InterpMethod::LinearLinear;
    else if ( endf_interp_code == 3 ) return InterpMethod::LinearLog;
    else if ( endf_interp_code == 4 ) return InterpMethod::LogLinear;
    else if ( endf_interp_code == 5 ) return InterpMethod::LogLog;
  }

  // Interpolation rules may also be given as strings
  else if ( rule == "const" || rule == "constant" )
    return InterpMethod::Constant;
  else if ( rule == "lin" || rule == "linlin" )
    return InterpMethod::LinearLinear;
  else if ( rule == "log" || rule == "loglog" )
    return InterpMethod::LogLog;
  // linear in energy, logarithmic in probability density
  else if ( rule == "linlog" )
    return InterpMethod::LinearLog;
  // logarithmic in energy, linear in probability density
  else if ( rule == "loglin" )
    return InterpMethod::LogLinear;
  else throw marley::Error( "Invalid interpolation rule '" + rule
    + "' given in the neutrino source specification" );

  // We shouldn't ever end up here, but return something just in case
  return InterpMethod::Constant;
}

//------------------------------------------------------------------------------
void marley::JSONConfig::prepare_neutrino_source( marley::Generator& gen ) const
{
  // Check whether the user provided their own estimate of the source PDF
  // maximum value. If they did, adopt that before building the source.
  // This is useful when automatic searches for the maximum don't work.
  // The user can put in a value manually to get rejection sampling to work.
  if ( json_.has_key("energy_pdf_max") ) {
    bool ok;
    const marley::JSON& max_spec = json_.at( "energy_pdf_max" );
    double user_max = max_spec.to_double( ok );
    if ( !ok ) handle_json_error( "energy_pdf_max", max_spec );
    else gen.set_default_E_pdf_max( user_max );
  }

  // Check whether the JSON configuration includes a neutrino source
  // specification
  if ( !json_.has_key("source") ) return;
  const marley::JSON& source_spec = json_.at( "source" );

  // If the neutrino source key has a null value, just return without doing
  // anything else
  if ( source_spec.is_null() ) {
    MARLEY_LOG_INFO() << "Null source specification detected. Skipping"
      << " neutrino source configuration.";
    return;
  }

  // Complain if the user didn't specify a source type
  if ( !source_spec.has_key("type") ) {
    throw marley::Error( "Missing \"type\" key in neutrino source"
      " specification." );
    return;
  }

  // Get the neutrino source type
  bool ok;
  std::string type = source_spec.at( "type" ).to_string( ok );
  if ( !ok ) handle_json_error( "source.type", source_spec.at("type") );

  // Complain if the user didn't specify a neutrino type
  if ( !source_spec.has_key("neutrino") ) {
    throw marley::Error( "Missing \"neutrino\" key in neutrino source"
      " specification." );
    return;
  }
  // Get the neutrino type
  std::string nu = source_spec.at( "neutrino" ).to_string( ok );
  if ( !ok ) handle_json_error( "source.neutrino", source_spec.at("neutrino") );

  // Particle Data Group code for the neutrino type produced by this source
  int pdg = neutrino_pdg( nu );

  std::unique_ptr< marley::NeutrinoSource > source;

  if ( type == "mono" || type == "monoenergetic" ) {
    double energy = source_get_double( "energy", source_spec, "monoenergetic" );
    source_check_positive( energy, "energy", "monoenergetic" );
    source = std::make_unique< marley::MonoNeutrinoSource >( pdg, energy );
    MARLEY_LOG_INFO() << "Created monoenergetic "
      << marley_utils::get_particle_symbol( pdg ) << " source with"
      << " neutrino energy = " << energy << " MeV";
  }
  else if ( type == "dar" || type == "decay-at-rest" ) {
    source = std::make_unique< marley::DecayAtRestNeutrinoSource >( pdg );
     MARLEY_LOG_INFO() << "Created muon decay-at-rest "
       << marley_utils::get_particle_symbol( pdg ) << " source";
  }
  else if ( type == "fd" || type == "fermi-dirac" || type == "fermi_dirac" ) {
    double Emin = source_get_double( "Emin", source_spec, "Fermi-Dirac" );
    double Emax = source_get_double( "Emax", source_spec, "Fermi-Dirac" );
    double temp = source_get_double( "temperature", source_spec,
      "Fermi-Dirac" );

    double eta = 0.;
    if ( source_spec.has_key("eta") ) {
      eta = source_get_double( "eta", source_spec, "Fermi-Dirac" );
    }

    source_check_nonnegative( Emin, "Emin", "Fermi-Dirac" );
    source_check_positive( temp, "temperature", "Fermi-Dirac" );

    if ( Emax <= Emin ) throw marley::Error( "Emax <= Emin for a Fermi-Dirac"
      " neutrino source" );

    source = std::make_unique< marley::FermiDiracNeutrinoSource >( pdg, Emin,
      Emax, temp, eta );
    MARLEY_LOG_INFO() << "Created Fermi-Dirac "
      << marley_utils::get_particle_symbol( pdg ) << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << Emin << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << Emax << " MeV";
    MARLEY_LOG_INFO() << "  temperature = " << temp << " MeV";
    MARLEY_LOG_INFO() << "  eta = " << eta;
  }
  else if ( type == "bf" || type == "beta" || type == "beta-fit" ) {
    double Emin = source_get_double( "Emin", source_spec, "beta-fit" );
    double Emax = source_get_double( "Emax", source_spec, "beta-fit" );
    double Emean = source_get_double( "Emean", source_spec, "beta-fit" );

    double beta = 4.5;
    if ( source_spec.has_key("beta") ) {
      beta = source_get_double( "beta", source_spec, "beta-fit" );
    }

    source_check_nonnegative( Emin, "Emin", "beta-fit" );
    source_check_positive( Emean, "Emean", "beta-fit" );

    if ( Emax <= Emin ) throw marley::Error( "Emax <= Emin for a beta-fit"
      " neutrino source" );

    source = std::make_unique< marley::BetaFitNeutrinoSource >( pdg, Emin,
      Emax, Emean, beta );
    MARLEY_LOG_INFO() << "Created beta-fit "
      << marley_utils::get_particle_symbol( pdg ) << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << Emin << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << Emax << " MeV";
    MARLEY_LOG_INFO() << "  average energy = " << Emean << " MeV";
    MARLEY_LOG_INFO() << "  beta = " << beta;
  }
  else if ( type == "hist" || type == "histogram" ) {

    std::vector< double > Es = get_vector( "E_bin_lefts", source_spec,
      "histogram" );
    std::vector< double > weights = get_vector( "weights", source_spec,
      "histogram" );

    if ( Es.size() != weights.size() ) throw marley::Error( "The sizes of the"
      " arrays of energy bin left edges and weights given for a histogram"
      " neutrino source are unequal." );

    double Emax = source_get_double( "Emax", source_spec, "histogram" );
    source_check_positive( Emax, "Emax", "histogram" );

    // Add Emax to the grid
    Es.push_back( Emax );

    // Set the probability density at E = Emax to be zero (this ensures
    // that no energies outside of the histogram will be sampled)
    weights.push_back( 0. );

    // Convert from bin weights to probability densities by dividing by the
    // width of each bin
    int jmax = Es.size() - 1;
    for ( int j = 0; j < jmax; ++j ) {

      double width = Es.at( j + 1 ) - Es.at( j );
      if ( width <= 0 ) throw marley::Error( "Invalid bin width"
        + std::to_string(width) + " encountered when creating a histogram"
        " neutrino source" );

      weights.at( j ) /= width;
    }

    // Create the source
    source = std::make_unique< marley::GridNeutrinoSource >( Es, weights, pdg,
      InterpMethod::Constant );
    MARLEY_LOG_INFO() << "Created histogram "
      << marley_utils::get_particle_symbol( pdg ) << " source";
  }
  else if ( type == "grid" ) {
    std::vector< double > energies = get_vector( "energies", source_spec,
      "grid" );
    std::vector< double > PDs = get_vector( "prob_densities", source_spec,
      "grid" );
    std::string rule = source_get( "rule", source_spec, "grid", "linlin" );

    InterpMethod method = get_interpolation_method( rule );

    source = std::make_unique< marley::GridNeutrinoSource >( energies, PDs,
      pdg, method );
    MARLEY_LOG_INFO() << "Created grid "
      << marley_utils::get_particle_symbol( pdg ) << " source";
  }
  else if ( !process_extra_source_types(type, source_spec, pdg, source) ) {
    throw marley::Error( "Unrecognized MARLEY neutrino source type '"
      + type + "'" );
  }

  // If the user has specified whether to weight the incident neutrino spectrum
  // by the reaction cross section(s), then set the weight_flux_ flag in the
  // new Generator object accordingly
  if ( source_spec.has_key("weight_flux") ) {
    bool ok = false;
    bool should_we_weight = source_spec.at( "weight_flux" ).to_bool( ok );
    if ( !ok ) handle_json_error( "source.weight_flux",
      source_spec.at("weight_flux") );
    gen.set_weight_flux( should_we_weight );
  }

  // Load the generator with the new source object
  gen.set_source( std::move(source) );

}

void marley::JSONConfig::prepare_target( marley::Generator& gen ) const {

  // Temporary storage for the list of target atoms and their atom fractions
  // in the possibly-composite neutrino target
  std::vector< marley::TargetAtom > atoms;
  std::vector< double > atom_fractions;

  // In the absence of a user-specified target, create one automatically from
  // the configured reactions. Each unique target atom involved in at least
  // one reaction will be included with equal weight.
  if ( !json_.has_key("target") ) {

    // If there aren't any configured reactions, just return without
    // configuring the target at all. This only happens in unusual
    // situations.
    const auto& reactions = gen.get_reactions();
    if ( reactions.empty() ) return;

    // Otherwise, store the set of target atoms involved in at least
    // one reaction
    std::set< marley::TargetAtom > temp_atom_set;
    for ( const auto& react : reactions ) {
      temp_atom_set.insert( react->atomic_target() );
    }

    // Add each target atom with equal weight to the target configuration
    for ( const auto& atom : temp_atom_set ) {
      atoms.push_back( atom );
      // This will be renormalized appropriately by the Target object itself
      atom_fractions.push_back( 1. );
    }
  }
  else {
    // If the user has specified a target composition explicitly, parse the
    // JSON object used to define it.
    const auto& tgt_spec = json_.at( "target" );
    if ( !tgt_spec.is_object() ) throw marley::Error( "Invalid neutrino target"
      " specification " + tgt_spec.dump_string() );

    if ( !tgt_spec.has_key("nuclides") ) throw marley::Error( "Missing \""
      "nuclides\" key in the neutrino target specification "
      + tgt_spec.dump_string() );

    const auto& n_spec = tgt_spec.at( "nuclides" );

    if ( !n_spec.is_array() ) throw marley::Error( "Invalid \"nuclides\""
      " array given in the neutrino target specification "
      + tgt_spec.dump_string() );

    if ( !tgt_spec.has_key("atom_fractions") ) throw marley::Error( "Missing \""
      "atom_fractions\" key in the neutrino target specification "
      + tgt_spec.dump_string() );

    const auto& af_spec = tgt_spec.at( "atom_fractions" );

    if ( !af_spec.is_array() ) throw marley::Error( "Invalid"
      " \"atom_fractions\" array given in the neutrino target specification "
      + tgt_spec.dump_string() );

    // Check that the two arrays used to specify the target are of equal length
    int num_nuclides = n_spec.length();
    if ( num_nuclides != af_spec.length() ) throw marley::Error(
      "Arrays of unequal length specified for the \"nuclides\" and \"atom"
      "_fractions\" keys in the neutrino target specification "
      + tgt_spec.dump_string() );

    // Check that at least one target atom is listed
    if ( num_nuclides < 1 ) throw marley::Error( "At least one target nuclide"
      " must be included in the neutrino target specification" );

    // Loop over each of the target atoms. Parse their information and add them
    // to the vectors that will be used to initialize the Target object.
    for ( int n = 0; n < num_nuclides; ++n ) {
      const auto& nuc = n_spec.at( n );
      bool ok = false;
      int nuc_pdg = nuc.to_long( ok );
      if ( ok ) atoms.emplace_back( nuc_pdg );
      // TODO: add support for string parsing
      //else if ( nuc.is_string() ) {
      //}
      else throw marley::Error( "Invalid target nuclide specifier "
        + nuc.dump_string() );

      // Parse and store the atom fraction. We already check for sane values
      // while initializing the Target object itself, so just make sure that the
      // conversion to a double worked out all right.
      const auto& frac_spec = af_spec.at( n );
      double frac = frac_spec.to_double( ok );
      if ( ok ) atom_fractions.push_back( frac );
      else throw marley::Error( "Invalid atom fraction "
        + frac_spec.dump_string() );
    }
  }

  // We're done. Create the new Target object and move it into the Generator.
  auto target = std::make_unique< marley::Target >( atoms, atom_fractions );
  if ( target->has_single_nuclide() ) {
    const marley::TargetAtom& ta = target->atom_fraction_map().cbegin()->first;
    MARLEY_LOG_INFO() << "Configured pure " << ta << " neutrino target";
  }
  else {
    MARLEY_LOG_INFO() << "Configured composite neutrino target with the"
      << " following nuclide fractions:\n" << *target;
  }
  gen.set_target( std::move(target) );
}

std::string marley::JSONConfig::source_get(const char* name,
  const marley::JSON& source_spec, const char* description,
  const char* default_str) const
{
  if ( !source_spec.has_key(name) ) {
    if ( default_str ) return default_str;
    else throw marley::Error( std::string("Missing source.") + name
      + " key for " + description + " source" );
  }
  bool ok;
  std::string result = source_spec.at( name ).to_string( ok );
  if ( !ok ) throw marley::Error( std::string("Invalid value given for source.")
    + name + " key for " + description + " source" );
  return result;
}

void marley::JSONConfig::update_logger_settings() const {

  using LogLevel = marley::Logger::LogLevel;
  marley::Logger& logger = marley::Logger::Instance();

  logger.clear_streams();

  if ( !json_.has_key("log") ) {
    // If the user hasn't specified a logger configuration, use the default,
    // which is logging at the INFO level to stdout.
    logger.add_stream( std::cout, LogLevel::INFO );
    logger.enable();
    return;
  }
  else {
    const marley::JSON& log_config = json_.at( "log" );

    if ( !log_config.is_array() ) {
      throw marley::Error( "The configuration given for the \"log\" key"
        " should be an array of JSON objects." );
    }

    auto elements = log_config.array_range();
    bool ok;
    // If the user has specified an empty list of logger files, then
    // disable the Logger and return immediately.
    if ( elements.begin() == elements.end() ) {
      logger.disable();
      return;
    }
    else logger.enable();

    // Loop over the list of log files and add them one-by-one to the Logger.
    for ( const auto& el : elements ) {
      // Get the file name for the new log file
      if ( !el.has_key("file") ) throw marley::Error( "Missing file name in"
        " a log file specification" );
      std::string file_name = el.at( "file" ).to_string( ok );
      if ( !ok ) throw marley::Error( "Invalid log file name \""
        + file_name + '\"' );

      // Use "info" as the default logging level
      LogLevel level = LogLevel::INFO;
      // Set the logging level for the current file to a non-default value
      // if the user has specified one. Complain if you get confused.
      if ( el.has_key("level") ) {
        std::string level_str = el.at( "level" ).to_string();
        if ( level_str == "error" ) level = LogLevel::ERROR;
        else if ( level_str == "warning" ) level = LogLevel::WARNING;
        else if ( level_str == "info" ) level = LogLevel::INFO;
        else if ( level_str == "debug" ) level = LogLevel::DEBUG;
        else if ( level_str == "disabled" ) level = LogLevel::DISABLED;
        else throw marley::Error( "Invalid logging level \""
          + el.dump_string() + '\"' );
      }

      // If the file name is "stdout", then add std::cout as a logging
      // stream. Do the same sort of thing for std::cerr. Otherwise, open the
      // requested file and add it to the logger streams.
      if ( file_name == "stdout" )
        logger.add_stream( std::cout, level );
      else if ( file_name == "stderr" )
        logger.add_stream( std::cerr, level );
      else {
        // If the user specified a value for the "overwrite" key, use it
        // to determine whether we should append to the file (false) or
        // overwrite it (true). Otherwise, assume we want to append to it.
        auto file_mode = std::ios::out;
        if ( el.has_key("overwrite") ) {

          marley::JSON ow = el.at( "overwrite" );

          bool overwrite = ow.to_bool( ok );
          if ( !ok ) throw marley::Error( "Invalid log file overwrite"
            " setting \"" + ow.dump_string() + '\"' );

          if ( overwrite ) file_mode |= std::ios::trunc;
          else file_mode |= std::ios::app;
        }
        else file_mode |= std::ios::app;

        auto outfile = std::make_shared< std::ofstream >( file_name,
          file_mode );
        if ( !outfile || (!outfile->good()) ) throw marley::Error( "Unable"
          " to open the log file \"" + file_name + "\"" );
        else logger.add_stream( outfile, level );
      }
    }
  }
}

void marley::JSONConfig::handle_json_error( const std::string& name,
  const marley::JSON& json )
{
  std::ostringstream message;
  message << "The JSON parameter \"" << name << "\" was set to the"
    << " invalid value " << json;
  throw marley::Error( message.str() );
}

bool marley::JSONConfig::process_extra_source_types(
  const std::string& type, const marley::JSON& source_spec, int pdg_code,
  std::unique_ptr<marley::NeutrinoSource>& source) const
{

#ifdef USE_ROOT
  if ( type != "th1" && type != "tgraph" ) return false;

  std::string tfile = source_get( "tfile", source_spec, type.c_str(), nullptr );
  std::string namecycle = source_get( "namecycle", source_spec, type.c_str(),
    nullptr );

  if ( type == "th1" ) {
    auto th1 = marley_root::get_root_object< TH1 >( tfile, namecycle );
    source = marley_root::make_root_neutrino_source( pdg_code, th1 );
    MARLEY_LOG_INFO() << "Created a TH1 "
      << marley_utils::neutrino_pdg_to_string( pdg_code )
      << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << source->get_Emin() << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << source->get_Emax() << " MeV";
    return true;
  }

  else if ( type == "tgraph" ) {
    auto tg = marley_root::get_root_object<TGraph>( tfile, namecycle );
    source = marley_root::make_root_neutrino_source( pdg_code, tg );
    MARLEY_LOG_INFO() << "Created a TGraph "
      << marley_utils::neutrino_pdg_to_string( pdg_code )
      << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << source->get_Emin() << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << source->get_Emax() << " MeV";
    return true;
  }
#endif

  return false;
}

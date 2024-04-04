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
#include <map>
#include <memory>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/Error.hh"
#include "marley/Reaction.hh"

namespace {

  constexpr int DUMMY_NUHEPMC_PROC_ID = 0;

  // G.R.4 and E.C.1
  struct NuHepMCProcess {
    NuHepMCProcess( int procID, std::string name, std::string description )
      : id_( procID ), name_( name ), desc_( description ) {}

    int id_;
    std::string name_;
    std::string desc_;
  };

  /// @todo Vary description for discrete/continuum processes?
  const std::map< marley::Reaction::ProcessType, NuHepMCProcess >
    ptype_to_nuhepmc_proc =
  {
    { marley::Reaction::ProcessType::NeutrinoCC_Discrete,
      { 100, "vCC", "inclusive charged-current neutrino-nucleus"
        " scattering" } },
    { marley::Reaction::ProcessType::AntiNeutrinoCC_Discrete,
      { 101, "anti-vCC", "inclusive charged-current antineutrino-nucleus"
        " scattering" } },
    { marley::Reaction::ProcessType::NC_Discrete,
      { 150, "NC", "inclusive neutral-current (anti)neutrino-nucleus"
        " scattering" } },
    { marley::Reaction::ProcessType::NuElectronElastic,
      { 700, "v-e", "(anti)neutrino-electron elastic scattering" } },
    { marley::Reaction::ProcessType::NeutrinoCC_Continuum,
      { 100, "vCC", "inclusive charged-current neutrino-nucleus"
        " scattering" } },
    { marley::Reaction::ProcessType::AntiNeutrinoCC_Continuum,
      { 101, "anti-vCC", "inclusive charged-current antineutrino-nucleus"
        " scattering" } },
    { marley::Reaction::ProcessType::NC_Continuum,
      { 150, "NC", "inclusive neutral-current (anti)neutrino-nucleus"
        " scattering" } },
  };

  // G.R.5
  std::map< int, std::pair< std::string, std::string > >
    vertex_status_map =
  {

    { marley_hepmc3::NUHEPMC_PRIMARY_VERTEX, { "Primary",
      "Represents the primary interaction" } },

    { marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX, { "HFDecay",
      "Represents a nuclear de-excitation step simulated using the"
      " Hauser-Feshbach treatment" } },

    { marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX, { "GammaDecay",
      "Represents a nuclear de-excitation step simulated using tabulated"
      " gamma-ray branching ratios" } },

  };

  // G.R.6 and P.R.1
  std::map< int, std::pair< std::string, std::string > >
    particle_status_map =
  {

    { marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, { "Final-state",
      "Undecayed physical particle" } },

    { marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, { "Projectile",
      "Incoming beam particle" } },

    { marley_hepmc3::NUHEPMC_TARGET_STATUS, { "Target",
      "Target particle struck by incoming beam particle in"
      " the primary interaction" } },

    { marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, { "UndecayedRemnant",
      "Nuclear remnant before simulation of nuclear de-excitations" } },

    { marley_hepmc3::NUHEPMC_INTERMEDIATE_RESIDUE_STATUS,
      { "IntermediateRemnant",
        "Nuclear remnant during simulation of nuclear de-excitations" } },

  };

  // G.C.1
  std::vector< std::string > nuhepmc_convention_vec = {
    "G.C.1", "G.C.4", "G.C.5", "G.C.6", "E.C.1", "E.C.2", "E.C.3"
  };

}

namespace marley_hepmc3 {

  int get_nuhepmc_proc_id( const marley::Reaction::ProcessType pt ) {
    auto itr = ptype_to_nuhepmc_proc.find( pt );
    if ( itr != ptype_to_nuhepmc_proc.end() ) {
      return itr->second.id_;
    }
    return DUMMY_NUHEPMC_PROC_ID;
  }

  void set_particle_charge( HepMC3::GenParticle& particle, int charge ) {
    bool added_ok = particle.add_attribute( "charge",
      std::make_shared< HepMC3::IntAttribute >(charge) );
    if ( !added_ok ) {
      throw marley::Error( "Failed to set particle charge in marley_hepmc3"
        "::set_particle_charge()" );
    }
  }

  int get_particle_charge( HepMC3::GenParticle& particle ) {
    // Return the charge stored in the particle attributes if it is set
    auto* q_ptr = particle.attribute< HepMC3::IntAttribute >( "charge" ).get();
    if ( q_ptr ) return q_ptr->value();
    // Otherwise, look up the charge based on the PDG code
    return marley_utils::get_particle_charge( particle.pid() );
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    const HepMC3::FourVector& mom4, int pdg, int status,
    double mass )
  {
    auto particle = std::make_shared< HepMC3::GenParticle >(
      mom4, pdg, status );
    if ( mass != DUMMY_PARTICLE_MASS ) {
      particle->set_generated_mass( mass );
    }
    return particle;
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, double px, double py, double pz, double E, int status,
    double mass )
  {
    HepMC3::FourVector mom4( px, py, pz, E );
    return make_particle( mom4, pdg, status, mass );
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, double px, double py, double pz, int status,
    double mass )
  {
    double E = marley_utils::real_sqrt( px*px + py*py + pz*pz + mass*mass );
    HepMC3::FourVector mom4( px, py, pz, E );
    return make_particle( mom4, pdg, status, mass );
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, int status, double mass )
  {
    HepMC3::FourVector mom4;
    if ( mass != DUMMY_PARTICLE_MASS ) {
      mom4.set_e( mass );
    }
    return make_particle( mom4, pdg, status, mass );
  }



  std::vector< std::shared_ptr< HepMC3::GenParticle > >
    get_particles_with_status( int status, HepMC3::GenEvent& ev )
  {
    const auto& particles = ev.particles();
    std::vector< std::shared_ptr< HepMC3::GenParticle > > found_particles;
    for ( auto& p : particles ) {
      if ( p->status() == status ) {
        found_particles.push_back( p );
      }
    }
    return found_particles;
  }

  std::vector< std::shared_ptr< HepMC3::GenVertex > >
    get_vertices_with_status( int status, HepMC3::GenEvent& ev )
  {
    const auto& vertices = ev.vertices();
    std::vector< std::shared_ptr< HepMC3::GenVertex > > found_vertices;
    for ( auto& v : vertices ) {
      if ( v->status() == status ) {
        found_vertices.push_back( v );
      }
    }
    return found_vertices;
  }

  std::shared_ptr< HepMC3::GenParticle >
    get_first_particle_with_status( int status, HepMC3::GenEvent& ev )
  {
    auto particle_ptrs = get_particles_with_status( status, ev );
    if ( !particle_ptrs.empty() ) return particle_ptrs.front();
    return nullptr;
  }

  std::shared_ptr< HepMC3::GenParticle > get_projectile(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, ev );
  }

  std::shared_ptr< HepMC3::GenParticle > get_target(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_TARGET_STATUS, ev );
  }

  std::shared_ptr< HepMC3::GenParticle > get_ejectile(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, ev );
  }

  std::shared_ptr< HepMC3::GenParticle > get_residue(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, ev );
  }

  // G.R.4
  void prepare_process_metadata( HepMC3::GenRunInfo& run_info ) {
    std::vector< int > proc_id_vec;
    for ( const auto& pair : ptype_to_nuhepmc_proc ) {
      int proc_id = pair.second.id_;
      proc_id_vec.push_back( proc_id );

      std::string attr_prefix = "NuHepMC.ProcessInfo["
        + std::to_string( proc_id ) + "].";

      run_info.add_attribute( attr_prefix + "Name",
        std::make_shared< HepMC3::StringAttribute >(pair.second.name_) );

      run_info.add_attribute( attr_prefix + "Description",
        std::make_shared< HepMC3::StringAttribute >(pair.second.desc_) );
    }

    run_info.add_attribute( "NuHepMC.ProcessIDs",
      std::make_shared< HepMC3::VectorIntAttribute >(proc_id_vec) );
  }

  // G.R.5
  void prepare_vertex_status_metadata( HepMC3::GenRunInfo& run_info ) {
    std::vector< int > status_vec;
    for ( const auto& pair : vertex_status_map ) {
      int status = pair.first;
      status_vec.push_back( status );

      std::string attr_prefix = "NuHepMC.VertexStatusInfo["
        + std::to_string( status ) + "].";

      run_info.add_attribute( attr_prefix + "Name",
        std::make_shared< HepMC3::StringAttribute >(pair.second.first) );

      run_info.add_attribute( attr_prefix + "Description",
        std::make_shared< HepMC3::StringAttribute >(pair.second.second) );
    }

    run_info.add_attribute( "NuHepMC.VertexStatusIDs",
      std::make_shared< HepMC3::VectorIntAttribute >(status_vec) );
  }

  // G.R.6
  void prepare_particle_status_metadata( HepMC3::GenRunInfo& run_info ) {
    std::vector< int > status_vec;
    for ( const auto& pair : particle_status_map ) {
      int status = pair.first;
      status_vec.push_back( status );

      std::string attr_prefix = "NuHepMC.ParticleStatusInfo["
        + std::to_string( status ) + "].";

      run_info.add_attribute( attr_prefix + "Name",
        std::make_shared< HepMC3::StringAttribute >(pair.second.first) );

      run_info.add_attribute( attr_prefix + "Description",
        std::make_shared< HepMC3::StringAttribute >(pair.second.second) );
    }

    run_info.add_attribute( "NuHepMC.ParticleStatusIDs",
      std::make_shared< HepMC3::VectorIntAttribute >(status_vec) );
  }

  // G.R.8
  // MARLEY currently doesn't use any non-standard PDG codes, so this step is
  // trivial
  void prepare_non_standard_pdg_code_metadata( HepMC3::GenRunInfo& run_info )
  {
    // Empty list
    std::vector< int > non_standard_PDGs;
    run_info.add_attribute( "NuHepMC.AdditionalParticleNumbers",
      std::make_shared< HepMC3::VectorIntAttribute >(non_standard_PDGs) );

    // No names or descriptions needed since the list is empty
  }

  // NOTE: the input flux-averaged total cross section is assumed to be in
  // natural units (MeV^{-2})
  void apply_nuhepmc_runinfo_conventions( HepMC3::GenRunInfo& run_info,
    const double flux_avg_xsec )
  {

    // G.C.1
    run_info.add_attribute( "NuHepMC.Conventions",
      std::make_shared< HepMC3::VectorStringAttribute >(
        nuhepmc_convention_vec )
    );

    // G.C.4
    run_info.add_attribute( "NuHepMC.Units.CrossSection.Unit",
      std::make_shared< HepMC3::StringAttribute >( "pb" )
    );

    run_info.add_attribute( "NuHepMC.Units.CrossSection.TargetScale",
      std::make_shared< HepMC3::StringAttribute >( "PerTargetAtom" )
    );

    // G.C.5
    double xsec_picobarn = flux_avg_xsec * marley_utils::hbar_c2
      * marley_utils::fm2_to_picobarn;
    run_info.add_attribute( "NuHepMC.FluxAveragedTotalCrossSection",
      std::make_shared< HepMC3::DoubleAttribute >( xsec_picobarn )
    );

    // G.C.6
    std::vector< std::string > marley_DOIs = {
      "10.1103/PhysRevC.103.044604",
      "10.1016/j.cpc.2021.108123"
    };

    run_info.add_attribute( "NuHepMC.Citations.Generator.DOI",
      std::make_shared< HepMC3::VectorStringAttribute >( marley_DOIs )
    );

    std::vector< std::string > marley_arXivs = {
      "2010.02393",
      "2101.11867"
    };

    run_info.add_attribute( "NuHepMC.Citations.Generator.arXiv",
      std::make_shared< HepMC3::VectorStringAttribute >( marley_arXivs )
    );

    std::vector< std::string > marley_INSPIREs = {
      "Gardiner:2020ulp",
      "Gardiner:2021qfr"
    };

    run_info.add_attribute( "NuHepMC.Citations.Generator.InspireHEP",
      std::make_shared< HepMC3::VectorStringAttribute >( marley_INSPIREs )
    );

  }

}

// Standard library includes
#include <iostream>
#include <string>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/NucleusDecayer.hh"

// De-excites the nuclear remnant(s) in every event of an existing HepMC3
// file and writes a new HepMC3 file with the decay products appended. This
// is meant for files produced by generators other than MARLEY (e.g. GENIE),
// so unlike the "marley" CLI's marley::EventFileReader/OutputFile classes
// (which require MARLEY-specific run-info attributes), this program reads
// and writes plain HepMC3 ASCII directly.
//
// A residue is found by marley::NucleusDecayer itself: it looks for HepMC3
// particles with status 27 (NuHepMC's "undecayed nuclear residue" code) in
// each event and de-excites every one it finds. Since the input event's
// GenRunInfo won't list a tool named "MARLEY" (it'll say "GENIE", or
// whatever the source generator was), NucleusDecayer automatically takes
// its "externally produced event" code path: it infers the excitation
// energy from how far the residue's mass sits above its isotope's
// ground-state mass, rather than reading it from MARLEY-native attributes.

int main( int argc, char* argv[] ) {

  if ( argc != 3 && argc != 4 ) {
    std::cout << "Usage: " << argv[0]
      << " INPUT_HEPMC3_FILE OUTPUT_HEPMC3_FILE [SEED]\n";
    return 1;
  }

  std::string input_file( argv[1] );
  std::string output_file( argv[2] );

  // Build a MARLEY Generator configured for de-excitation only, following
  // the same pattern MARLEY's own "marley decay" command uses (see
  // src/app/cmd_decay.cc): a JSON config with "reactions" and "source" both
  // explicitly null makes JSONConfig::create_generator() skip everything
  // related to simulating a primary neutrino interaction (cross sections,
  // target/projectile setup, etc.) and just build the mass table, nuclear
  // structure database, and random number generator that NucleusDecayer
  // needs.
  //
  // This is deliberately not marley::Generator's bare default constructor:
  // that constructor leaves an internal "weighter" object unset, which
  // crashes inside set_up_run_info() below. Routing through JSONConfig
  // (even with an empty config) sets it up correctly.
  marley::JSON json = marley::JSON::object();
  json["reactions"] = marley::JSON( nullptr );
  json["source"] = marley::JSON( nullptr );

  if ( argc == 4 ) {
    // An explicit seed makes the de-excitation reproducible across runs.
    // Without one, MARLEY seeds itself from the system clock.
    json["seed"] = std::stoll( argv[3] );
  }

  marley::JSONConfig config( json );
  marley::Generator gen = config.create_generator();
  gen.set_up_run_info();

  HepMC3::ReaderAscii reader( input_file );
  HepMC3::WriterAscii writer( output_file );

  HepMC3::GenEvent event;
  long n_events = 0;
  long n_particles_added = 0;

  // Note: read_event()'s own return value is NOT a reliable end-of-file
  // signal for HepMC3::ReaderAscii -- MARLEY's own EventFileReader class
  // (src/EventFileReader.cc) never trusts it alone either; it separately
  // checks the reader's failed() state. Without this, the loop below would
  // keep "succeeding" past the last real event, reprocessing a stale event
  // object forever.
  while ( !reader.failed() ) {

    reader.read_event( event );
    if ( reader.failed() ) break;

    size_t n_particles_before = event.particles().size();

    // NucleusDecayer finds every status-27 (undecayed residue) particle in
    // the event on its own and de-excites each one, appending the decay
    // products (gammas, nucleons, light fragments) as new particles and
    // vertices directly into "event".
    marley::NucleusDecayer decayer;
    try {
      decayer.process_event( event, gen );
    }
    catch ( const marley::Error& error ) {
      std::cerr << "MARLEY de-excitation failed on event " << n_events
        << ": " << error.what() << '\n';
      return 1;
    }

    n_particles_added += event.particles().size() - n_particles_before;

    writer.write_event( event );
    ++n_events;
  }

  reader.close();
  writer.close();

  std::cout << "Processed " << n_events << " event(s), added "
    << n_particles_added << " de-excitation product particle(s).\n";

  return 0;
}

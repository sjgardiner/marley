# Executable examples

For a variety of applications, it may be helpful to create external C++
programs that make use of one or more MARLEY classes. To aid users in building
such programs, a Bash script called *marley-config* is generated and placed in
the build/bin/ folder whenever MARLEY is successfully built. The three programs
in this folder may each be compiled using the marley-config script by executing
a command of the form

```
g++ -o EXECUTABLE_NAME $(marley-config --cflags --libs) SOURCE_NAME
```

where EXECUTABLE_NAME is the file name for the compiled executable and
SOURCE_NAME is the name of a C++ source file in this folder. The syntax shown
above assumes that the build/bin/ folder has previously been added to the system
PATH (e.g., by sourcing the setup_marley.sh script in the root source code
folder).

The four example programs provided in this folder are

* efr.cc: Reads a single MARLEY output file specified as the first command-line
  argument. The flux-averaged total cross section is printed to stdout,
  followed by each event in the HepMC3 text format.

* evgen.cc: Generates 10 MARLEY events using the job configuration file
  /home/config.js. This file name is currently hard-coded. Each event is
  printed to stdout in the HepMC3 text format after being generated.

* mass_40Ar.cc: Prints a brief message to stdout with two pieces of
  information: (1) whether or not MARLEY was built with ROOT support
  enabled, and (2) the atomic mass of 40Ar in MeV/c^2

* marley_deexcite_hepmc.cc: Reads a plain HepMC3 file produced by another
  generator (e.g. GENIE) whose events contain nuclear remnants tagged with
  NuHepMC status 27 ("undecayed residue"), runs MARLEY's de-excitation on
  each one, and writes a new HepMC3 file with the decay products appended.
  Unlike efr.cc, it does not go through marley::EventFileReader, since that
  class expects MARLEY-native run-info attributes that a file from another
  generator won't have. Usage:
  `marley_deexcite_hepmc INPUT_HEPMC3_FILE OUTPUT_HEPMC3_FILE [SEED]`

An example of using the marley-config script in a Makefile may be seen in
examples/marg4/Makefile.

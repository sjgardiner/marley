CXX=g++
CXXFLAGS=-g -O3 -std=c++14 -I. -Wall -Wextra -Wpedantic -Werror -Wno-error=unused-parameter
USE_ROOT=yes

OBJ = marley_utils.o meta_numerics.o TMarleyParticle.o TMarleyEvent.o
OBJ += TMarleyGenerator.o TMarleyReaction.o TMarleyNuclearReaction.o
OBJ += TMarleyElectronReaction.o TMarleyGamma.o TMarleyLevel.o
OBJ += TMarleyDecayScheme.o TMarleyMassTable.o TMarleyStructureDatabase.o
OBJ += TMarleyConfigFile.o TMarleyNuclearPhysics.o
OBJ += TMarleyBackshiftedFermiGasModel.o TMarleySphericalOpticalModel.o
OBJ += TMarleyNeutrinoSource.o TMarleyKinematics.o TMarleyDecayChannel.o

ifdef USE_ROOT
# Adding the g++ compiler option -DUSE_ROOT to the CXXFLAGS
# variable allows you to use conditional compilation
# via the preprocessor directive #ifdef USE_ROOT.
# Currently none of the core MARLEY classes use such
# preprocessor directives, but the example executable
# react does.
#
# The root_dict.o object file should be added
# to the list of prerequisites for any executable
# that uses TTrees containing TMarleyEvents
# or TMarleyParticles
CXXFLAGS += `root-config --cflags` -DUSE_ROOT
OBJ_DICT = root_dict.o
LDFLAGS=`root-config --libs`
endif

all: parse react validate check

%.o: %.c
	$(CXX) -c -o $@

parse: $(OBJ) parse.o
	$(CXX) -o $@ $^

react: $(OBJ) $(OBJ_DICT) react.o
	$(CXX) -o $@ $^ $(LDFLAGS)

validate: $(OBJ) $(OBJ_DICT) validate.o
	$(CXX) -o $@ $^ $(LDFLAGS)

check: $(OBJ) $(OBJ_DICT) check.o
	$(CXX) -o $@ $^ $(LDFLAGS)

brs: $(OBJ) $(OBJ_DICT) brs.o
	$(CXX) -o $@ $^ $(LDFLAGS)

brs2: $(OBJ) $(OBJ_DICT) brs2.o
	$(CXX) -o $@ $^ $(LDFLAGS)

brs3: $(OBJ) $(OBJ_DICT) brs3.o
	$(CXX) -o $@ $^ $(LDFLAGS)

brs4: $(OBJ) $(OBJ_DICT) brs4.o
	$(CXX) -o $@ $^ $(LDFLAGS)

hf: $(OBJ) $(OBJ_DICT) hf.o
	$(CXX) -o $@ $^ $(LDFLAGS)

plots: $(OBJ) $(OBJ_DICT) plots.o
	$(CXX) -o $@ $^ $(LDFLAGS)

dump: $(OBJ) $(OBJ_DICT) dump.o
	$(CXX) -o $@ $^ $(LDFLAGS)

nu_source_plot: $(OBJ) nu_source_plot.o
	$(CXX) -o $@ $^ $(LDFLAGS)

exs: $(OBJ) $(OBJ_DICT) exs.o
	$(CXX) -o $@ $^ $(LDFLAGS)

# Add more header files to the prerequisites for
# root_dict.o if you would like to store other
# MARLEY classes in ROOT trees. All such classes
# currently use a single automatically-generated
# dictionary source file root_dict.cc
# 
# The commands invoked to create root_dict.o
# do the following things:
# 1. Remove old dictionary files
# 2. Create new dictionary files, enabling ROOT
#    i/o by adding the '+' suffix to each prerequisite
#    header file (.hh extension)
# 3. Compile the dictionary source file
root_dict.o: TMarleyParticle.hh TMarleyEvent.hh
	rm -f root_dict.cc root_dict.h
	rootcint root_dict.cc -c $(subst .hh,.hh+,$^)
	$(CXX) $(CXXFLAGS) -c -o root_dict.o root_dict.cc

.PHONY: clean

clean:
	rm -f *.o parse react validate brs brs2 brs3 brs4 root_dict.cc
	rm -f *.o root_dict.h nu_source_plot dump

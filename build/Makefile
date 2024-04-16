# Names of various MARLEY directories
TOP_DIR = $(shell cd .. && pwd)
BUILD_DIR = $(TOP_DIR)/build
DATA_DIR = $(TOP_DIR)/data
INCLUDE_DIR = $(TOP_DIR)/include
SRC_DIR = $(TOP_DIR)/src
HEPMC3_DIR = $(TOP_DIR)/hepmc3

# Default version of the C++ standard to use for compilation
CXX_STD=c++17

UNKNOWN_REV=unknown version

# If the target is "debug", then include debugging information and build the
# marley executable with optimization turned off. Prepend the default -std
# option so that it will be overridden if the user has manually specified one
# in CXXFLAGS.
ifeq ($(MAKECMDGOALS),debug)
  override CXXFLAGS := -std=$(CXX_STD) $(CXXFLAGS) -O0 -g
else
  # Otherwise, use full optimization and do not include debugging info. Prepend
  # the default -O3 and -std options here in case the user wants to manually
  # override them in CXXFLAGS.
  override CXXFLAGS := -O3 -std=$(CXX_STD) $(CXXFLAGS)
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHARED_LIB_SUFFIX=dylib
else ifeq ($(UNAME_S),Linux)
  SHARED_LIB_SUFFIX=so
else
  $(warning Unrecognized operating system encountered.)
  SHARED_LIB_SUFFIX=so
endif

# If available, use git to determine the hash of the current version
# of MARLEY. If not, set it to the string "unknown version"

# Check if git is available on the system path.
GIT := $(shell command -v git 2> /dev/null)

ifndef GIT
  GIT_REVISION := $(UNKNOWN_REV)
else
  # Also check that we're inside of a folder managed by a git repository. If
  # we're not, then git rev-parse won't work
  GIT_REVPARSE_CODE := $(shell git rev-parse 2> /dev/null && echo "$$?")
  ifeq ($(GIT_REVPARSE_CODE),0)

    GIT_REVISION := $(shell git rev-parse --short HEAD)

    # Verify that the working tree and index are identical to HEAD. If they
    # are not (because there are uncommitted changes to at least one tracked
    # file), then append "-dirty" to the git revision variable.
    GIT_DIFF_INDEX_CODE := $(shell git diff-index --quiet HEAD \
      2> /dev/null && echo "$$?")

    ifneq ($(GIT_DIFF_INDEX_CODE),0)
      GIT_REVISION := $(GIT_REVISION)-dirty
    endif

  else
    GIT_REVISION := $(UNKNOWN_REV)
  endif

endif

# If the .VERSION file exists, use its contents as the
# MARLEY version number. The .VERSION file should only
# be added to the source tree in tagged releases.
ifneq (,$(wildcard ../.VERSION))
  MARLEY_VERSION := $(shell cat ../.VERSION)

  # Define a link to a tarball on GitHub for the current tagged release
  TARBALL_LINK = "<a href=\"https://github.com/MARLEY-MC/marley/$\
    archive/v$(MARLEY_VERSION).tar.gz\">here</a> or"

  VERSION_PREFIX="v"

else

  # Just link to the repository in the doxygen documentation if we're not
  # working with a tagged release
  TARBALL_LINK = ""

  # Also use the git revision as the version number in this case
  MARLEY_VERSION := $(GIT_REVISION)

endif

# Define the MARLEY_VERSION preprocessor macro
override CXXFLAGS += -DMARLEY_VERSION="\"$(MARLEY_VERSION)\""

# Define the MARLEY_GIT_REVISION preprocessor macro
override CXXFLAGS += -DMARLEY_GIT_REVISION="\"$(MARLEY_GIT_REVISION)\""

SHARED_LIB_NAME := MARLEY
SHARED_LIB := lib$(SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)
ROOT_SHARED_LIB_NAME := MARLEY_ROOT
ROOT_SHARED_LIB := lib$(ROOT_SHARED_LIB_NAME).$(SHARED_LIB_SUFFIX)

TEST_EXECUTABLE = martest
TEST_OBJECTS = $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/tests/*.cc)))

all: marley
debug: marley
test: $(TEST_EXECUTABLE)

# Skip lots of initialization if all we want is "make clean/uninstall"
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),uninstall)

  # Use g++ as the default compiler
  CXX ?= g++
  CXXFLAGS += -I$(INCLUDE_DIR) -Wall -Wextra -Wpedantic
  CXXFLAGS += -Wno-error=unused-parameter -Wcast-align

  # Add extra compiler flags for recognized compilers (currently just gcc
  # and clang)
  CXXVERSION = $(shell $(CXX) --version)
  COMPILER_VERSION := $(word 3, $(CXXVERSION))
  ifneq (,$(findstring clang,$(CXXVERSION)))
    # clang
    $(info Compiling using version $(COMPILER_VERSION) of clang)

    # The ROOT headers trigger clang's no-keyword-macro warning, so
    # disable it. Also disable (for now) warnings for braces around
    # initialization of subobjects (overkill in the meta_numerics header)
    CXXFLAGS += -Wno-keyword-macro -Wno-missing-braces
  else
    ifneq (,$(or $(findstring GCC,$(CXXVERSION)), $(findstring g++,$(CXXVERSION))))
      # gcc
      $(info Compiling using version $(COMPILER_VERSION) of GCC)
      ifneq (,$(findstring $(COMPILER_VERSION), 4.9.))
        # g++ 4.9 gives many false positives for -Wshadow, so disable it
        # for now.
        override CXXFLAGS += -Wno-shadow
      endif
      # Linking to ROOT libraries can be problematic on distributions (e.g.,
      # Ubuntu) that set the g++ flag -Wl,--as-needed by default (see
      # http://www.bnikolic.co.uk/blog/gnu-ld-as-needed.html for details), so
      # disable this behavior on Linux.
      ifneq ($(UNAME_S),Darwin)
        override CXXFLAGS += -Wl,--no-as-needed
      endif
    endif
  endif
  
  OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(SRC_DIR)/*.cc)))
  OBJECTS := $(filter-out marley.o marley_root.o, $(OBJECTS))
  OBJECTS := $(filter-out marsum.o OutputFileRoot.o, $(OBJECTS))
  OBJECTS := $(filter-out OutputFilePlainRoot.o, $(OBJECTS))
  OBJECTS := $(filter-out MacroEventFileReader.o, $(OBJECTS))

  # Get information about the GNU Scientific Library installation
  # (required as of MARLEY v1.1.0)
  GSLCONFIG := $(shell which gsl-config)
  ifeq (, $(shell which gsl-config))
    define err_message
Could not find a working GNU Scientific Library installation.
Please ensure that the gsl-config executable is on your PATH and try again.
    endef
    $(error "$(err_message)")
  endif

  # Add the appropriate compilation flags to use the GNU Scientific Library
  GSL_CXXFLAGS := $(shell $(GSLCONFIG) --cflags)
  GSL_LDFLAGS := $(shell $(GSLCONFIG) --libs)

  # Get information about the HepMC3 installation
  ifndef IGNORE_HEPMC3
    HEPMC3CONFIG := $(shell which HepMC3-config)
    ifeq (, $(shell which HepMC3-config))
      define hepmc3_builtin_message
    Could not find a working HepMC3 installation.
    A built-in HepMC3 library will be used instead.
      endef
      $(info "$(hepmc3_builtin_message)")
    else
      # Add the appropriate compilation flags to use the HepMC3 reference
      # implementation
      HEPMC3_INCDIR := $(shell $(HEPMC3CONFIG) --includedir)
      HEPMC3_CXXFLAGS := $(shell $(HEPMC3CONFIG) --cflags)
      HEPMC3_LDFLAGS := $(shell $(HEPMC3CONFIG) --libs)
      FOUND_HEPMC3 := TRUE
    endif
  else
    $(info Ignoring any HepMC3 installations that may be present.)
  endif

  ifndef FOUND_HEPMC3
    # Use the built-in HepMC3 library
    HEPMC3_INCDIR := $(HEPMC3_DIR)/include
    HEPMC3_CXXFLAGS := -I$(HEPMC3_INCDIR)
    HEPMC3_LDFLAGS := -L$(BUILD_DIR) -lHepMC3

    HEPMC3_OBJECTS := $(notdir $(patsubst %.cc,%.o,$(wildcard $(HEPMC3_DIR)/*.cc)))
    HEPMC3_SHARED_LIB := libHepMC3.$(SHARED_LIB_SUFFIX)

$(HEPMC3_SHARED_LIB): $(HEPMC3_OBJECTS)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) \
        $(HEPMC3_CXXFLAGS) -fPIC -shared $(ROOT_LDFLAGS) $(GSL_LDFLAGS) \
        -o $@ $^

  endif

  # The user may force the Makefile to ignore ROOT entirely by defining
  # IGNORE_ROOT="yes" (or any non-empty string) on the command line
  # invocation of make.
  ifndef IGNORE_ROOT
    ROOTCONFIG := $(shell command -v root-config 2> /dev/null)
    # prefer rootcling as the dictionary generator executable name, but use
    # rootcint if you can't find it
    ROOTCLING := $(shell command -v rootcling 2> /dev/null)
    ifndef ROOTCLING
      ROOTCLING := $(shell command -v rootcint 2> /dev/null)
    endif
    ROOT := $(shell command -v root 2> /dev/null)

    ifndef ROOTCONFIG
      $(info WARNING: Could not find a valid ROOT installation.)
      $(info MARLEY will be built without ROOT support.)
      USE_ROOT = no
      CXXFLAGS += -std=$(CXX_STD)
    else
      ROOT_VERSION := $(shell $(ROOTCONFIG) --version)
      $(info Found ROOT version $(ROOT_VERSION) in $(ROOT))
      $(info MARLEY will be built with ROOT support.)
      override CXXFLAGS += -DUSE_ROOT
      MAYBE_MARSUM = marsum mroot
      ROOT_CXXFLAGS := $(shell $(ROOTCONFIG) --cflags)

      # If ROOT was built with a later C++ standard, switch to building MARLEY
      # with it as well, just in case. Later -std options take precedence, so
      # just tack on a new one to keep things simple.
      ifneq (, $(findstring c++2a, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++2a
        override CXXFLAGS += -std=$(CXX_STD)
      endif
      ifneq (, $(findstring c++20, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++20
        override CXXFLAGS += -std=$(CXX_STD)
      endif
      ifneq (, $(findstring c++2b, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++2b
        override CXXFLAGS += -std=$(CXX_STD)
      endif
      ifneq (, $(findstring c++23, $(ROOT_CXXFLAGS)))
        override CXX_STD = c++23
        override CXXFLAGS += -std=$(CXX_STD)
      endif

      ROOT_LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
      ROOT_LIBDIR := $(shell $(ROOTCONFIG) --libdir)
      ROOT_LDFLAGS += -L$(ROOT_LIBDIR) -lCore -lRIO -lHist -lTree
      ifeq ($(UNAME_S),Linux)
        ROOT_LDFLAGS += -rdynamic
      endif
      ROOT_DICT_INCLUDES := -I$(INCLUDE_DIR) -I$(HEPMC3_INCDIR) \
        -I$(HEPMC3_INCDIR)/HepMC3/Data \
        $(notdir $(wildcard $(HEPMC3_INCDIR)/HepMC3/Data/*.h)) \
        marley/marley_linkdef.hh
      USE_ROOT = yes
      ROOT_OBJ_DICT = marley_root_dict.o

      # If we're linking against ROOT, then we need to move the OutputFile.o,
      # OutputFileAscii.o, EventFileReader.o objects from the MARLEY library to the
      # MARLEY_ROOT library
      OBJECTS := $(filter-out OutputFile.o OutputFileAscii.o, $(OBJECTS))
      OBJECTS := $(filter-out EventFileReader.o JSONConfig.o, $(OBJECTS))

      ROOT_SHARED_LIB_LDFLAGS := -l$(ROOT_SHARED_LIB_NAME)
      ROOT_SHARED_LIB_OBJECTS = marley_root.o JSONConfig.o OutputFile.o
      ROOT_SHARED_LIB_OBJECTS += OutputFileAscii.o OutputFileRoot.o
      ROOT_SHARED_LIB_OBJECTS += OutputFilePlainRoot.o
      ROOT_SHARED_LIB_OBJECTS += EventFileReader.o $(ROOT_OBJ_DICT)
      #ROOT_SHARED_LIB_OBJECTS += MacroEventFileReader.o
  
$(ROOT_OBJ_DICT):
	$(RM) marley_root_dict*.*
	$(ROOTCLING) -f marley_root_dict.cc -c $(ROOT_DICT_INCLUDES)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	  -I$(INCLUDE_DIR) -I$(HEPMC3_INCDIR)/HepMC3/Data -fPIC \
          -o $(ROOT_OBJ_DICT) -c marley_root_dict.cc
	$(RM) marley_root_dict.cc
  
$(ROOT_SHARED_LIB): $(SHARED_LIB) $(ROOT_SHARED_LIB_OBJECTS)
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -fPIC -shared \
	  -o $(ROOT_SHARED_LIB) -L. -l$(SHARED_LIB_NAME) \
	  $(ROOT_SHARED_LIB_OBJECTS) $(ROOT_LDFLAGS) $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS)

    endif
  else
    $(info Ignoring any ROOT installations that may be present.)
    $(info MARLEY will be built without ROOT support.)
    USE_ROOT = no
  endif
endif
endif

# Set up default variables for the install/uninstall targets and for setting
# the executable rpath
prefix = /usr
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
datadir = $(prefix)/share
libdir = $(exec_prefix)/lib
incdir = $(prefix)/include

MARLEY_LIBS := $(SHARED_LIB)
ifeq ($(USE_ROOT),yes)
  MARLEY_LIBS += $(ROOT_SHARED_LIB)

  # If we're building the tests, then link to some extra ROOT libraries
  ifneq (,$(findstring test,$(MAKECMDGOALS)))
    ROOT_LDFLAGS += -lGraf -lGpad
  endif
endif

# Causes GNU make to auto-delete the object files when the build is complete
.INTERMEDIATE: $(OBJECTS) $(TEST_OBJECTS) $(ROOT_SHARED_LIB_OBJECTS) $(HEPMC3_OBJECTS) marley.o marsum.o

%.o: $(SRC_DIR)/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $@ -c $^

%.o: $(SRC_DIR)/tests/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $@ -c $^

%.o: $(HEPMC3_DIR)/%.cc
	$(CXX) $(ROOT_CXXFLAGS) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) \
	-I$(INCLUDE_DIR) -fPIC -o $@ -c $^

$(SHARED_LIB): $(HEPMC3_SHARED_LIB) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) \
	-fPIC -shared -o $@ $(OBJECTS)

marsum: $(MARLEY_LIBS) marsum.o
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) -Wl,-rpath -Wl,$(libdir):$(shell pwd) marsum.o

mroot: $(MARLEY_LIBS)
	cp $(SRC_DIR)/scripts/mroot .

# We use a temporary backup file here so that the invocation of sed is
# compatible with both the GNU/Linux and BSD/macOS versions.
# See https://stackoverflow.com/a/22084103/4081973 for details.
marley-config: $(MARLEY_LIBS)
	$(RM) marley-config
	cp $(SRC_DIR)/scripts/marley-config.in marley-config
	sed -i.bak -e '/^##/d' -e "s|@@VERSION@@|\"$(MARLEY_VERSION)\"|g" \
	  -e "s|@@GIT_REVISION@@|\"$(GIT_REVISION)\"|g" \
	  -e "s|@@CXX_STD@@|\"$(CXX_STD)\"|g" \
	  -e "s|@@HEPMC3_CFLAGS@@|\"$(HEPMC3_CXXFLAGS)\"|g" \
	  -e "s|@@HEPMC3_LIBS@@|\"$(HEPMC3_LDFLAGS)\"|g" \
	  -e "s|@@USE_ROOT@@|\"$(USE_ROOT)\"|g" marley-config
	$(RM) marley-config.bak

marley: $(MARLEY_LIBS) marley.o marley-config $(MAYBE_MARSUM)
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) -Wl,-rpath -Wl,$(libdir):$(shell pwd) marley.o

$(TEST_EXECUTABLE): $(TEST_OBJECTS) $(MARLEY_LIBS)
	$(CXX) $(CXXFLAGS) $(GSL_CXXFLAGS) $(HEPMC3_CXXFLAGS) -o $@ -L. \
	  -l$(SHARED_LIB_NAME) $(ROOT_SHARED_LIB_LDFLAGS) $(ROOT_LDFLAGS) \
	  $(GSL_LDFLAGS) $(HEPMC3_LDFLAGS) -o $@ $(TEST_OBJECTS)

marg4: $(MARLEY_LIBS)
	$(RM) ../examples/marg4/build/marg4
	cd ../examples/marg4/build && $(MAKE)
	cp ../examples/marg4/build/marg4 .
	$(RM) ../examples/marg4/build/marg4

marprint: $(MARLEY_LIBS)
	$(RM) ../examples/executables/build/marprint
	cd ../examples/executables/build && $(MAKE) marprint
	cp ../examples/executables/build/marprint .
	$(RM) ../examples/executables/build/marprint

mardumpxs: $(MARLEY_LIBS)
	$(RM) ../examples/executables/build/mardumpxs
	cd ../examples/executables/build && $(MAKE) mardumpxs
	cp ../examples/executables/build/mardumpxs .
	$(RM) ../examples/executables/build/mardumpxs

.PHONY: docs clean install uninstall

doxygen:
	export MARLEY_VERSION=$(VERSION_PREFIX)$(MARLEY_VERSION) \
	export TARBALL_LINK=$(TARBALL_LINK) \
	&& cd ../docs && $(MAKE) doxygen

docs:
	export MARLEY_VERSION=$(VERSION_PREFIX)$(MARLEY_VERSION) \
	export TARBALL_LINK=$(TARBALL_LINK) \
	&& cd ../docs && $(MAKE) html

clean:
	$(RM) *.$(SHARED_LIB_SUFFIX) *.o marley_root_dict*.* marley
	$(RM) -rf *.dSYM marsum mroot $(TEST_EXECUTABLE) marg4
	$(RM) -rf marprint mardumpxs marley-config ../doxygen/html/*
	$(RM) -rf ../docs/_build/*

install: marley
	mkdir -p $(DESTDIR)$(bindir)
	mkdir -p $(DESTDIR)$(libdir)
	mkdir -p $(DESTDIR)$(incdir)/marley
	mkdir -p $(DESTDIR)$(datadir)/marley
	cp marley $(MAYBE_MARSUM) $(DESTDIR)$(bindir)
	cp $(SHARED_LIB) $(DESTDIR)$(libdir)
	cp $(ROOT_SHARED_LIB) $(DESTDIR)$(libdir) 2> /dev/null || true
	cp marley_root_dict_rdict.pcm $(DESTDIR)$(libdir) 2> /dev/null || true
	cp -r ../react $(DESTDIR)$(datadir)/marley
	cp -r ../structure $(DESTDIR)$(datadir)/marley
	cp -r ../examples $(DESTDIR)$(datadir)/marley
	cp -r ../include/marley $(DESTDIR)$(incdir)
	ldconfig

uninstall:
	$(RM) $(DESTDIR)$(bindir)/marley
	$(RM) $(DESTDIR)$(bindir)/marsum
	$(RM) $(DESTDIR)$(bindir)/mroot
	$(RM) $(DESTDIR)$(libdir)/$(SHARED_LIB)
	$(RM) $(DESTDIR)$(libdir)/$(ROOT_SHARED_LIB)
	$(RM) $(DESTDIR)$(libdir)/marley_root_dict_rdict.pcm
	$(RM) -r $(DESTDIR)$(datadir)/marley
	$(RM) -r $(DESTDIR)$(incdir)/marley
	ldconfig

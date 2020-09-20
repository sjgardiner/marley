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
#include <string>
#include <unordered_map>

namespace marley {

  // Forward-declare some classes
  class Fragment;
  class JSON;

  /// @brief Singleton lookup table for particle and atomic masses
  class MassTable {

    public:

      /// @brief Deleted copy constructor
      MassTable(const MassTable&) = delete;

      /// @brief Deleted move constructor
      MassTable(MassTable&&) = delete;

      /// @brief Deleted copy assignment operator
      MassTable& operator=(const MassTable&) = delete;

      /// @brief Deleted move assignment operator
      MassTable& operator=(MassTable&&) = delete;

      /// @brief Get a const reference to the singleton instance of the
      /// MassTable
      static const MassTable& Instance();

      /// @brief Get the mass of a particle
      /// @param pdg_code PDG code identifying the type of particle
      /// @return %Particle mass (MeV)
      double get_particle_mass(int pdg_code) const;

      /// @brief Get the mass of an atom
      /// @param pdg_code PDG code identifying the nucleus of the atom
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Atomic mass (MeV)
      double get_atomic_mass(int pdg_code, bool theory_ok = true) const;

      /// @brief Get the mass of an atom
      /// @param Z Atomic number of the atom
      /// @param A Mass number of the atom
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Atomic mass (MeV)
      double get_atomic_mass(int Z, int A, bool theory_ok = true) const;

      /// @brief Get the separation energy for emission of a nuclear fragment
      /// from a nucleus
      /// @param Z Atomic number for the mother nucleus
      /// @param A Mass number for the mother nucleus
      /// @param pdg PDG code for the emitted nuclear fragment
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Separation energy (MeV)
      double get_fragment_separation_energy(int Z, int A, int pdg,
        bool theory_ok = true) const;

      /// @brief Get the separation energy for emission of a nuclear fragment
      /// from a nucleus
      /// @param nuc_pdg PDG code for the mother nucleus
      /// @param frag_pdg PDG code for the emitted nuclear fragment
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Separation energy (MeV)
      double get_fragment_separation_energy(int nuc_pdg, int frag_pdg,
        bool theory_ok = true) const;

      /// @brief Get the binding energy of a nucleus
      /// @param Z Atomic number
      /// @param A Mass number
      /// @param theory_ok Whether to use a theoretical (liquid drop model)
      /// atomic mass during the calculation if an experimental atomic mass
      /// cannot be found. If theory_ok is false, then a marley::Error will be
      /// thrown if an experimental atomic mass cannot be found.
      /// @return Binding energy (MeV)
      double get_binding_energy(int Z, int A, bool theory_ok = true) const;

      /// @brief Get the <a href="https://en.wikipedia.org/wiki/Mass_excess">
      /// mass excess</a> of a nucleus
      /// @param Z Atomic number
      /// @param A Mass number
      /// @param theory_ok Whether to use a theoretical (liquid drop model)
      /// atomic mass during the calculation if an experimental atomic mass
      /// cannot be found. If theory_ok is false, then a marley::Error will be
      /// thrown if an experimental atomic mass cannot be found.
      /// @return Mass excess (MeV)
      double get_mass_excess(int Z, int A, bool theory_ok = true) const;

      /// @brief Calculate a theoretical
      /// <a href="https://en.wikipedia.org/wiki/Mass_excess">mass excess</a>
      /// for a nucleus using the liquid drop model
      /// @param Z Atomic number
      /// @param A Mass number
      /// @return Mass excess (MeV)
      double liquid_drop_model_mass_excess(int Z, int A) const;

      /// @brief Calculate a theoretical atomic mass using the liquid drop model
      /// @param Z Atomic number
      /// @param A Mass number
      /// @return Atomic mass (MeV)
      double liquid_drop_model_atomic_mass(int Z, int A) const;

      /// @brief Get the approximate excitation energy threshold for emission
      /// of a particular nuclear fragment
      /// @param Zi Proton number of the initial nucleus
      /// @param Ai Mass number of the initial nucleus
      /// @param f Fragment to be emitted
      double fragment_emission_threshold(const int Zi, const int Ai,
        const marley::Fragment& f) const;

      /// @brief Computes the lowest excitation energy at which one of the
      /// nuclear fragments considered by the HauserFeshbachDecay class may be
      /// emitted
      /// @param Zi Proton number of the initial nucleus
      /// @param Ai Mass number of the initial nucleus
      double unbound_threshold(const int Zi, const int Ai) const;

      /// @brief Computes the lowest excitation energy at which one of the
      /// nuclear fragments considered by the HauserFeshbachDecay class may be
      /// emitted
      /// @param initial_nucleus_pdg PDG code of the initial nucleus
      double unbound_threshold(const int initial_nucleus_pdg) const;

    protected:

      /// @brief Create the singleton MassTable object
      MassTable();

    private:

      // Helper function that converts the input JSON array into
      // (PDG code, mass) pairs and stores them in map_to_use
      void assign_masses(const marley::JSON& obj_array,
        const std::string& array_key,
        std::unordered_map<int, double>& map_to_use);

      // Function used internally by the mass table. Returns an atomic mass and
      // loads the boolean value exp with true if it is an experimental value
      // from the lookup table (and therefore has units of micro-amu) or false
      // if it is a theoretical value computed using the liquid drop model (and
      // has units of MeV).  This method of looking up masses is used to avoid
      // unnecessary unit conversions that can lead to losses of precision.
      double lookup_atomic_mass(int nucleus_pid, bool& exp,
        bool theory_ok = true) const;
      double lookup_atomic_mass(int Z, int A, bool& exp,
        bool theory_ok = true) const;

      // Lookup table for particle masses. Keys are PDG particle
      // ID numbers, values are masses in micro-amu.
      std::unordered_map<int, double> particle_masses_;

      // Lookup table for atomic masses. Keys are PDG particle
      // ID numbers for the nuclei, values are masses in micro-amu.
      std::unordered_map<int, double> atomic_masses_;

      // Factor to use when converting from micro-amu to MeV
      static constexpr double micro_amu_ = 0.000931494061; // MeV/uAMU

      /// Name of the JSON data file containing the masses
      static const std::string data_file_name_;
  };

}

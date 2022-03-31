/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
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
#include "marley/Fragment.hh"
#include "marley/FileManager.hh"
#include "marley/JSON.hh"
#include "marley/MassTable.hh"
#include "marley/StructureDatabase.hh"
#include "marley/marley_utils.hh"

// Initialize the static data file name
const std::string marley::MassTable::data_file_name_ = "mass_table.js";

marley::MassTable::MassTable() {

  // Instantiate the file manager and use it to find
  // the mass table data file
  const auto& fm = marley::FileManager::Instance();
  std::string full_mt_file_name
    = fm.find_file( data_file_name_ );

  if ( full_mt_file_name.empty() ) {
    throw marley::Error("Could not find the MARLEY mass table"
      " data file " + data_file_name_ + ". Please ensure that"
      " the folder containing it is on the MARLEY search path."
      " If needed, the folder can be appended to the MARLEY_SEARCH_PATH"
      " environment variable.");
  }

  MARLEY_LOG_INFO() << "Loading particle and atomic masses from "
    << full_mt_file_name;

  // Read in the mass table from a JSON data file
  auto json_table = marley::JSON::load_file( full_mt_file_name );

  // Store the JSON entries in the relevant unordered maps
  if ( !json_table.has_key("particle_masses") ) {
    throw marley::Error("Problem reading the mass table data file "
      + data_file_name_ + ". Missing \"particle_masses\" JSON array.");
  }
  const auto& pm_json = json_table.at("particle_masses");
  assign_masses( pm_json, "particle_masses", this->particle_masses_ );

  if ( !json_table.has_key("atomic_masses") ) {
    throw marley::Error("Problem reading the mass table data file "
      + data_file_name_ + ". Missing \"atomic_masses\" JSON array.");
  }
  const auto& am_json = json_table.at("atomic_masses");
  assign_masses( am_json, "atomic_masses", this->atomic_masses_ );

}

const marley::MassTable& marley::MassTable::Instance() {

  // Create the mass table using a static variable. This ensures
  // that the singleton instance is only created once.
  static std::unique_ptr<marley::MassTable>
    the_instance(new marley::MassTable());

  // Return a reference to the singleton instance
  return *the_instance;
}

marley::MassTable& marley::MassTable::NonConstInstance() {
  return const_cast< marley::MassTable& >( marley::MassTable::Instance() );
}

double marley::MassTable::liquid_drop_model_atomic_mass(int Z, int A) const {
  return liquid_drop_model_mass_excess(Z, A) + micro_amu_*1e6*A;
}

double marley::MassTable::get_particle_mass(int particle_id) const {
  int id = particle_id;
  // The lookup table only includes entries for particles (as opposed to
  // antiparticles), so flip the sign of the input particle id for the
  // lookup if it represents an antiparticle.
  if (id < 0) id *= -1;
  // Find the particle's mass in the lookup table, and convert its
  // value from micro-amu to MeV
  return micro_amu_ * particle_masses_.at(id);
}

double marley::MassTable::get_atomic_mass(int nucleus_pid, bool theory_ok) const
{
  bool exp;
  double mass = lookup_atomic_mass(nucleus_pid, exp,
    theory_ok);
  if (exp) mass *= micro_amu_;
  return mass;
}

double marley::MassTable::lookup_atomic_mass(int nucleus_pid, bool& exp,
  bool theory_ok) const
{
  // Find the atom's mass (in micro-amu) in the lookup table using its
  // nucleus's particle ID number. If it can't be found, either return a
  // theoretical mass using the liquid drop model or throw an error.
  // std::unordered_map<int, double>::iterator search
  auto search = atomic_masses_.find(nucleus_pid);

  // If the mass was found in the lookup table, return it and flag it as an
  // experimental value
  if (search != atomic_masses_.end()) {
    exp = true;
    return search->second;
  }
  // Otherwise, return a theoretical estimate using the liquid drop model or
  // throw an error depending on whether the user has indicated that using a
  // theoretical estimate is acceptable. For either case, set the experimental
  // flag to false.
  else {
    exp = false;

    int Z = marley_utils::get_particle_Z(nucleus_pid);
    int A = marley_utils::get_particle_A(nucleus_pid);

    if (theory_ok) {
      return liquid_drop_model_atomic_mass(Z, A);
    }
    else throw marley::Error(std::string("Entry for Z = ")
      + std::to_string(Z) + " and A = " + std::to_string(A)
      + " not found in the MARLEY atomic mass table.");
  }
}

double marley::MassTable::lookup_atomic_mass(int Z, int A, bool& exp,
  bool theory_ok) const
{
  int nucleus_pid = marley_utils::get_nucleus_pid(Z, A);

  auto search = atomic_masses_.find(nucleus_pid);

  // If the mass was found in the lookup table, return it and flag it as an
  // experimental value
  if (search != atomic_masses_.end()) {
    exp = true;
    return search->second;
  }
  // Otherwise, return a theoretical estimate using the liquid drop model or
  // throw an error depending on whether the user has indicated that using a
  // theoretical estimate is acceptable. For either case, set the experimental
  // flag to false.
  else {
    exp = false;

    if (theory_ok) {
      return liquid_drop_model_atomic_mass(Z, A);
    }
    else throw marley::Error(std::string("Entry for Z = ")
      + std::to_string(Z) + " and A = " + std::to_string(A)
      + " not found in the MARLEY atomic mass table.");
  }
}

double marley::MassTable::get_binding_energy(int Z, int A, bool theory_ok) const
{
  int N = A - Z;
  double m_hydrogen_1 = atomic_masses_.at(1000010010);
  double mn = particle_masses_.at(marley_utils::NEUTRON);

  bool exp;
  double mN = lookup_atomic_mass(Z, A, exp, theory_ok);

  // Experimental masses are given in micro-amu, while our liquid drop model
  // estimates are given in MeV, so adjust the calculation appropriately
  // depending on which unit we are using for mN.
  if (exp) return micro_amu_ * (Z*m_hydrogen_1 + N*mn - mN);
  else return micro_amu_ * (Z*m_hydrogen_1 + N*mn) - mN;
}

double marley::MassTable::get_mass_excess(int Z, int A, bool theory_ok) const
{
  bool exp;
  double mN = lookup_atomic_mass(Z, A, exp, theory_ok);
  if (exp) return micro_amu_*(mN - A*1e6);
  else return mN - micro_amu_*A*1e6;
}

double marley::MassTable::get_atomic_mass(int Z, int A, bool theory_ok) const {
  bool exp;
  double mass = lookup_atomic_mass(Z, A, exp, theory_ok);
  if (exp) mass *= micro_amu_;
  return mass;
}

double marley::MassTable::get_fragment_separation_energy(int nuc_pdg,
  int frag_pdg, bool theory_ok) const
{
  int Zi = marley_utils::get_particle_Z( nuc_pdg );
  int Ai = marley_utils::get_particle_A( nuc_pdg );
  return this->get_fragment_separation_energy( Zi, Ai, frag_pdg, theory_ok );
}

double marley::MassTable::get_fragment_separation_energy(int Z, int A, int pid,
  bool theory_ok) const
{
  int Zx = marley_utils::get_particle_Z(pid);
  int Zf = Z - Zx;
  int Af = A - marley_utils::get_particle_A(pid);

  double extra_mass = Zx*particle_masses_.at(marley_utils::ELECTRON)
    + particle_masses_.at(pid);

  bool exp_i, exp_f;
  double m_atom_initial = lookup_atomic_mass(Z, A, exp_i, theory_ok);
  double m_atom_final = lookup_atomic_mass(Zf, Af, exp_f, theory_ok);

  if (exp_i) {
    if (exp_f) return micro_amu_*(m_atom_final - m_atom_initial + extra_mass);
    else return m_atom_final + micro_amu_*(extra_mass - m_atom_initial);
  }
  else if (exp_f) return micro_amu_*(m_atom_final + extra_mass)
    - m_atom_initial;
  else return micro_amu_*extra_mass + m_atom_final - m_atom_initial;
}

/// @details The liquid drop model parameters used here are based on those
/// given in <a href="http://dx.doi.org/10.1016/j.nuclphysa.2008.06.005">A. J.
/// Koning, et al., Nucl. Phys. A 810 (2008) pp. 13-76</a> for use with the
/// back-shifted Fermi gas nuclear level density model.
double marley::MassTable::liquid_drop_model_mass_excess(int Z, int A) const {

  // Liquid drop model parameters (taken from paper by A. J. Koning, et al.)
  static constexpr double Mn = 8.07144; // MeV
  static constexpr double MH = 7.28899; // MeV
  static constexpr double a1 = 15.677; // MeV
  static constexpr double a2 = 18.56; // MeV
  static constexpr double kappa = 1.79;
  static constexpr double c3 = 0.717; // MeV
  static constexpr double c4 = 1.21129; // MeV

  int N = A - Z;

  double kappa_term = kappa * std::pow((N - Z) / static_cast<double>(A), 2);
  double c1 = a1 * (1 - kappa_term);
  double c2 = a2 * (1 - kappa_term);

  double Evol = -c1 * A;
  double Esur = c2 * std::pow(A, 2.0/3.0);
  double Ecoul = (c3 / std::pow(A, 1.0/3.0) - c4 / A) * std::pow(Z, 2);

  double delta_LDM = 0;
  bool z_odd = Z % 2;
  bool n_odd = N % 2;

  // delta_LDM will be zero if the nucleus is odd-even
  if (z_odd && n_odd) delta_LDM = 11/std::sqrt(A);
  else if (!z_odd && !n_odd) delta_LDM = -11/std::sqrt(A);

  return Mn * N + MH * Z + Evol + Esur + Ecoul + delta_LDM;
}

double marley::MassTable::fragment_emission_threshold(const int Zi,
  const int Ai, const marley::Fragment& f) const
{
  // Separation energy for the fragment
  double Sa = this->get_fragment_separation_energy( Zi, Ai, f.get_pid() );
  return Sa;
}

// Returns the smallest nuclear fragment emission threshold (separation energy)
// for a particular initial nucleus. All nuclear fragments recognized by the
// StructureDatabase class are considered.
double marley::MassTable::unbound_threshold(const int Zi, const int Ai) const
{
  // Before looking up the separation energies, start by setting the unbound
  // threshold to infinity.
  double unbound_threshold = std::numeric_limits<double>::max();
  MARLEY_LOG_DEBUG() << "unbound_threshold = " << unbound_threshold << '\n';

  // Loop over each available nuclear fragment. If it has a smaller separation
  // energy than the current value of unbound_threshold, update the stored value
  for ( const auto& pair : marley::StructureDatabase::fragments() ) {
    const marley::Fragment& f = pair.second;
    double thresh = this->fragment_emission_threshold( Zi, Ai, f );
    MARLEY_LOG_DEBUG() << f.get_pid() << " emission threshold = "
      << thresh << '\n';
    if ( thresh < unbound_threshold ) unbound_threshold = thresh;
    MARLEY_LOG_DEBUG() << "unbound_threshold = " << unbound_threshold << '\n';
  }

  return unbound_threshold;
}

double marley::MassTable::unbound_threshold(const int initial_nucleus_pdg) const
{
  int Zi = marley_utils::get_particle_Z( initial_nucleus_pdg );
  int Ai = marley_utils::get_particle_A( initial_nucleus_pdg );
  return this->unbound_threshold( Zi, Ai );
}

void marley::MassTable::assign_masses(const marley::JSON& obj_array,
  const std::string& array_key, std::unordered_map<int, double>& map_to_use)
{
  if ( !obj_array.is_array() ) {
    throw marley::Error("The \"" + array_key + "\" key in the"
      " mass data file " + data_file_name_
      + " does not refer to a JSON array.");
  }
  auto elements = obj_array.array_range();
  if ( elements.begin() == elements.end() ) {
    throw marley::Error("The \"" + array_key + "\" array in"
      " the mass data file " + data_file_name_ + " is empty.");
  }
  bool pdg_ok, mass_ok;
  for (const auto& el : elements) {

    if ( !el.has_key("pdg") ) {
      throw marley::Error("Missing pdg code for an element of the"
      " \"" + array_key + "\" array in the mass data file " + data_file_name_);
    }
    else if ( !el.has_key("mass") ) {
      throw marley::Error("Missing mass for an element of the"
      " \"" + array_key + "\" array in the mass data file " + data_file_name_);
    }

    auto pdg_json = el.at("pdg");
    int pdg = pdg_json.to_long( pdg_ok );
    if ( !pdg_ok ) throw marley::Error(std::string("Invalid PDG code \"")
      + pdg_json.dump_string() + "\" given in the \"" + array_key + "\" array"
      " in the mass data file " + data_file_name_);

    auto mass_json = el.at("mass");
    double mass = mass_json.to_double( mass_ok );
    if ( !mass_ok ) throw marley::Error(std::string("Invalid mass value \"")
      + mass_json.dump_string() + "\" given in the \"" + array_key + "\" array"
      " in the mass data file " + data_file_name_);

    map_to_use[ pdg ] = mass;
  }

}

void marley::MassTable::set_particle_mass( int pdg_code, double mass ) {
  // Convert the input mass (assumed to be in MeV) into micro-amu to match
  // the convention for internal storage in the MassTable class
  this->particle_masses_[ pdg_code ] = mass / marley_utils::micro_amu;
  MARLEY_LOG_INFO() << "Set " << marley_utils::get_particle_symbol( pdg_code )
    << " mass to " << mass << " MeV";
}

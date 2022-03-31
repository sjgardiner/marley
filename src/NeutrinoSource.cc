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

#include <limits>

#include "marley/Generator.hh"
#include "marley/NeutrinoSource.hh"

marley::NeutrinoSource::NeutrinoSource(int particle_id) {
  if (!pdg_is_allowed(particle_id)) throw marley::Error(
    "Creating a neutrino source object that produces"
    " particles with PDG ID number " + std::to_string(particle_id)
    + " is not allowed.");
  else pid_ = particle_id;
}

/// PDG particle IDs of each neutrino that could possibly be produced by a
/// NeutrinoSource object
const std::set<int> marley::NeutrinoSource::pids_ = {
  marley_utils::ELECTRON_NEUTRINO,
  marley_utils::ELECTRON_ANTINEUTRINO,
  marley_utils::MUON_NEUTRINO,
  marley_utils::MUON_ANTINEUTRINO,
  marley_utils::TAU_NEUTRINO,
  marley_utils::TAU_ANTINEUTRINO,
  marley_utils::FERMIONIC_DM
};

double marley::NeutrinoSource::sample_incident_neutrino(int& pdg,
  marley::Generator& gen) const
{
  static double max = marley_utils::UNKNOWN_MAX;
  pdg = pid_;
  return gen.rejection_sample([this](double E)
    -> double { return this->pdf(E); }, get_Emin(), get_Emax(), max);
}

marley::FermiDiracNeutrinoSource::FermiDiracNeutrinoSource(int particle_id,
  double Emin, double Emax, double temp, double eta)
  : NeutrinoSource(particle_id), Emin_(Emin), Emax_(Emax), temperature_(temp),
  eta_(eta), C_(1.)
{
  // Normalize the source spectrum (not strictly necessary, but having the
  // spectrum approximately normalized makes the default rejection sampling
  // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
  double integral = marley_utils::num_integrate(
    [this](double E) -> double { return this->pdf(E); }, Emin_, Emax_);

  // Update the normalization constant, thereby normalizing this object's
  // pdf in the process.
  C_ /= integral;
}

double marley::FermiDiracNeutrinoSource::pdf(double E) const {
  if (E < Emin_ || E > Emax_) return 0.;
  else return (C_ / std::pow(temperature_, 4)) * (std::pow(E, 2)
    / (1 + std::exp((E / temperature_) - eta_)));
}

marley::BetaFitNeutrinoSource::BetaFitNeutrinoSource(int particle_id,
  double Emin, double Emax, double Emean, double beta)
  : NeutrinoSource(particle_id), Emin_(Emin), Emax_(Emax), Emean_(Emean),
  beta_(beta), C_(1.)
{
  if ( beta_ < 1. ) {
    throw marley::Error("For a \"beta-fit\" neutrino source, a value of"
      " the fit parameter beta < 1 is unphysical. It causes the neutrino energy"
      " distribution to diverge near zero. Please choose another value and"
      " try again.");
  }

  // Normalize the source spectrum (not strictly necessary, but having the
  // spectrum approximately normalized makes the default rejection sampling
  // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
  double integral = marley_utils::num_integrate(
    [this](double E) -> double { return this->pdf(E); }, Emin_, Emax_);

  // Update the normalization constant, thereby normalizing this object's
  // pdf in the process.
  C_ /= integral;
}

double marley::BetaFitNeutrinoSource::pdf(double E) const {
  if (E < Emin_ || E > Emax_) return 0.;
  // Guard against NaNs in the std::pow factor below
  else if ( E == 0. && beta_ < 1. ) return 0.;
  else return C_ * std::pow(E / Emean_, beta_ - 1.)
    * std::exp(-beta_ * E / Emean_);
}

marley::FunctionNeutrinoSource::FunctionNeutrinoSource(int particle_id,
  double Emin, double Emax, std::function<double(double)> prob_dens_func)
  : NeutrinoSource(particle_id), Emin_(Emin), Emax_(Emax)
{
  // Normalize the supplied spectrum (not strictly necessary, but having the
  // spectrum approximately normalized makes the default rejection sampling
  // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
  double integral = marley_utils::num_integrate(prob_dens_func, Emin, Emax);
  probability_density_ = [prob_dens_func, integral](double E)
    -> double { return prob_dens_func(E) / integral; };
}

marley::DecayAtRestNeutrinoSource::DecayAtRestNeutrinoSource(int particle_id)
  : NeutrinoSource(particle_id)
{
  int abs_pid = std::abs(particle_id);
  if ( abs_pid != marley_utils::ELECTRON_NEUTRINO &&
    abs_pid != marley_utils::MUON_NEUTRINO )
  {
    throw marley::Error("Invalid projectile "
      + marley_utils::get_particle_symbol(particle_id) + " requested"
      " for a muon decay-at-rest neutrino source");
  }
}

double marley::DecayAtRestNeutrinoSource::pdf(double E) const {
  if (E < Emin_ || E > Emax_) return 0.;
  // Note that both of these source spectra are normalized to 1
  // on the energy interval [0., m_mu_ / 2.]
  /// @todo Refine the approximate decay-at-rest Michel spectra to use more
  /// exact expressions.
  int abs_pid = std::abs(pid_);
  if (abs_pid == marley_utils::ELECTRON_NEUTRINO)
    return 96. * std::pow(E, 2) * m_mu_to_the_minus_four_
      * (m_mu_ - 2*E);
  // Spectrum for muon antineutrinos
  else return 16. * std::pow(E, 2) * m_mu_to_the_minus_four_
    * (3*m_mu_ - 4*E);
}

/// @brief Method called near the end of construction to verify that the
/// newly-created grid source object is valid.
void marley::GridNeutrinoSource::check_for_errors() {
  size_t grid_size = grid_.size();
  if (grid_size < 2) throw marley::Error(std::string("Grid with")
    + " less than 2 gridpoints passed to the constructor of"
    + " marley::GridNeutrinoSource.");

  double sum_of_PDs = 0.;
  for (size_t j = 0; j < grid_size; ++j) {
    auto& pair = grid_.at(j);
    if (pair.first < 0.) throw marley::Error(std::string("All energy")
      + " values used in a marley::GridNeutrinoSource"
      + " object must be nonnegative");

    // Prevent actually sampling an energy value of zero by advancing to
    // the next representable double value.
    else if (pair.first == 0.) pair.first = std::nextafter(0.,
      marley_utils::infinity);
    if (pair.second < 0.) throw marley::Error(std::string("All PDF")
      + " values used in a marley::GridNeutrinoSource object must be"
      + " nonnegative");
    else sum_of_PDs += pair.second;
  }

  if (sum_of_PDs <= 0.) throw marley::Error(std::string("All probability")
    + " density grid point values are zero for the neutrino source");
}

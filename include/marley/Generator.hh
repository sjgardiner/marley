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
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

// MARLEY includes
#include "marley/GammaStrengthFunctionModel.hh"
#include "marley/NeutrinoSource.hh"
#include "marley/NuclearReaction.hh"
#include "marley/LevelDensityModel.hh"
#include "marley/OpticalModel.hh"
#include "marley/Parity.hh"
#include "marley/ProjectileDirectionRotator.hh"
#include "marley/RotationMatrix.hh"
#include "marley/StructureDatabase.hh"
#include "marley/Target.hh"
#include "marley/marley_utils.hh"

namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
}

namespace marley {

  class ChebyshevInterpolatingFunction;
  class JSON;
  class JSONConfig;

  /// @brief The MARLEY Event generator
  class Generator {

    // Allow the JSONConfig class to access the private
    // constructors of Generator
    friend JSONConfig;

    public:

      /// @brief Create a Generator using default settings
      Generator();

      /// @brief Enumeration of the possible methods for handling CRPA strenght leakage
      /// below the unbound threshold
      enum class CRPADiscreteMode { IGNORE, MIRROR };

      /// @brief Create an Event using the NeutrinoSource, Target, Reaction,
      /// and StructureDatabase objects owned by this Generator
      std::shared_ptr< HepMC3::GenEvent > create_event();

      /// @brief Get the seed used to initialize this Generator
      inline uint_fast64_t get_seed() const;

      /// @brief Reseeds the Generator
      void reseed(uint_fast64_t seed);

      /// @brief Use a string to set this Generator's internal state
      /// @details This function is typically used to restore a Generator
      /// to a state saved using get_state_string().
      void seed_using_state_string(const std::string& state_string);

      /// @brief Get a string that represents the current internal state of
      /// this Generator
      std::string get_state_string() const;

      /// @brief Sample a random number uniformly on either [min, max) or
      /// [min, max]
      /// @param min Lower bound of the sampling interval
      /// @param max Upper bound of the sampling interval
      /// @param inclusive Whether the upper bound should be included (true)
      /// or excluded (false) from the possible sampling outcomes
      double uniform_random_double(double min, double max, bool inclusive);

      /// @brief Sample from a given 1D probability density function f(x) on
      /// the interval [xmin, xmax] using a simple rejection method
      /// @param f Probability density function to use for sampling
      /// @param xmin Lower bound of the sampling interval
      /// @param xmax Upper bound of the sampling interval
      /// @param max_search_tolerance Tolerance to use when finding the
      /// maximum of f(x) using <a href="http://tinyurl.com/ntqkfck">Brent's
      /// method</a>
      /// @return Sampled value of x
      double rejection_sample(const std::function<double(double)>& f,
        double xmin, double xmax, double& fmax, double safety_factor = 1.01,
        double max_search_tolerance = DEFAULT_REJECTION_SAMPLING_TOLERANCE_);

      /// @brief Sample from a given 1D cumulative density function cdf(x) on
      /// the interval [xmin, xmax] using bisection
      /// @param cdf Cumulative density function to use for sampling
      /// @param xmin Lower bound of the sampling interval
      /// @param xmax Upper bound of the sampling interval
      /// @return Sampled value of x
      double inverse_transform_sample(
        const marley::ChebyshevInterpolatingFunction& cdf,
        double xmin, double xmax, double bisection_tolerance = 1e-12);

      /// @brief Sample from a given 1D probability density function f(x) on
      /// the interval [xmin, xmax] using an inverse transform technique
      /// @param f Probability density function to use for sampling
      /// @param xmin Lower bound of the sampling interval
      /// @param xmax Upper bound of the sampling interval
      /// @return Sampled value of x
      double inverse_transform_sample(const std::function<double(double)>& f,
        double xmin, double xmax, double bisection_tolerance = 1e-12);

      /// @brief Get a reference to the StructureDatabase owned by this
      /// Generator
      marley::StructureDatabase& get_structure_db();

      /// @brief Get a const reference to the vector of Reaction objects
      /// owned by this Generator
      inline const std::vector< std::unique_ptr<marley::Reaction> >&
        get_reactions() const;

      /// @brief Take ownership of a new Reaction
      /// @param reaction A pointer to the new Reaction to use
      void add_reaction(std::unique_ptr<marley::Reaction> reaction);

      /// @brief Clear the vector of Reaction objects owned by this Generator
      void clear_reactions();

      /// @brief Sample a Reaction and an energy for the reacting neutrino
      /// @param[out] E Total energy of the neutrino undergoing the reaction
      /// @return Reference to the sampled Reaction owned by this Generator
      marley::Reaction& sample_reaction(double& E);

      /// @brief Probability density function that describes the distribution
      /// of reacting neutrino energies
      /// @details This function computes the cross-section weighted neutrino
      /// flux (normalized to unity between source_->E_min and
      /// source_->E_max) including cross-section contributions from all
      /// Reactions owned by this Generator. For the distribution of
      /// <i>incident</i> neutrino energies, use marley::NeutrinoSource::pdf()
      /// @param E Total energy of the reacting neutrino
      /// @return Probability density (MeV<sup> -1</sup>)
      double E_pdf(double E);

      /// @brief Get a const reference to the NeutrinoSource owned by this
      /// Generator
      /// @details Throws a marley::Error if this Generator does not own a
      /// NeutrinoSource object.
      const marley::NeutrinoSource& get_source() const;

      /// @brief Take ownership of a new NeutrinoSource, replacing any
      /// existing source owned by this Generator
      /// @param source A pointer to the new NeutrinoSource to use
      void set_source(std::unique_ptr<marley::NeutrinoSource> source);

      /// @brief Take ownership of a new Target, replacing any
      /// existing target owned by this Generator
      /// @param source A pointer to the new Target to use
      void set_target(std::unique_ptr<marley::Target> target);

      /// @brief Get a const reference to the Target owned by this
      /// Generator
      /// @details Throws a marley::Error if this Generator does not own a
      /// Target object.
      const marley::Target& get_target() const;

      /// @brief Sample from an arbitrary probability distribution (defined
      /// here as any object that implements an operator()(std::mt19937_64&)
      /// function)
      /// @detail This template function is based on
      /// https://stackoverflow.com/a/9154394/4081973
      template <class RandomNumberDistribution>
        inline auto sample_from_distribution(RandomNumberDistribution& rnd)
        -> decltype( std::declval<RandomNumberDistribution&>().operator()(
        std::declval<std::mt19937_64&>()) )
      {
        return rnd(rand_gen_);
      }

      /// @brief Sample from an arbitrary probability distribution (defined
      /// here as any object that implements an
      /// operator()(std::mt19937_64&, const ParamType&) function) using the
      /// parameters params
      template <class RandomNumberDistribution, typename ParamType>
        inline auto sample_from_distribution(RandomNumberDistribution& rnd,
        const ParamType& params) -> decltype(
        std::declval<RandomNumberDistribution&>().operator()(
        std::declval<std::mt19937_64&>(), std::declval<const ParamType&>() ) )
      {
        return rnd(rand_gen_, params);
      }

      /// @brief Sets the direction of the incident neutrinos to use when
      /// generating events
      /// @param dir_vec Vector that points in the direction of the incident
      /// neutrinos
      /// @note The dir_vec passed to this function does not need to be
      /// normalized, but it must have at least one nonzero element or a
      /// marley::Error will be thrown.
      void set_neutrino_direction(const std::array<double, 3>& dir_vec);

      /// @brief Gets the direction of the incident neutrinos that is used when
      /// generating events
      inline const std::array<double, 3>& neutrino_direction();

      /// @brief Sets the value of the weight_flux flag
      /// @details This is potentially dangerous. Use only if you know
      /// what you are doing.
      void set_weight_flux(bool should_we_weight );

      /// @brief Sets the value of the do_deexcitations flag.
      /// @details This should almost always be true. Use only if you
      /// know what you are doing.
      inline void set_do_deexcitations( bool do_them );

      /// @brief Returns the value of the crpa_discrete_mode flag.
      inline CRPADiscreteMode crpa_discrete_mode() const;

      /// @brief Sets the value of the crpa_discrete_mode flag.
      inline void set_crpa_discrete_mode( CRPADiscreteMode mode );

      /// @brief Computes the flux-averaged total cross section for all
      /// enabled neutrino reactions, taking target atom fractions into
      /// account as appropriate
      /// @details If flux weighting is disabled (via a call to
      /// set_weight_flux()) then this function will return zero
      /// @return Total cross section (MeV<sup> -2</sup>)
      double flux_averaged_total_xs() const;

      /// @brief Computes the total cross section at fixed energy for all
      /// configured reactions involving a particular target atom.
      /// @details Atom fractions in the owned Target are ignored by this
      /// function.
      /// @note This function is not used as part of the normal MARLEY
      /// workflow. It serves as part of an API that enables MARLEY to
      /// be interfaced with an external flux driver.
      /// @param pdg_a The PDG code for the projectile
      /// @param KEa The kinetic energy of the projectile (MeV)
      /// @param pdg_atom The nuclear PDG code for the atomic target
      /// @return Total cross section (MeV<sup> -2</sup> / atom)
      double total_xs(int pdg_a, double KEa, int pdg_atom) const;

      /// @brief Computes the abundance-weighted total cross section at fixed
      /// energy for all configured reactions
      /// @details Atom fractions in the owned Target are used to perform
      /// the weighting by nuclide abundance
      /// @note This function is not used as part of the normal MARLEY
      /// workflow. It exposes the abundance-weighted total cross section
      /// for use by the mardumpxs command-line tool
      /// (see examples/executables/mardumpxs.cc)
      /// @param pdg_a The PDG code for the projectile
      /// @param KEa The kinetic energy of the projectile (MeV)
      /// @return Abundance-weighted total cross section (MeV<sup> -2</sup> / atom)
      double total_xs(int pdg_a, double KEa) const;

      /// @brief Creates an event object for a fixed projectile species,
      /// kinetic energy, and atomic target
      /// @details If no energetically-accessible reaction is available for
      /// the given input parameters, then a marley::Error will be thrown.
      /// @note This function is not used as part of the normal MARLEY
      /// workflow. It serves as part of an API that enables MARLEY to
      /// be interfaced with an external flux driver.
      /// @param pdg_a The PDG code for the projectile
      /// @param KEa The kinetic energy of the projectile (MeV)
      /// @param pdg_atom The nuclear PDG code for the atomic target
      /// @param dir_vec Direction three-vector of the projectile
      std::shared_ptr< HepMC3::GenEvent > create_event( int pdg_a, double KEa,
        int pdg_atom, const std::array<double, 3>& dir_vec );

      /// @brief Provides access to the owned ProjectileDirectionRotator
      inline marley::ProjectileDirectionRotator& get_rotator()
        { return rotator_; }

      inline const std::string& json_config() const { return json_config_; }

      inline void set_run_info( const std::shared_ptr< HepMC3::GenRunInfo >&
        run_info ) { run_info_ = run_info; }

      /// Convert a string to a CRPA discrete mode value
      static CRPADiscreteMode crpa_discrete_mode_from_string(
        const std::string& str );

      /// Convert a CRPA discrete mode value to a string
      static std::string string_from_crpa_discrete_mode(
        CRPADiscreteMode mode );

    private:

      /// @brief Create a Generator using default settings except for
      /// a given initial seed
      /// @param seed The initial seed to use for this Generator
      Generator(uint_fast64_t seed);

      /// @brief Helper function that updates the normalization factor to
      /// use in E_pdf()
      void normalize_E_pdf();

      /// @brief Print the MARLEY logo (called once during construction of
      /// the first Generator object) to any active Logger streams
      void print_logo();

      /// @brief Seed for the random number generator
      uint_fast64_t seed_;

      /// @brief 64-bit Mersenne Twister random number generator
      std::mt19937_64 rand_gen_;

      /// @brief Default stopping tolerance for rejection sampling
      static constexpr double DEFAULT_REJECTION_SAMPLING_TOLERANCE_ = 1e-8;

      /// @brief Normalizaton factor for E_pdf()
      double norm_ = 1.;

      /// @brief NeutrinoSource used to sample reacting neutrino energies
      std::unique_ptr<marley::NeutrinoSource> source_;

      /// @brief Representation of the neutrino target material (may be
      /// composite)
      std::unique_ptr<marley::Target> target_;

      /// @brief StructureDatabase used to simulate nuclear de-excitations
      /// when creating Event objects
      std::unique_ptr<marley::StructureDatabase> structure_db_;

      /// @brief Reaction(s) used to sample reacting neutrino energies
      std::vector< std::unique_ptr<marley::Reaction> > reactions_;

      /// @brief Total cross section values to use as weights for sampling a
      /// Reaction in sample_reaction()
      /// @details These total cross section values are energy dependent. They
      /// are therefore updated with every call to E_pdf().
      std::vector<double> total_xs_values_;

      /// @brief Discrete distribution used for Reaction sampling
      std::discrete_distribution<size_t> r_index_dist_;

      /// @brief Whether the generator should weight the incident
      /// neutrino spectrum by the reaction cross section(s)
      /// @details Don't change this unless you understand what you
      /// are doing!
      bool weight_flux_ = true;

      /// @brief Current estimate of the maximum value of E_pdf(),
      /// the probability density used for sampling reacting neutrino
      /// energies.
      /// @details This variable is used to double-check the validity
      /// of the estimated maximum during rejection sampling
      double E_pdf_max_;

      /// @brief When E_pdf(), the probability density used for
      /// sampling reacting neutrino energies, is changed or
      /// initialized for the first time, E_pdf_max_ is set
      /// to this value.
      double E_PDF_MAX_DEFAULT_ = marley_utils::UNKNOWN_MAX;

      /// @brief Sets the default value of E_pdf_max_
      /// @details This function is used only in cases where the user
      /// has provided their own estimate of the PDF maximum (likely
      /// because automatic estimation failed)
      inline void set_default_E_pdf_max( double def_max ) {
        E_PDF_MAX_DEFAULT_ = def_max;
      }

      /// @brief Flag manipulated by JSONConfig to prevent premature
      /// normalization of E_pdf() during construction of a Generator
      /// object
      bool dont_normalize_E_pdf_ = false;

      /// @brief Helper object that rotates events as needed to achieve
      /// the configured projectile direction
      marley::ProjectileDirectionRotator rotator_;

      /// @brief Flag that controls whether de-excitations are applied
      /// after two-two scatters are simulated
      /// @details This should be set to true except under unusual
      /// circumstances
      bool do_deexcitations_ = true;

      /// @brief Helper map used by the methods to convert a CRPADiscreteMode value from a std::string
      static std::map< CRPADiscreteMode, std::string > crpa_discrete_mode_string_map_;

      /// @brief Flag that controls how we deal with the CRPA strength that
      /// leaks under the unbound threshold
      CRPADiscreteMode crpa_discrete_mode_;

      /// @brief Information describing the generator configuration at the
      /// start of a run. This is saved together with the events in HepMC3
      /// files.
      std::shared_ptr< HepMC3::GenRunInfo > run_info_;

      /// @brief String representation of the JSON configuration used to
      /// initialize the Generator. This will be saved to the output RunInfo
      /// to assist in resuming an interrupted job.
      std::string json_config_;

      /// @brief Computes the total cross section at fixed energy for all
      /// configured reactions involving a particular target atom.
      /// @details Atom fractions in the owned Target are ignored by this
      /// function.
      /// @note This function is not used as part of the normal MARLEY
      /// workflow. It serves as part of an API that enables MARLEY to
      /// be interfaced with an external flux driver.
      /// @note Null pointer values for index_vec and/or xsec_vec are valid
      /// and will disable storage of the relevant information
      /// @param pdg_a The PDG code for the projectile
      /// @param KEa The kinetic energy of the projectile (MeV)
      /// @param pdg_atom The nuclear PDG code for the atomic target
      /// @param index_vec Pointer to a vector that will be loaded
      /// with the indices of all reactions that can handle the
      /// input initial state
      /// @param xsec_vec Pointer to a vector that will be loaded
      /// with the total cross sections of all reactions that can handle the
      /// input initial state
      /// @return Total cross section (MeV<sup> -2</sup>)
      double total_xs(int pdg_a, double KEa, int pdg_atom,
        std::vector<size_t>* index_vec, std::vector<double>* xsec_vec) const;

      /// @brief Initializes the owned GenRunInfo object that will be used
      /// to associate run metadata with the output events
      void set_up_run_info();

      /// @brief Add final pieces of metadata (e.g., the RNG state) to an
      /// otherwise complete event
      void finish_event_metadata( HepMC3::GenEvent& ev );

      /// @brief Sets the saved string copy of the JSON configuration used to
      /// initialize the generator
      void set_json_config( const marley::JSON& jc );
  };

  // Inline function definitions
  inline uint_fast64_t Generator::get_seed() const { return seed_; }

  inline const std::vector<std::unique_ptr<marley::Reaction> >&
    Generator::get_reactions() const { return reactions_; }

  inline const std::array<double, 3>& Generator::neutrino_direction()
    { return rotator_.projectile_direction(); }

  inline void Generator::set_do_deexcitations( bool do_them )
    { do_deexcitations_ = do_them; }

  inline Generator::CRPADiscreteMode Generator::crpa_discrete_mode() const
    { return crpa_discrete_mode_; }

  inline void Generator::set_crpa_discrete_mode( CRPADiscreteMode mode )
    { crpa_discrete_mode_ = mode; }
}

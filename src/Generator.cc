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
#include <chrono>
#include <cmath>
#include <limits>
#include <string>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/ChebyshevInterpolatingFunction.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/Logger.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/Reaction.hh"
#include "marley/StructureDatabase.hh"

// The default constructor uses the system time as the seed and a
// default-constructured monoenergetic neutrino source. No reactions are
// defined, so the user must call add_reaction() at least once before using a
// default-constructed Generator.
marley::Generator::Generator()
  : seed_( std::chrono::system_clock::now().time_since_epoch().count() ),
  source_( new marley::MonoNeutrinoSource ),
  structure_db_( new marley::StructureDatabase )
{
  print_logo();
  reseed( seed_ );
}

// The seed-only constructor is like the default constructor, but it uses
// a specific initial seed.
marley::Generator::Generator( uint_fast64_t seed )
  : seed_( seed ), source_( new marley::MonoNeutrinoSource ),
  structure_db_( new marley::StructureDatabase )
{
  print_logo();
  reseed( seed_ );
}

// Helper map used to convert a CRPADiscreteMode value from a std::string
std::map< marley::Generator::CRPADiscreteMode, std::string > marley::Generator::crpa_discrete_mode_string_map_ =
{
  { marley::Generator::CRPADiscreteMode::IGNORE, "ignore" },
  { marley::Generator::CRPADiscreteMode::MIRROR, "mirror" },
  { marley::Generator::CRPADiscreteMode::ACCUMULATE, "accumulate" }
};

// Print the MARLEY logo to the logger stream(s) if you haven't already.
void marley::Generator::print_logo() {
  static bool printed_logo = false;
  if ( !printed_logo ) {
    MARLEY_LOG_INFO() << '\n' << marley_utils::marley_logo
      << "\nDon't worry about a thing,\n'Cause every little thing"
      << " gonna be all right.\n-- Bob, \"Three Little Birds\"\n\n"
      << "Model of Argon Reaction Low Energy Yields\n"
      << "version " << MARLEY_VERSION << '\n';
    printed_logo = true;
  }
}

std::shared_ptr< HepMC3::GenEvent > marley::Generator::create_event() {

  // (0) Initialize the run information if it has not been set up yet
  if ( !run_info_ ) {
    this->set_up_run_info();
  }

  // (1) Select a reacting neutrino energy and reaction using the
  // flux-weighted total cross section(s)
  double E_nu;
  marley::Reaction& r = sample_reaction( E_nu );

  // (2) Create the prompt two-two scattering event using the
  // sampled reaction object
  int pdg_a = source_->get_pid();
  std::shared_ptr< HepMC3::GenEvent > ev = r.create_event( pdg_a, E_nu, *this );

  // E.C.2 and E.C.3
  // Save total and reaction cross sections as metadata in the event.
  // NOTE: the expected value for NuHepMC3 E.C.2 is the total cross section for
  // the selected projectile and target, so we loop over relevant reactions
  // manually to avoid imposing target fraction weighting to this particular
  // calculation.
  double totXS = 0.;
  for ( const auto& temp_r : reactions_ ) {

    // Skip reactions which involve a different target atom from the selected
    // one
    if ( temp_r->atomic_target().pdg() != r.atomic_target().pdg() ) continue;

    // Skip reactions which involve a different projectile
    if ( temp_r->pdg_a() != r.pdg_a() ) continue;

    double temp_xsec = temp_r->total_xs( pdg_a, E_nu );
    if ( temp_xsec > 0. ) {
      totXS += temp_xsec;
    }

  }

  // Convert to picobarns
  totXS *= marley_utils::hbar_c2 * marley_utils::fm2_to_picobarn;

  ev->add_attribute( "TotXS",
    std::make_shared< HepMC3::DoubleAttribute >( totXS )
  );

  double procXS = r.total_xs( pdg_a, E_nu ) * marley_utils::hbar_c2
    * marley_utils::fm2_to_picobarn;

  ev->add_attribute( "ProcXS",
    std::make_shared< HepMC3::DoubleAttribute >( procXS )
  );

  // (3) If needed, de-excite the final-state residue
  if ( do_deexcitations_ ) {
    marley::NucleusDecayer nd;
    nd.process_event( *ev, *this );
  }

  // (4) If needed, rotate the event to match the desired projectile direction
  rotator_.process_event( *ev, *this );

  // (5) Finish adding metadata to the event object
  this->finish_event_metadata( *ev );

  // Return the completed event object
  return ev;
}

void marley::Generator::seed_using_state_string(
  const std::string& state_string)
{
  // TODO: add error handling here (check that state_string is valid)
  std::stringstream strstr( state_string );
  strstr >> rand_gen_;
}

void marley::Generator::reseed( uint_fast64_t seed ) {
  // This is an attempt to do a decent job of seeding the random number
  // generator, but optimally accomplishing this can be tricky (see, for
  // example, http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
  seed_ = seed;
  std::seed_seq seed_sequence{ seed_ };
  rand_gen_.seed( seed_sequence );

  MARLEY_LOG_INFO() << "Seeded random number generator with " << seed_;
}

std::string marley::Generator::get_state_string() const {
  std::stringstream ss;
  ss << rand_gen_;
  return ss.str();
}

void marley::Generator::normalize_E_pdf() {

  // If we're not ready to do the normalization, then just return without doing
  // anything. JSONConfig may set this flag to prevent premature calls to
  // normalize_E_pdf() as the Generator is being constructed
  if ( dont_normalize_E_pdf_ ) return;

  // This function is called whenever the reacting neutrino energy PDF changes,
  // so reset the estimated maximum PDF value to its default. We will update
  // this during rejection sampling.
  E_pdf_max_ = E_PDF_MAX_DEFAULT_;

  // Treat monoenergetic sources differently since they can cause
  // problems for the standard numerical integration check
  if ( source_->get_Emin() == source_->get_Emax() ) {
    // Set the normalization factor back to one. It's used
    // in the call to E_pdf() below, so we need to do this before
    // we assign it a different value.
    norm_ = 1.0; //
    // Now norm_ is assigned to be the product of the total cross section times
    // the source PDF at energy Emin
    norm_ = E_pdf( source_->get_Emin() );
    if ( norm_ <= 0. || std::isnan(norm_) ) {
      throw marley::Error("The total cross section for all defined reactions"
        " is <= 0 or NaN for the neutrino energy defined in a monoenergetic"
        " source. Please verify that your neutrino source produces particles"
        " above threshold for at least one reaction.");
    }
  }
  else {
    // Reset the normalization factor to its default of one until
    // we can calculate the new value. This prevents strange things
    // from happening when we lose precision due to an abnormally
    // high or low norm_ value from a previous source or reaction
    // definition.
    norm_ = 1.;

    // Update the normalization factor for use with the reacting neutrino
    // energy probability density function
    norm_ = marley_utils::num_integrate( [this](double E)
      -> double { return this->E_pdf(E); }, source_->get_Emin(),
      source_->get_Emax() );

    if ( norm_ <= 0. || std::isnan(norm_) ) {
      throw marley::Error( "The integral of the cross-section-weighted"
        " neutrino flux is <= 0 or NaN. Please verify that your neutrino"
        " source spectrum produces significant flux above the reaction"
        " threshold(s)." );
    }
  }
}

// Sample a random double uniformly between min and max using the class
// member random number generator rand_gen_. The inclusive flag
// determines whether or not max is included in the range. That is,
// when inclusive == false, the sampling is done on the interval [min, max),
// while inclusive == true uses [min, max].
double marley::Generator::uniform_random_double( double min, double max,
  bool inclusive )
{
  // Defaults to sampling from [0,1). We will always
  // explicitly supply the upper and lower bounds to
  // this distribution, so we won't worry about the
  // default setting.
  static std::uniform_real_distribution<double> udist;

  double max_to_use;

  if ( inclusive ) { // sample from [min, max]

    // Find the double value that comes immediately after max. This allows us
    // to sample uniformly on [min, max] rather than [min,max). This trick
    // comes from http://tinyurl.com/n3ocg3p.
    max_to_use = std::nextafter( max, std::numeric_limits<double>::max() );
  }
  else { // sample from [min, max)
    max_to_use = max;
  }

  std::uniform_real_distribution<double>::param_type params( min, max_to_use );

  // Sample a random double from this distribution
  return udist( rand_gen_, params );
}

/// @details The rejection method used by this function consists of the
/// following steps:
/// <ol><li>Find the maximum of the function f(x) on [xmin, xmax]</li>
/// <li>Sample an x value uniformly over f(x)'s domain</li>
/// <li>Sample a y value uniformly over [0, max(f(x))]</li>
/// <li>If y <= f(x), accept the sampled x value</li>
/// <li>If y > f(x), reject the sampled x value, and return
/// to step 2 to try again</li> </ol>
/// Note that f(x) does not need to be normalized, but its range must be
/// nonnegative. In the first step, an iterative method (<a
/// href="http://tinyurl.com/ntqkfck">Brent's method</a>) is
/// used to find the maximum of f(x). The iterations will continue until two
/// successive iterations agree within max_search_tolerance on the location of
/// the maximum of f(x). To avoid problems with functions that yield double
/// values that are small compared to a typical value of max_search_tolerance
/// (say, max_search_tolerance = 1e-8, while many neutrino cross sections of
/// interest for MARLEY are less than 1e-40 cm^2), MARLEY normalizes all
/// probability density functions to unity before using rejection sampling.
double marley::Generator::rejection_sample(
  const std::function<double(double)>& f, double xmin, double xmax,
  double& fmax, double safety_factor, double max_search_tolerance )
{
  // If we were passed the value marley_utils::UNKNOWN_MAX for fmax, then this
  // signals that we need to search for the function maximum ourselves.
  // Otherwise, we'll assume that the value passed over is good.
  if ( fmax == marley_utils::UNKNOWN_MAX  ) {
    // This variable will be loaded with the value of x
    // that corresponds to the maximum of f(x).
    // We don't actually use this, but currently it's
    // a required parameter of marley_utils::maximize
    double x_at_max;

    // Maximize the function and multiply by a safety factor just
    // in case we didn't quite find the exact peak
    fmax = marley_utils::maximize( f, xmin, xmax, max_search_tolerance,
      x_at_max ) * safety_factor;
  }

  double x, y, val;

  do {
    // Sample x value uniformly from [xmin, xmax]
    x = uniform_random_double( xmin, xmax, true );

    // Sample y uniformly from [0, fmax]
    y = uniform_random_double( 0., fmax, true );

    val = f( x );
    if ( val > fmax ) {

      MARLEY_LOG_WARNING() << "PDF value f(x) = "
      << val << " at x = " << x << " exceeded the estimated maximum"
      << " fmax = " << fmax << " during rejection sampling.";

      fmax = val * safety_factor;
      MARLEY_LOG_WARNING() << "A new estimate fmax = " << val * safety_factor
        << " will now be adopted.";
    }
  }
  // Keep sampling until you get a y value less than f(x)
  // (the probability density function evaluated at the sampled value of x)
  while ( y > val );

  return x;
}

double marley::Generator::E_pdf( double E ) {

  // Initialize the return value to zero
  double pdf = 0.;

  // Sum all of the reaction total cross sections, saving
  // each individual value along the way. Take weighting
  // by atom fraction in the target material into account.
  for ( size_t j = 0, s = reactions_.size(); j < s; ++j ) {

    // Get the current reaction
    const auto& react = reactions_.at( j );

    // Compute the total cross section for the current reaction for a single
    // target atom
    double tot_xs = react->total_xs( source_->get_pid(), E );

    // If the target_ member has not been initialized, don't bother doing any
    // weighting by atom fraction (equivalent to a weight of unity for all
    // target atoms)
    if ( target_ ) {
      // If it has been configured, then apply the appropriate atom fraction
      // weight from the target as appropriate.
      tot_xs *= target_->atom_fraction( react->atomic_target() );
    }

    // Cache the atom-fraction-weighted total cross section for sampling a
    // reaction mode later
    total_xs_values_.at( j ) = tot_xs;

    // Add the weighted total cross section value to the total
    pdf += tot_xs;
  }

  // Normally, we want to fold the flux with the reaction cross section(s)
  // in order to obtain the distribution of reacting neutrino energies
  if ( weight_flux_ ) {
    // Multiply the total cross section by the neutrino spectrum
    // from the source object to get the (unnormalized) PDF
    // for sampling reacting neutrino energies.
    pdf *= source_->pdf( E );
  }
  else {
    // If the user has specifically requested it, don't weight the
    // energy PDF by the cross section(s), as long as at least one of them
    // is non-vanishing
    if ( pdf <= 0. ) return 0.;
    pdf = source_->pdf( E );
  }

  //  Divide by the normalization factor (computed when this source
  //  was made available to the Generator) to obtain the normalized PDF.
  return pdf / norm_;
}

marley::Reaction& marley::Generator::sample_reaction( double& E ) {
  if ( reactions_.empty() ) throw marley::Error( "Cannot sample"
    " a reaction in marley::Generator::sample_reaction(). The vector of"
    " marley::Reaction objects owned by this generator is empty." );

  // Store the "old" value of E_pdf_max_, i.e., the one we had before calling
  // rejection_sample(). This will be used to check for problems.
  double old_max = E_pdf_max_;

  // TODO: protect against source_ changing E_min or E_max after you compute
  // the normalization factor norm_ in marley::Generator::init()
  E = rejection_sample( [this](double E_nu)
    -> double { return this->E_pdf(E_nu); }, source_->get_Emin(),
    source_->get_Emax(), E_pdf_max_ );

  // If the value of max changed after the call to rejection_sample() and the
  // old value wasn't UNKNOWN_MAX, then the rejection sampling routine must
  // have encountered a PDF value that was larger than our estimated maximum.
  // Alert the user about this and advise them to change the configuration
  // appropriately to avoid a biased reacting neutrino energy distribution.
  static bool issued_long_error_message = false;
  if ( old_max != marley_utils::UNKNOWN_MAX
    && old_max != E_pdf_max_ )
  {
    if ( !issued_long_error_message ) {
      MARLEY_LOG_ERROR() << "Estimation of the maximum PDF value failed when"
        << " using a rejection method to sample reacting neutrino energies.\n"
        << "This may occur when, e.g., an incident neutrino flux"
        << " is used that includes multiple sharp peaks.\n"
        << "To avoid biasing the energy distribution, please rerun the"
        << " simulation after adding the following line to the MARLEY job"
        << " configuration file:\n"
        << "    energy_pdf_max: " << E_pdf_max_ << ",\n"
        << "If this error message persists after raising energy_pdf_max to a"
        << " relatively high value, please contact the MARLEY developers for"
        << " troubleshooting help.";
      issued_long_error_message = true;
    }
    else {
      MARLEY_LOG_ERROR() << "The maximum PDF value for sampling reacting"
       << " neutrino energies was exceeded again. The new estimated maximum is"
       << "\n    energy_pdf_max: " << E_pdf_max_ << ',';
    }
  }

  // The atom-fraction-weighted total cross section values have already been
  // updated by the final call to E_pdf() during rejection sampling, so we can
  // now sample a reaction type using our discrete distribution object.
  std::discrete_distribution<size_t>::param_type
    params( total_xs_values_.begin(), total_xs_values_.end() );
  size_t r_index = r_index_dist_( rand_gen_, params );
  return *reactions_.at( r_index );
}

const marley::NeutrinoSource& marley::Generator::get_source() const {
  if ( source_ ) return *source_;
  else throw marley::Error( "Error in marley::Generator::get_source()."
    " The member variable source_ == nullptr." );
}

const marley::Target& marley::Generator::get_target() const {
  if ( target_ ) return *target_;
  else throw marley::Error( "Error in marley::Generator::get_target()."
    " The member variable target_ == nullptr." );
}

void marley::Generator::set_source(
  std::unique_ptr<marley::NeutrinoSource> source )
{
  // If we're passed a nullptr, then don't bother to do anything
  if ( source ) {
    // Transfer ownership of the neutrino source the generator, leaving
    // the std::unique_ptr passed to this function null afterwards.
    source_.reset( source.release() );

    // Don't bother to renormalize if there are no reactions defined yet
    if ( reactions_.empty() ) return;

    // Update the neutrino energy probability density function based on the
    // new source spectrum
    this->normalize_E_pdf();
  }
}

void marley::Generator::add_reaction(
  std::unique_ptr<marley::Reaction> reaction )
{
  // If we're passed a nullptr, then don't bother to do anything
  if ( reaction ) {

    // Transfer ownership to a new unique_ptr in the reactions vector, leaving
    // the original empty
    reactions_.push_back( std::move(reaction) );

    // Add a new entry in the reaction cross sections vector
    total_xs_values_.push_back( 0. );

    // TODO: consider adding a check to see whether source_ is non-null.
    // Right now, this shouldn't be possible, but an explicit check might
    // be good.

    // Update the neutrino energy probability density function by including the
    // cross section for the new reaction
    normalize_E_pdf();
  }
}

void marley::Generator::clear_reactions() {
  reactions_.clear();
  total_xs_values_.clear();
  // Reset the normalization factor to 1. We don't need it until we define
  // one or more new reactions.
  norm_ = 1.;
}

marley::StructureDatabase& marley::Generator::get_structure_db() {
  if ( structure_db_ ) return *structure_db_;
  else throw marley::Error( "Error in marley::Generator::get_structure_db()."
    " The member variable structure_db_ == nullptr." );
}

void marley::Generator::set_neutrino_direction(
  const std::array<double, 3>& dir_vec )
{
  rotator_.set_projectile_direction( dir_vec );

  const auto& normalized_dir_vec = rotator_.projectile_direction();

  // Print a log message announcing the change of direction
  std::string dir_msg("Incident neutrino direction: (");
  for ( size_t i = 0; i < 3; ++i ) {
    dir_msg += std::to_string( normalized_dir_vec[i] );
    if ( i < 2 ) dir_msg += ", ";
  }
  MARLEY_LOG_INFO() << dir_msg << ')';
}

void marley::Generator::set_weight_flux( bool should_we_weight ) {
  weight_flux_ = should_we_weight;
}

double marley::Generator::inverse_transform_sample(
  const std::function<double(double)>& f, double xmin, double xmax,
  double bisection_tolerance )
{
  // Build an approximate CDF corresponding to the integral of the input PDF.
  // Use a polynomial approximant at Chebyshev points to do it.
  /// @todo Remove hard-coded number of points here
  marley::ChebyshevInterpolatingFunction func(f, xmin, xmax,
    DEFAULT_N_CHEBYSHEV);
  auto cdf = func.cdf();

  // Now that we have a CDF to use for sampling, delegate the rest of the
  // action to the overloaded version of this function.
  return this->inverse_transform_sample(cdf, xmin, xmax, bisection_tolerance);
}

double marley::Generator::inverse_transform_sample(
  const marley::ChebyshevInterpolatingFunction& cdf, double xmin, double xmax,
  double bisection_tolerance )
{
  // Sample a probability value uniformly on [0, 1]
  double prob = uniform_random_double(0., 1., true);

  // If we chose an endpoint, we're done, so just return the appropriate one
  if ( prob == 0. ) return xmin;
  else if ( prob == 1. ) return xmax;

  // A properly normalized CDF should evaluate to unity at x = xmax. We enforce
  // this here so that the user doesn't have to do it in advance.
  double norm = cdf.evaluate( xmax );

  // Find the x value corresponding to the sampled probability via bisection
  // (slow but robust)
  double a = xmin;
  double b = xmax;
  while ( (b - a) > bisection_tolerance ) {
    double midpoint = (a + b) / 2.;
    double mid_cdf = cdf.evaluate( midpoint ) / norm;
    // If the CDF at the midpoint exactly matches our sampled
    // probability, we're done. Just return the midpoint.
    if ( mid_cdf == prob ) return midpoint;
    // Otherwise, shrink the bisection interval and try again
    else if ( mid_cdf > prob ) b = midpoint;
    else a = midpoint; // mid_cdf < prob
  }

  // Return the midpoint of the bisection interval as our sampled x value
  double x = (a + b) / 2.;
  return x;
}

double marley::Generator::flux_averaged_total_xs() const {
  // If we've disabled weighting the neutrino energy PDF
  // by the total cross section, just return zero
  if ( !weight_flux_ ) return 0.;

  double avg_total_xs = 0.;

  // For a monoenergetic source, don't bother to do the full
  // integral
  double Emin = source_->get_Emin();
  double Emax = source_->get_Emax();
  if ( Emin == Emax ) {
    avg_total_xs = norm_ / source_->pdf( Emin );
  }
  else {
    double source_norm = marley_utils::num_integrate(
      [this](double Ev) -> double { return this->source_->pdf(Ev); },
      Emin, Emax );

    // Use the precomputed integral of the reacting neutrino energy PDF
    avg_total_xs = norm_ / source_norm;
  }
  return avg_total_xs;
}

void marley::Generator::set_target( std::unique_ptr<marley::Target> target )
{
  // If we're passed a nullptr, then don't bother to do anything
  if ( target ) {
    // Transfer ownership of the neutrino target to the generator, leaving
    // the std::unique_ptr passed to this function null afterwards.
    target_.reset( target.release() );

    // Don't bother to renormalize if there are no reactions defined yet
    if ( reactions_.empty() ) return;

    // Update the neutrino energy probability density function based on the
    // new target composition
    this->normalize_E_pdf();
  }
}

double marley::Generator::total_xs( int pdg_a, double KEa, int pdg_atom ) const
{
  return this->total_xs( pdg_a, KEa, pdg_atom, nullptr, nullptr );
}

double marley::Generator::total_xs( int pdg_a, double KEa, int pdg_atom,
  std::vector<size_t>* index_vec, std::vector<double>* xsec_vec ) const
{
  double xsec_sum = 0.;

  if ( index_vec ) index_vec->clear();
  if ( xsec_vec ) xsec_vec->clear();

  for ( size_t j = 0u; j < reactions_.size(); ++j ) {

    const auto& r = reactions_.at( j );

    // Skip reactions which involve a different target atom
    if ( pdg_atom != r->atomic_target().pdg() ) continue;

    // Skip reactions which involve a different projectile
    if ( pdg_a != r->pdg_a() ) continue;

    // If the cross section is non-vanishing, store information about it
    // (as appropriate) and add it to the sum
    double xsec = r->total_xs( pdg_a, KEa );
    if ( xsec > 0. ) {
      xsec_sum += xsec;
      if ( index_vec ) index_vec->push_back( j );
      if ( xsec_vec ) xsec_vec->push_back( xsec );
    }

  }

  return xsec_sum;
}


std::shared_ptr< HepMC3::GenEvent > marley::Generator::create_event(
  int pdg_a, double KEa, int pdg_atom, const std::array<double, 3>& dir_vec )
{
  // (0) Initialize the run information if it has not been set up yet
  if ( !run_info_ ) {
    this->set_up_run_info();
  }

  // (1) Sample a reaction mode from all configured reactions that can handle
  // the given initial-state parameters
  std::vector<size_t> indices;
  std::vector<double> xsecs;
  double tot_xsec = this->total_xs( pdg_a, KEa, pdg_atom, &indices, &xsecs );

  if ( xsecs.empty() || tot_xsec <= 0. ) throw marley::Error(
    "Cannot create an event for a projectile with kinetic energy = "
    + std::to_string(KEa) + " MeV and PDG code " + std::to_string(pdg_a)
    + " striking an atom with PDG code " + std::to_string(pdg_atom)
    + ". The total cross section vanishes for all configured reactions" );

  // The total cross section values and indices in the full reactions_ vector
  // have already been loaded into temporary vectors, so we can immediately use
  // those to sample a reaction using a discrete distribution.
  std::discrete_distribution< size_t > react_dist( xsecs.begin(), xsecs.end() );
  size_t sampled_index = react_dist( rand_gen_ );
  auto& r = reactions_.at( indices.at(sampled_index) );

  // (2) Create the prompt two-two scattering event using the sampled reaction
  // object
  std::shared_ptr< HepMC3::GenEvent > ev = r->create_event( pdg_a, KEa, *this );

  // E.C.2 and E.C.3
  // Save total and reaction cross sections as metadata in the event
  double totXS = tot_xsec * marley_utils::hbar_c2
    * marley_utils::fm2_to_picobarn;

  ev->add_attribute( "TotXS",
    std::make_shared< HepMC3::DoubleAttribute >( totXS )
  );

  double procXS = xsecs.at( sampled_index ) * marley_utils::hbar_c2
    * marley_utils::fm2_to_picobarn;

  ev->add_attribute( "ProcXS",
    std::make_shared< HepMC3::DoubleAttribute >( procXS )
  );

  // Do the usual post-processing

  // (3) If needed, de-excite the final-state residue
  if ( do_deexcitations_ ) {
    marley::NucleusDecayer nd;
    nd.process_event( *ev, *this );
  }

  // (4) If needed, rotate the event to match the desired projectile direction
  // Set the incident neutrino direction for this event
  static marley::ProjectileDirectionRotator my_rotator;
  my_rotator.set_projectile_direction( dir_vec );

  // Rotate the coordinate system of the event if needed
  my_rotator.process_event( *ev, *this );

  // (5) Finish adding metadata to the event object
  this->finish_event_metadata( *ev );

  // Return the completed event object
  return ev;
}

// Compute the abundance-weighted total reaction cross section at fixed
// projectile kinetic energy
double marley::Generator::total_xs( int pdg_a, double KEa ) const {

  // Initialize the return value to zero
  double tot_xsec = 0.;

  // Sum all of the reaction total cross sections. Take weighting by atom
  // fraction in the target material into account.
  for ( const auto& react : reactions_ ) {

    // Compute the total cross section for the current reaction for a single
    // target atom
    double xsec = react->total_xs( pdg_a, KEa );

    // If the target_ member has not been initialized, don't bother doing any
    // weighting by atom fraction (equivalent to a weight of unity for all
    // target atoms)
    if ( target_ ) {
      // If it has been configured, then apply the appropriate atom fraction
      // weight from the target as appropriate.
      xsec *= target_->atom_fraction( react->atomic_target() );
    }

    // Add the weighted total cross section value to the total
    tot_xsec += xsec;
  }

  return tot_xsec;
}

void marley::Generator::set_up_run_info() {

  // G.R.1
  run_info_ = std::make_shared< HepMC3::GenRunInfo >();

  // G.R.2
  run_info_->add_attribute( "NuHepMC.Version.Major",
    std::make_shared< HepMC3::IntAttribute >(
      marley_hepmc3::NUHEPMC_MAJOR_VERSION )
  );

  run_info_->add_attribute( "NuHepMC.Version.Minor",
    std::make_shared< HepMC3::IntAttribute >(
      marley_hepmc3::NUHEPMC_MINOR_VERSION )
  );

  run_info_->add_attribute( "NuHepMC.Version.Patch",
    std::make_shared< HepMC3::IntAttribute >(
      marley_hepmc3::NUHEPMC_PATCH_VERSION )
  );

  // G.R.3
  run_info_->tools().emplace_back(
    HepMC3::GenRunInfo::ToolInfo{ "MARLEY", MARLEY_VERSION,
      MARLEY_GIT_REVISION }
  );

  // G.R.4
  marley_hepmc3::prepare_process_metadata( *run_info_ );

  // G.R.5
  marley_hepmc3::prepare_vertex_status_metadata( *run_info_ );

  // G.R.6
  marley_hepmc3::prepare_particle_status_metadata( *run_info_ );

  // G.R.7
  run_info_->set_weight_names( { "CV" } );

  // G.R.8
  marley_hepmc3::prepare_non_standard_pdg_code_metadata( *run_info_ );

  // G.C.1, G.C.4, G.C.5, G.C.6
  double avg_xsec = this->flux_averaged_total_xs(); // MeV^{-2}
  marley_hepmc3::apply_nuhepmc_runinfo_conventions( *run_info_, avg_xsec );

  // Save the information needed to restore an interrupted MARLEY job
  run_info_->add_attribute( "MARLEY.RNGseed",
    std::make_shared< HepMC3::StringAttribute >( std::to_string(seed_) )
  );

  run_info_->add_attribute( "MARLEY.JSONconfig",
    std::make_shared< HepMC3::StringAttribute >( json_config_ )
  );

}

void marley::Generator::finish_event_metadata( HepMC3::GenEvent& ev ) {

  // Associate the owned run information with the event
  ev.set_run_info( run_info_ );

  // Add the generator state after the event was completed as a string
  // attribute. This allows resuming an interrupted job from where it left off.
  ev.add_attribute( "MARLEY.GeneratorState",
    std::make_shared< HepMC3::StringAttribute >( this->get_state_string() )
  );

  // E.R.5
  // TODO: revisit spatial position when MARLEY is interfaced with a
  // detector geometry simulation
  const std::vector< double > lab_pos = { 0., 0., 0. };
  ev.add_attribute( "LabPos",
    std::make_shared< HepMC3::VectorDoubleAttribute >( lab_pos )
  );

}

void marley::Generator::set_json_config( const marley::JSON& jc ) {
  json_config_ = jc.dump_string();
}

// Convert a string to a CRPA discrete mode value
marley::Generator::CRPADiscreteMode marley::Generator::crpa_discrete_mode_from_string(
  const std::string& str )
{
  for ( const auto& pair : crpa_discrete_mode_string_map_ ) {
    if ( str == pair.second ) return pair.first;
  }
  throw marley::Error( "The string \"" + str + "\" was not recognized"
    " as a valid CRPA discrete mode setting" );
}

// Convert a CRPA discrete mode value to a string
std::string marley::Generator::string_from_crpa_discrete_mode(  marley::Generator::CRPADiscreteMode mode ) 
{
  auto it = crpa_discrete_mode_string_map_.find( mode );
  if ( it != crpa_discrete_mode_string_map_.end() ) return it->second;
  else throw marley::Error( "Unrecognized CRPADiscreteMode value encountered in"
    " marley::Generator::string_from_crpa_discrete_mode " );
}
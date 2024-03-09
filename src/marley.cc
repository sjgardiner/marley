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
#include <csignal>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"

// ROOT includes
#ifdef USE_ROOT
  #include "TFile.h"
  #include "TInterpreter.h"
  #include "TParameter.h"
  #include "TROOT.h"
  #include "TTree.h"
#endif

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Generator.hh"
#include "marley/Logger.hh"
#include "marley/JSONConfig.hh"
#include "marley/OutputFile.hh"

// Flag that will alert the main event generation loop that the user has
// interrupted the program's execution (probably by pressing ctrl+c). This is a
// global variable (out of necessity, since it must be modified by our signal
// handler), but its scope is limited to this source file by the use of the
// static keyword in its declaration.
volatile static std::sig_atomic_t interrupted = false;

// Function that will be used to handle a SIGINT signal in our event generation
// loop. The method used here is taken from the second answer at
// http://tinyurl.com/jm6zmwg
// Note that, since we don't use the signal code, we leave the integer argument
// of the function nameless. This prevents "unused parameter" warnings from the
// compiler.
void signal_handler( int )
{
  interrupted = true;
}

namespace {

  constexpr int DEFAULT_STATUS_UPDATE_INTERVAL = 100;

  // Show a number using one decimal digit without scientific notation. Used to
  // print certain numbers in this way without affecting the settings currently
  // in use for std::cout.
  std::string format_number( double number ) {
    static std::stringstream temp_stream;
    static bool configured = false;
    if (!configured) {
      temp_stream << std::fixed << std::setprecision(1);
    }

    // Clear the std::stringstream of any of its previous contents This can be
    // done more elegantly by move-assigning a default-constructed
    // std::stringstream object like this:
    // temp_stream = std::stringstream();
    // However, g++ versions < 5 don't implement the move assignment operator
    // for std::stringstream
    // (see discussion at http://stackoverflow.com/q/27152263/4081973).
    // The approach below is a reasonable workaround.
    temp_stream.str("");
    temp_stream.clear();

    temp_stream << number;

    return temp_stream.str();
  }

  void print_help(std::string executable_name) {
    constexpr char help_message1[] = "Usage: ";
    constexpr char help_message2[] = " [OPTION...] CONFIG_FILE\n"
      "\n"
      "  -h, --help     Print this help message\n"
      "  -v, --version  Print version and exit\n";

    std::cout << help_message1 + executable_name + help_message2;
    std::cout << "\nMARLEY home page: <http://www.marleygen.org>\n";
    std::cout << "E-mail bug reports to: <support@marleygen.org>\n";
    exit(0);
  }

  void print_version() {
    std::cout << "MARLEY (Model of Argon Reaction Low Energy Yields) "
      << MARLEY_VERSION << '\n';

    std::cout << "Copyright (C) 2016-2023 Steven Gardiner\n";
    std::cout << "License: GNU GPL version 3 "
      << "<http://opensource.org/licenses/GPL-3.0>\n";
    std::cout << "This is free software: you are free to change and"
      << " redistribute it.\n";
    exit( 0 );
  }

  // Use this instead of std::put_time to allow this executable to be built
  // with g++ 4.9 (issue fixed in 5.0). See discussion here:
  // http://stackoverflow.com/a/14137287
  std::string put_time( std::tm* time, const char* format )
  {
    constexpr size_t TIME_STR_SIZE = 100;
    std::string time_str( TIME_STR_SIZE, ' ' );
    // Pass a pointer to the first character in the time string as the char*
    // argument of std::strftime(). This is well-defined behavior because the
    // storage of std::string is guaranteed by the C++11 standard to be
    // contiguous, and we have pre-allocated storage inside the string of the
    // necessary size using the constructor above.
    std::strftime( &time_str.front(), TIME_STR_SIZE, format, time );
    marley_utils::trim_right_inplace( time_str );
    return time_str;
  }

  // Formats the lines of the status display shown at the bottom of the
  // screen when this executable is running
  std::string makeStatusLines( long ev_count, long num_events,
    long num_old_events,
    std::chrono::system_clock::time_point start_time_point,
    const std::vector< std::shared_ptr<marley::OutputFile> >& output_files )
  {
    std::chrono::system_clock::time_point current_time_point
      = std::chrono::system_clock::now();

    // Compute the average event generation rate (events / s) up to this point
    auto elapsed_time = std::chrono::duration_cast<
      marley_utils::seconds<double> >( current_time_point - start_time_point );
    long events_since_start = ev_count - num_old_events;
    double elapsed_seconds = elapsed_time.count();
    double avg_event_rate = events_since_start / elapsed_seconds;

    std::ostringstream temp_oss;
    // Print a status message showing the current number of events
    double percent_complete = static_cast<double>( ev_count )
      / num_events * 100.;
    temp_oss << "\nEvent Count = " << ev_count << "/" << num_events
      << " (" << format_number( percent_complete ) << "% complete, "
      << format_number( avg_event_rate ) << " events / s)\033[K\n";

    // Print timing information
    temp_oss << "Elapsed time: "
      << marley_utils::elapsed_time_string(start_time_point,
      current_time_point) << " (Estimated total run time: ";

    marley_utils::seconds< double > estimated_total_time =
      ( current_time_point - start_time_point )
      * ( static_cast<double>(num_events - num_old_events)
      / (ev_count - num_old_events) );

    temp_oss << marley_utils::duration_to_string
      < marley_utils::seconds<double> >( estimated_total_time )
      << ")\033[K\n";

    for ( const auto& file : output_files ) {
      temp_oss << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string( file->bytes_written(), 2 )
        << "\033[K\n";
    }

    std::time_t estimated_end_time = std::chrono::system_clock::to_time_t(
      start_time_point + std::chrono::duration_cast
      < std::chrono::system_clock::duration >( estimated_total_time ) );

    temp_oss << "MARLEY is estimated to terminate on "
      << put_time( std::localtime(&estimated_end_time), "%c %Z" ) << '\n';

    // Move up an extra line in the terminal for each output file because we're
    // displaying the amount of data written to disk for each of them on a
    // separate line.
    for ( size_t i = 0; i < output_files.size(); ++i ) temp_oss << "\033[F";

    // Move up four lines
    temp_oss << "\033[F\033[F\033[F\033[F";

    return temp_oss.str();
  }

  // Inserts an updated copy of the status lines after each
  // newline is written to std::cout (e.g., by the marley::Logger
  // class). Based on https://stackoverflow.com/a/22043916
  class StatusInserter : public std::streambuf {
    public:
      StatusInserter( std::streambuf* dest, const long& ev_count,
        const long& num_events, const long& num_old_events,
        const std::chrono::system_clock::time_point& start_time_point,
        const std::vector< std::shared_ptr<marley::OutputFile> >& output_files )
        : std::streambuf(), myDest_( dest ), myIsAtStartOfLine_( true ),
        do_status_( true ), ev_count_( ev_count ), num_events_( num_events ),
        num_old_events_( num_old_events ),
        start_time_point_( start_time_point ), output_files_( output_files )
        {}

      inline void set_do_status( bool do_it ) { do_status_ = do_it; }

    protected:
      std::streambuf* myDest_;
      bool myIsAtStartOfLine_;
      bool do_status_;
      const long& ev_count_;
      const long& num_events_;
      const long& num_old_events_;
      const std::chrono::system_clock::time_point& start_time_point_;
      const std::vector< std::shared_ptr<marley::OutputFile> >& output_files_;

      int overflow( int ch ) override {
        int retval = 0;
        if ( ch != traits_type::eof() ) {
          if ( do_status_ && myIsAtStartOfLine_ ) {
            std::string status = makeStatusLines(ev_count_,
              num_events_, num_old_events_, start_time_point_, output_files_);
            myDest_->sputn( status.data(), status.size() );
          }
          myIsAtStartOfLine_ = ch == '\n';
          if ( myIsAtStartOfLine_ ) {
            std::string erase( "\033[K" );
            myDest_->sputn( erase.data(), erase.size() );
          }
          retval = myDest_->sputc( ch );
        }
        return retval;
      }
  };

}

int main( int argc, char* argv[] ) {

  std::streambuf* cout_default_buf = std::cout.rdbuf();
  std::streambuf* cerr_default_buf = std::cerr.rdbuf();

  // Disable automatic logging of marley::Error objects.
  // We will handle this manually below in the catch block.
  marley::Error::set_logging_status( false );

  try {

    // Initialize the logger. Use a default log level of INFO. This will
    // quickly be overwritten while parsing the JSON configuration file.
    marley::Logger::Instance().add_stream( std::cout,
      marley::Logger::LogLevel::INFO );

    // TODO: if you need command line parsing beyond the trivial stuff used
    // here, consider using a header-only option parser library like cxxopts
    // (https://github.com/jarro2783/cxxopts)
    std::string config_file_name;

    // If the user has not supplied any command-line
    // arguments, display the standard help message
    // and exit
    if ( argc <= 1 ) {
      print_help( argv[0] );
    }
    // The first command-line argument does not begin with a
    // hyphen, so assume that it is the configuration file name.
    else if ( std::string(argv[1]).substr(0, 1) != "-" ) {
      config_file_name = argv[1];
    }
    // The first command-line argument begins with a hyphen,
    // so treat it as an option and react accordingly.
    // All of the current options cannot be combined, so
    // just parse the first one and ignore the others.
    else {
      std::string option = argv[1];
      if ( option == "-h" || option == "--help" ) {
        print_help( argv[0] );
      }
      else if ( option == "-v" || option == "--version" ) {
        print_version();
      }
      else if ( option == "--marley" ) {
        std::cout << marley_utils::marley_pic;
        exit( 0 );
      }
      else {
       std::cout << argv[0] << ": unrecognized "
         << "command line option '" <<  option << "'\n";
       print_help( argv[0] );
      }
    }

    // Parse the objects from the JSON-based configuration file
    marley::JSON json = marley::JSON::load_file( config_file_name );

    // Process them to get the generator configuration
    marley::JSONConfig jc( json );

    // Get the time that the program was started
    std::chrono::system_clock::time_point start_time_point
      = std::chrono::system_clock::now();

    std::time_t start_time = std::chrono::system_clock::to_time_t(
      start_time_point );

    std::cout << "\nMARLEY started on "
      << put_time( std::localtime(&start_time), "%c %Z" ) << '\n';

    long num_old_events = 0;

    marley::JSON ex_set = json.get_object( "executable_settings" );

    // Desired number of events to be generated in this run. This
    // will be read back from the ROOT file and overwritten
    // if we're doing a continuation run.
    long num_events = ex_set.get_long( "events", 1e3 );

    // Update the status messages at the bottom of the screen
    // after this many events have been generated. The user may
    // set a non-default value from the job configuration file.
    int status_update_interval = DEFAULT_STATUS_UPDATE_INTERVAL;
    if ( ex_set.has_key("status_update_interval") ) {
      const auto& sui = ex_set.at( "status_update_interval" );

      bool ok;
      int sui_value = sui.to_long( ok );

      // Check for settings that are not positive integers
      if ( !ok || sui_value < 1 ) {
        throw marley::Error( "Invalid value " + sui.dump_string()
          + " given for the \"status_update_interval\" key in the"
          " job configuration file" );
      }
      else status_update_interval = sui_value;
    }

    std::vector< std::shared_ptr<marley::OutputFile> > output_files;

    if ( ex_set.has_key("output") ) {
      marley::JSON output_set = ex_set.at( "output" );
      if ( !output_set.is_array() ) throw marley::Error( "The"
        " \"output\" key in the executable settings must have a value that"
        " is a JSON array." );
      else for ( const auto& el : output_set.array_range() ) {
        output_files.push_back( marley::OutputFile::make_OutputFile(el) );
      }
    }
    else {
      // If the user didn't specify anything for the output key, then
      // by default write to a single ASCII-format file.
      std::string out_config_str = "{ format: \"ascii\","
        " file: \"events.hepmc3\", mode: \"overwrite\", force: false }";
      auto out_config = marley::JSON::load( out_config_str );

      output_files.push_back( marley::OutputFile::make_OutputFile(out_config) );
    }

    // This std::unique_ptr to a Generator object will be initialized below
    // with the proper configuration.
    std::unique_ptr< marley::Generator > gen;

    // Check if we need to resume from one of the files. If so, go ahead and
    // do it. If more than one file is requested for the resume, complain.
    bool need_to_resume = false;
    for ( auto& file : output_files ) {
      if ( file->mode_is_resume() ) {
        if ( need_to_resume ) throw marley::Error( "Only one file may be used"
          " to resume a previous run." );
        else {
          need_to_resume = true;
          bool resume_ok = file->resume( gen, num_old_events );
          if ( !resume_ok ) throw marley::Error( "Failed to resume previous run"
            " from the file \"" + file->name() + '\"' );
        }
      }
    }

    // If we didn't resume a run from any of the files, use the current
    // configuration
    if ( !need_to_resume ) gen = std::make_unique<marley::Generator>(
      jc.create_generator() );

    // Use the signal handler defined above to deal with
    // SIGINT signals (e.g., ctrl+c interruptions initiated
    // by the user). This will allow us to terminate the
    // loop gracefully, leaving a valid event tree file, etc.
    std::signal( SIGINT, signal_handler );

    // Generate all of the requested events. End the loop early
    // if the user interrupts execution (e.g., via ctrl+C)
    long ev_count = 1 + num_old_events;

    // Make std::cout use our "status inserter" std::streambuf
    // object so that the status lines get automatically updated
    // with every newline
    StatusInserter my_status_inserter( cout_default_buf, ev_count,
      num_events, num_old_events, start_time_point, output_files );
    std::cout.rdbuf( &my_status_inserter );
    std::cerr.rdbuf( &my_status_inserter );

    // Compute the flux-averaged total cross section for all
    // enabled reactions
    double avg_tot_xs = gen->flux_averaged_total_xs(); // MeV^(-2)

    // Update the start time to just before we begin the event loop.
    // This will help us get the best estimate for the remaining time
    // that the program will run.
    start_time_point = std::chrono::system_clock::now();
    start_time = std::chrono::system_clock::to_time_t( start_time_point );

    // Main generator event loop
    for (; ev_count <= num_events && !interrupted; ++ev_count) {

      // Create an event using the generator object
      auto event = gen->create_event();
      // E.R.2
      event->set_event_number( ev_count );

      for ( const auto& file : output_files ) {
        file->write_event( event.get() );
      }

      // Print status messages about simulation progress after every
      // status_update_interval events have been generated
      if ( (ev_count - num_old_events) % status_update_interval == 1
        || ev_count == num_events || status_update_interval == 1 )
      {
        // Temporarily disable the auto-printing of the status lines by
        // the StatusInserter class. This will avoid duplication of the
        // information.
        my_status_inserter.set_do_status( false );

        // Print a status message showing the current number of events
        std::cout << makeStatusLines( ev_count, num_events, num_old_events,
          start_time_point, output_files );

        // Re-enable the auto-printing of the status lines now that we've
        // printed them manually
        my_status_inserter.set_do_status( true );
      }
    }

    // Restore the default std::streambuf to std::cout
    std::cout.rdbuf( cout_default_buf );
    std::cerr.rdbuf( cerr_default_buf );

    std::cout << "\033[E";

    for (const auto& file : output_files) {
      std::cout << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string( file->bytes_written() )
        << "\033[K\n";
    }

    // Display the time that the program terminated
    std::chrono::system_clock::time_point end_time_point
      = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(
      end_time_point);

    if (!interrupted) std::cout << "MARLEY terminated normally on ";
    else std::cout << "MARLEY was interrupted by the user on ";
    std::cout << put_time(std::localtime(&end_time), "%c %Z")
      << "\033[K\033[E\033[K\n\033[K";

    return 0;
  }

  catch ( const std::exception& error ) {
    // Restore the old std::streambuf for std::cout
    std::cout.rdbuf( cout_default_buf );
    std::cerr.rdbuf( cerr_default_buf );

    // Flush the Logger, then rethrow the error for the system to handle.
    // This will probably result in std::terminate() being called.
    auto& log = marley::Logger::Instance();
    log.flush();

    // We've flushed the logger, so we can now print the exception's error
    // message without messing up the logger output.
    MARLEY_LOG_ERROR() << error.what();
    //throw error;
  }

  // We shouldn't ever get here unless there was an error
  return 1;
}

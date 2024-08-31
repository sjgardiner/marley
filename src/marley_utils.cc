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

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <string>
#include <vector>

#include "marley/marley_utils.hh"
#include "marley/Integrator.hh"
#include "marley/Error.hh"

// Strings to use for latex table output of ENSDF data
std::string marley_utils::latex_table_1 = "\\documentclass[12pt]{article}\n"
  "\n"
  "\\usepackage{amsmath}\n"
  "\\usepackage{booktabs}\n"
  "\\usepackage[justification=justified,\n"
  "width=\\textwidth,\n"
  "labelformat=empty]{caption}\n"
  "\\usepackage[top=1in, bottom=1in, left=0.25in, right=0.25in]{geometry}\n"
  "\\usepackage{isotope}\n"
  "\\usepackage{longtable}\n"
  "\\usepackage{multirow}\n"
  "\\usepackage{siunitx}\n"
  "\n"
  "\\newcommand{\\ExtraRowSpace}{1cm}\n"
  "\n"
  "\\begin{document}\n"
  "\\begin{center}\n"
  "\\begin{longtable}{\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "c\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "%  table-column-width = 2cm]\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "}\n"
  "\\caption";

std::string marley_utils::latex_table_2 = "\\toprule\n"
  "%{\\centering\\textbf{Level Energy (keV)}}\n"
  "%& {\\centering\\textbf{Spin-Parity}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}\n"
  "%& {\\centering\\textbf{Final Energy (keV)}} \\\\\n"
  "{\\multirow{3}{2cm}{\\centering\\textbf{Level Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{Spin-Parity}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}}\n"
  "& {\\multirow{3}{2.5cm}{\\centering\\textbf{Final Level Energy (keV)}}} \\\\\n"
  "& & & & \\\\\n"
  "& & & & \\\\\n"
  "\\midrule\n"
  "\\endfirsthead\n"
  "\\caption[]";

std::string marley_utils::latex_table_3 = "\\toprule\n"
  "%{\\centering\\textbf{Level Energy (keV)}}\n"
  "%& {\\centering\\textbf{Spin-Parity}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}\n"
  "%& {\\centering\\textbf{Final Energy (keV)}} \\\\\n"
  "{\\multirow{3}{2cm}{\\centering\\textbf{Level Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{Spin-Parity}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}}\n"
  "& {\\multirow{3}{2.5cm}{\\centering\\textbf{Final Level Energy (keV)}}} \\\\\n"
  "& & & & \\\\\n"
  "& & & & \\\\\n"
  "\\midrule\n"
  "\\noalign{\\vspace{-\\ExtraRowSpace}}\n"
  "\\endhead\n"
  "\\bottomrule \\multicolumn{4}{r}{\\textit{Continued on next page}} \\\\\n"
  "\\endfoot\n"
  "\\bottomrule\n"
  "\\endlastfoot\n"
  "% Table data\n"
  "\n";

std::string marley_utils::latex_table_4 = "\\end{longtable}\n"
  "\\end{center}\n"
  "\\end{document}";

// This implementation of the complex gamma function is based on the
// Lanczos approximation and its Python implementation given
// on Wikipedia (https://en.wikipedia.org/wiki/Lanczos_approximation)
// The C++ version given here is taken almost verbatim from
// http://bytes.com/topic/c/answers/576697-c-routine-complex-gamma-function
std::complex<double> marley_utils::gamma(std::complex<double> z)
{
  // Initialize some constants used in the algorithm. The "static
  // const" keywords ensure that these constants are initialized only
  // once (not reinitialized each time this function is called)
  static const int g=7;
  static const double pi =
  3.1415926535897932384626433832795028841972;
  static const double p[g+2] = {0.99999999999980993, 676.5203681218851,
  -1259.1392167224028, 771.32342877765313, -176.61502916214059,
  12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
  1.5056327351493116e-7};

  if (std::real(z) < 0.5) {
    return pi / (std::sin(pi*z)*gamma(1.0-z));
  }

  z -= 1.0;

  std::complex<double> x = p[0];

  for (int j = 1; j < g + 2; ++j) {
    x += p[j]/(z+std::complex<double>(j,0));
  }

  std::complex<double> t = z + (g + 0.5);

  return std::sqrt(2*pi) * std::pow(t, z + 0.5) * std::exp(-t) * x;
}

// This function is a modified version of a public-domain implementation of
// Brent's algorithm for minimizing a function. You can download the original
// source code from http://tinyurl.com/hhoy3ky.

// The return value of minimize is the minimum of the function f.
// The location where f takes its minimum is returned in the variable minLoc.
// Notation and implementation based on Chapter 5 of Richard Brent's book
// "Algorithms for Minimization Without Derivatives".
double marley_utils::minimize(const std::function<double(double)> f, // [in] objective function to minimize
  double leftEnd,     // [in] smaller value of bracketing interval
  double rightEnd,    // [in] larger value of bracketing interval
  double epsilon,     // [in] stopping tolerance
  double& minLoc)     // [out] location of minimum
{
    double d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
    static const double c = 0.5*(3.0 - std::sqrt(5.0));
    static const double SQRT_DBL_EPSILON = std::sqrt(DBL_EPSILON);

    double& a = leftEnd;
    double& b = rightEnd;
    double& x = minLoc;

    v = w = x = a + c*(b - a);
    d = e = 0.0;
    fv = fw = fx = f(x);

    // Check stopping criteria
    while (m = 0.5*(a + b),
      tol = SQRT_DBL_EPSILON*fabs(x) + epsilon,
      t2 = 2.0*tol,
      fabs(x - m) > t2 - 0.5*(b - a))
    {
        p = q = r = 0.0;
        if (fabs(e) > tol)
        {
            // fit parabola
            r = (x - w)*(fx - fv);
            q = (x - v)*(fx - fw);
            p = (x - v)*q - (x - w)*r;
            q = 2.0*(q - r);
            (q > 0.0) ? p = -p : q = -q;
            r = e; e = d;
        }
        if (fabs(p) < fabs(0.5*q*r) && p < q*(a - x) && p < q*(b - x))
        {
            // A parabolic interpolation step
            d = p/q;
            u = x + d;
            // f must not be evaluated too close to a or b
            if (u - a < t2 || b - u < t2)
                d = (x < m) ? tol : -tol;
        }
        else
        {
            // A golden section step
            e = (x < m) ? b : a;
            e -= x;
            d = c*e;
        }
        // f must not be evaluated too close to x
        if (fabs(d) >= tol)
            u = x + d;
        else if (d > 0.0)
            u = x + tol;
        else
            u = x - tol;
        fu = f(u);
        // Update a, b, v, w, and x
        if (fu <= fx)
        {
            (u < x) ? b = x : a = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        }
        else
        {
            (u < x) ? a = u : b = u;
            if (fu <= fw || w == x)
            {
                v = w; fv = fw;
                w = u; fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u; fv = fu;
            }
        }
    }
    return  fx;
}

// We can maximize a function using the same technique by minimizing its
// opposite
double marley_utils::maximize(const std::function<double(double)> f, // [in] objective function to maximize
  double leftEnd,     // [in] smaller value of bracketing interval
  double rightEnd,    // [in] larger value of bracketing interval
  double epsilon,     // [in] stopping tolerance
  double& maxLoc)     // [out] location of maximum
{
  double result = minimize([&f](double x) -> double { return -1.0*f(x); },
    leftEnd, rightEnd, epsilon, maxLoc);
  return -1.0*result;
}

// Numerically integrate a given function f (that takes a
// double argument to integrate over and returns a double)
// over the interval [a,b] using Clenshaw-Curtis quadrature
double marley_utils::num_integrate(const std::function<double(double)> &f,
  double a, double b)
{
  /// @todo remove the hard-coded number of Chebyshev points here in
  /// favor of adaptive integration.
  static marley::Integrator integrator(50);
  return integrator.num_integrate(f, a, b);
}

// Solves a quadratic equation of the form A*x^2 + B*x + C = 0
// while attempting to minimize errors due to the inherent limitations
// of floating-point arithmetic. The variables solPlus and solMinus
// are loaded with the two solutions x = (-B ± sqrt(B^2 - 4*A*C))/(2*A).
// Unsurprisingly, solPlus corresponds to the choice of the plus sign,
// while solMinus corresponds to the choice of the minus sign.
void marley_utils::solve_quadratic_equation(double A, double B,
  double C, double &solPlus, double &solMinus) {

  // Restructure the calculation to avoid some potentially bad cancellations
  // See, e.g., http://www.petebecker.com/js/js200010.html for details
  double c = C/A;
  double b = B/(2*A);

  double discr = b*b - c;

  // Find both solutions of the quadratic equation while avoiding
  // an extra subtraction (which can potentially lead to catastrophic
  // loss of precision) between -b and the square root of the
  // discriminant.
  if (b > 0) {

    solMinus = -b - marley_utils::real_sqrt(discr);
    solPlus = c/solMinus;
  }
  else {

    solPlus = -b + marley_utils::real_sqrt(discr);
    solMinus = c/solPlus;
  }

}

// Takes the square root of nonnegative arguments. Returns zero otherwise.
// This sqrt implementation is a quick fix for cases where roundoff errors
// give an argument to std::sqrt() that is slightly negative, causing it
// to return NaN.

// TODO: consider other ways of resolving this problem
double marley_utils::real_sqrt(double num) {
  if (num < 0) {
    return 0;
  }
  else {
    return std::sqrt(num);
  }
}

// Prompts the user with a question that requires a yes/no answer. Returns
// true if the answer was yes, or false if the answer was no.
bool marley_utils::prompt_yes_no(const std::string& message) {
  std::cout << '\n';
  std::string response;
  while (std::cout << message << " [y/n]? "
    && std::getline(std::cin, response)
    && !(response == "y" || response == "n" || response == "Y"
    || response == "N"));
  if (response == "y" || response == "Y") return true;
  else return false;
}

// For a given atomic number Z and mass number A, return a matching ENSDF nucid
std::string marley_utils::nuc_id(int Z, int A) {
  // Check to make sure Z and A have acceptable values
  if (Z < 0 || A < 1 || A > 999) throw marley::Error(
    std::string("The atomic number Z = ") + std::to_string(Z)
    + " and the mass number A = " + std::to_string(A)
    + " do not correspond to a valid ENSDF nucid.");

  // Create a three-character string representing the mass number
  std::string atomic_mass_number;
  if (A < 10) {
    atomic_mass_number = "  " + std::to_string(A);
  }
  else if (A < 100) {
    atomic_mass_number = " " + std::to_string(A);
  }
  else {
    atomic_mass_number = std::to_string(A);
  }

  // Get the element symbol as a string
  std::string symbol = element_symbols.at(Z);

  // Make the symbol completely uppercase
  to_uppercase_inplace(symbol);

  // If the symbol is only one character, pad the string so that
  // it is two characters long
  if (symbol.length() == 1) symbol += " ";

  return atomic_mass_number + symbol;
}

// Function that takes a number of bytes and returns a string
// representing the amount of memory in more readable units
std::string marley_utils::num_bytes_to_string(double bytes,
  unsigned precision)
{
  double divisor;
  std::string SI_prefix;

  if (bytes < 1e3) return std::to_string(static_cast<int>(bytes)) + " B";
  else if (bytes < 1e6) {
    divisor = 1e3;
    SI_prefix = "k";
  }
  else if (bytes < 1e9) {
    divisor = 1e6;
    SI_prefix = "M";
  }
  else if (bytes < 1e12) {
    divisor = 1e9;
    SI_prefix = "G";
  }
  else {
    divisor = 1e12;
    SI_prefix = "T";
  }

  std::ostringstream out;
  out << std::fixed;
  out.precision(precision);
  out << bytes/divisor << " " << SI_prefix << "B";
  return out.str();
}

// This function exploits the observation (given in the first answer at
// http://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code)
// that the difference of two std::chrono::system_clock::time_point objects
// can be assigned to a std::chrono::system_clock::duration.
//
// Function that takes two std::system_clock::time_point objects and returns
// a string (in the format days hours:minutes:seconds) representing the time
// between them
std::string marley_utils::elapsed_time_string(
  std::chrono::system_clock::time_point &start_time,
  std::chrono::system_clock::time_point &end_time)
{
  std::chrono::system_clock::duration time_elapsed
    = end_time - start_time;

  return marley_utils::duration_to_string
    <std::chrono::system_clock::duration>(time_elapsed);
}

// Efficiently read in an entire file as a std::string
// This function was taken from
// http://insanecoding.blogspot.in/2011/11/how-to-read-in-file-in-c.html
std::string marley_utils::get_file_contents(std::string filename) {

  // Read from the file using a stream. Ignore windows line ending changes
  std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);

  if (in) {
    std::string contents;

    // Determine the size of the file. Allocate sufficient memory
    // so that the string can hold the entire file without
    // resizing itself
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);

    // Slurp in the entire file
    in.read(&contents[0], contents.size());
    in.close();
    return contents;
  }
  throw marley::Error("Could not read from file " + filename);
}

// Advance to the next line of an ifstream that either matches (match == true)
// or does not match (match == false) a given regular expression
std::string marley_utils::get_next_line(std::ifstream &file_in,
  const std::regex &rx, bool match)
{

  std::string line;
  while (!file_in.eof() && file_in.good()) {
    // Get the next line of the file
    std::getline(file_in, line);

    // Check to see if the new line fulfills the search criteria
    if (std::regex_match(line, rx) == match) {
      // If it does, return it
      return line;
    }
    // If not, keep looking
  }

  // If the end of the file is encountered before a suitable
  // line is found, return an empty string
  return std::string("");
}

// Version of get_next_line that stores the number of lines
// scanned into num_lines
std::string marley_utils::get_next_line(std::ifstream &file_in,
  const std::regex &rx, bool match, int& num_lines)
{

  num_lines = 0;

  std::string line;
  while (!file_in.eof() && file_in.good()) {
    // Get the next line of the file
    std::getline(file_in, line);

    ++num_lines;

    // Check to see if the new line fulfills the search criteria
    if (std::regex_match(line, rx) == match) {
      // If it does, return it
      return line;
    }
    // If not, keep looking
  }

  // If the end of the file is encountered before a suitable
  // line is found, return an empty string
  return std::string("");
}

/// @brief Convert a string to a neutrino PDG code
/// @param str String to attempt to convert
/// @param[out] pdg PDG code of the requested neutrino
/// @return Whether the conversion was successful (true) or not (false)
bool marley_utils::string_to_neutrino_pdg(const std::string& str, int& pdg) {
  if (str == "ve") pdg = marley_utils::ELECTRON_NEUTRINO;
  else if (str == "vebar") pdg = marley_utils::ELECTRON_ANTINEUTRINO;
  else if (str == "vu") pdg = marley_utils::MUON_NEUTRINO;
  else if (str == "vubar") pdg = marley_utils::MUON_ANTINEUTRINO;
  else if (str == "vt") pdg = marley_utils::TAU_NEUTRINO;
  else if (str == "vtbar") pdg = marley_utils::TAU_ANTINEUTRINO;
  else {
    pdg = 0;
    return false;
  }
  return true;
}

std::string marley_utils::neutrino_pdg_to_string(int pdg) {
  if (pdg == marley_utils::ELECTRON_NEUTRINO)
    return std::string("ve");
  else if (pdg == marley_utils::ELECTRON_ANTINEUTRINO)
    return std::string("vebar");
  else if (pdg == marley_utils::MUON_NEUTRINO)
    return std::string("vu");
  else if (pdg == marley_utils::MUON_ANTINEUTRINO)
    return std::string("vubar");
  else if (pdg == marley_utils::TAU_NEUTRINO)
    return std::string("vt");
  else if (pdg == marley_utils::TAU_ANTINEUTRINO)
    return std::string("vtbar");
  else return std::string("?");
}

// Trim an ENSDF nucid string and make two-letter element symbols have a
// lowercase last letter. Currently, no checking is done to see if the
// string is a valid nucid.
std::string marley_utils::nucid_to_symbol(std::string nucid) {
  if (nucid.length() != 5) {
    // Split the string into "A" and "element name" pieces
    std::smatch m;
    std::regex_search(nucid, m, std::regex("[0-9]+"));
    std::string z_str = m.str();
    std::string e_str = m.suffix().str();

    // If the element name has more than one letter,
    // then make the last one lowercase.
    if (e_str.length() > 1)
      e_str.back() = tolower(e_str.back());
    return z_str + e_str;
  }
  // The nucid is 5 characters long, so getting our
  // desired format is a lot easier
  nucid.back() = tolower(nucid.back());
  return marley_utils::trim_copy(nucid);
}

// Assign a helicity value based on the PDG code and check its validity
int marley_utils::get_particle_helicity( const int pdg ) {

  // Unphysical, just used as a placeholder
  constexpr int DUMMY_HELICITY = 0;
  // Helicity value for antineutrinos (right-handed)
  constexpr int RIGHT_HANDED = 1;
  // Helicity value for neutrinos (left-handed)
  constexpr int LEFT_HANDED = -1;

  if ( pdg == marley_utils::ELECTRON_NEUTRINO
    || pdg == marley_utils::MUON_NEUTRINO
    || pdg == marley_utils::TAU_NEUTRINO )
  {
    // All Standard Model neutrinos are left-handed
    return LEFT_HANDED;
  }
  else if ( pdg == marley_utils::ELECTRON_ANTINEUTRINO
    || pdg == marley_utils::MUON_ANTINEUTRINO
    || pdg == marley_utils::TAU_ANTINEUTRINO )
  {
    // All Standard Model antineutrinos are right-handed
    return RIGHT_HANDED;
  }
  else throw marley::Error( "Handling of particles with PDG code = "
    + std::to_string(pdg) + " is unimplemented in marley_utils::"
    "get_particle_helicity()" );

  return DUMMY_HELICITY;
}

// Converts an ENSDF nucid to an atomic number. Currently, no checking is done
// to see if the string is a valid nucid.
int marley_utils::nucid_to_Z(std::string nucid) {
  // String that will be loaded with the element symbol
  std::string e_str;

  if (nucid.length() != 5) {
    // Extract the element name from the string
    std::smatch m;
    std::regex_search(nucid, m, std::regex("[0-9]+"));
    e_str = m.suffix().str();

    // If the element name has more than one letter,
    // then make the last one lowercase.
    if (e_str.length() > 1)
      e_str.back() = tolower(e_str.back());
  }
  else {
    // The nucid is 5 characters long, so getting our
    // desired format is a lot easier
    nucid.back() = tolower(nucid.back());
    e_str = nucid.substr(nucid.size() - 2); // Get the last two characters
    marley_utils::trim_right_inplace(e_str); // Trims the string if needed
  }
  return atomic_numbers.at(e_str);
}

  // Lookup table for element symbols (keys are atomic numbers Z,
  // values are symbols on the periodic table). The symbol "Nn" is
  // used for a neutron to match the ENSDF convention.
  const std::unordered_map<int, std::string> marley_utils::element_symbols = {
    { 0, "Nn"},
    { 1, "H" },
    { 2, "He" },
    { 3, "Li" },
    { 4, "Be" },
    { 5, "B" },
    { 6, "C" },
    { 7, "N" },
    { 8, "O" },
    { 9, "F" },
    { 10, "Ne" },
    { 11, "Na" },
    { 12, "Mg" },
    { 13, "Al" },
    { 14, "Si" },
    { 15, "P" },
    { 16, "S" },
    { 17, "Cl" },
    { 18, "Ar" },
    { 19, "K" },
    { 20, "Ca" },
    { 21, "Sc" },
    { 22, "Ti" },
    { 23, "V" },
    { 24, "Cr" },
    { 25, "Mn" },
    { 26, "Fe" },
    { 27, "Co" },
    { 28, "Ni" },
    { 29, "Cu" },
    { 30, "Zn" },
    { 31, "Ga" },
    { 32, "Ge" },
    { 33, "As" },
    { 34, "Se" },
    { 35, "Br" },
    { 36, "Kr" },
    { 37, "Rb" },
    { 38, "Sr" },
    { 39, "Y" },
    { 40, "Zr" },
    { 41, "Nb" },
    { 42, "Mo" },
    { 43, "Tc" },
    { 44, "Ru" },
    { 45, "Rh" },
    { 46, "Pd" },
    { 47, "Ag" },
    { 48, "Cd" },
    { 49, "In" },
    { 50, "Sn" },
    { 51, "Sb" },
    { 52, "Te" },
    { 53, "I" },
    { 54, "Xe" },
    { 55, "Cs" },
    { 56, "Ba" },
    { 57, "La" },
    { 58, "Ce" },
    { 59, "Pr" },
    { 60, "Nd" },
    { 61, "Pm" },
    { 62, "Sm" },
    { 63, "Eu" },
    { 64, "Gd" },
    { 65, "Tb" },
    { 66, "Dy" },
    { 67, "Ho" },
    { 68, "Er" },
    { 69, "Tm" },
    { 70, "Yb" },
    { 71, "Lu" },
    { 72, "Hf" },
    { 73, "Ta" },
    { 74, "W" },
    { 75, "Re" },
    { 76, "Os" },
    { 77, "Ir" },
    { 78, "Pt" },
    { 79, "Au" },
    { 80, "Hg" },
    { 81, "Tl" },
    { 82, "Pb" },
    { 83, "Bi" },
    { 84, "Po" },
    { 85, "At" },
    { 86, "Rn" },
    { 87, "Fr" },
    { 88, "Ra" },
    { 89, "Ac" },
    { 90, "Th" },
    { 91, "Pa" },
    { 92, "U" },
    { 93, "Np" },
    { 94, "Pu" },
    { 95, "Am" },
    { 96, "Cm" },
    { 97, "Bk" },
    { 98, "Cf" },
    { 99, "Es" },
    { 100, "Fm" },
    { 101, "Md" },
    { 102, "No" },
    { 103, "Lr" },
    { 104, "Rf" },
    { 105, "Db" },
    { 106, "Sg" },
    { 107, "Bh" },
    { 108, "Hs" },
    { 109, "Mt" },
    { 110, "Ds" },
    { 111, "Rg" },
    { 112, "Cn" },
    //{ 113, "Uut" },
    { 114, "Fl" },
    //{ 115, "Uup" },
    { 116, "Lv" },
    //{ 117, "Uus" },
    //{ 118, "Uuo" },
  };

  // TODO: consider other ways of adding support for fast reverse lookups
  // in the element_symbols map rather than reproducing its contents here
  const std::unordered_map<std::string, int> marley_utils::atomic_numbers = {
    {"Nn", 0    },
    {"H" , 1    },
    {"He", 2    },
    {"Li", 3    },
    {"Be", 4    },
    {"B" , 5    },
    {"C" , 6    },
    {"N" , 7    },
    {"O" , 8    },
    {"F" , 9    },
    {"Ne", 10   },
    {"Na", 11   },
    {"Mg", 12   },
    {"Al", 13   },
    {"Si", 14   },
    {"P" , 15   },
    {"S" , 16   },
    {"Cl", 17   },
    {"Ar", 18   },
    {"K" , 19   },
    {"Ca", 20   },
    {"Sc", 21   },
    {"Ti", 22   },
    {"V" , 23   },
    {"Cr", 24   },
    {"Mn", 25   },
    {"Fe", 26   },
    {"Co", 27   },
    {"Ni", 28   },
    {"Cu", 29   },
    {"Zn", 30   },
    {"Ga", 31   },
    {"Ge", 32   },
    {"As", 33   },
    {"Se", 34   },
    {"Br", 35   },
    {"Kr", 36   },
    {"Rb", 37   },
    {"Sr", 38   },
    {"Y" , 39   },
    {"Zr", 40   },
    {"Nb", 41   },
    {"Mo", 42   },
    {"Tc", 43   },
    {"Ru", 44   },
    {"Rh", 45   },
    {"Pd", 46   },
    {"Ag", 47   },
    {"Cd", 48   },
    {"In", 49   },
    {"Sn", 50   },
    {"Sb", 51   },
    {"Te", 52   },
    {"I" , 53   },
    {"Xe", 54   },
    {"Cs", 55   },
    {"Ba", 56   },
    {"La", 57   },
    {"Ce", 58   },
    {"Pr", 59   },
    {"Nd", 60   },
    {"Pm", 61   },
    {"Sm", 62   },
    {"Eu", 63   },
    {"Gd", 64   },
    {"Tb", 65   },
    {"Dy", 66   },
    {"Ho", 67   },
    {"Er", 68   },
    {"Tm", 69   },
    {"Yb", 70   },
    {"Lu", 71   },
    {"Hf", 72   },
    {"Ta", 73   },
    {"W" , 74   },
    {"Re", 75   },
    {"Os", 76   },
    {"Ir", 77   },
    {"Pt", 78   },
    {"Au", 79   },
    {"Hg", 80   },
    {"Tl", 81   },
    {"Pb", 82   },
    {"Bi", 83   },
    {"Po", 84   },
    {"At", 85   },
    {"Rn", 86   },
    {"Fr", 87   },
    {"Ra", 88   },
    {"Ac", 89   },
    {"Th", 90   },
    {"Pa", 91   },
    {"U" , 92   },
    {"Np", 93   },
    {"Pu", 94   },
    {"Am", 95   },
    {"Cm", 96   },
    {"Bk", 97   },
    {"Cf", 98   },
    {"Es", 99   },
    {"Fm", 100  },
    {"Md", 101  },
    {"No", 102  },
    {"Lr", 103  },
    {"Rf", 104  },
    {"Db", 105  },
    {"Sg", 106  },
    {"Bh", 107  },
    {"Hs", 108  },
    {"Mt", 109  },
    {"Ds", 110  },
    {"Rg", 111  },
    {"Cn", 112  },
    //{ "Uut", 113 },
    {"Fl", 114  },
    //{ "Uup", 115 },
    {"Lv", 116  },
    //{ "Uus", 117 },
    //{ "Uuo", 118 },
  };

  const std::string marley_utils::marley_logo = "╔╦╗╔═╗╦═╗╦  ╔═╗╦ ╦\n"
                                                "║║║╠═╣╠╦╝║  ║╣ ╚╦╝\n"
                                                "╩ ╩╩ ╩╩╚═╩═╝╚═╝ ╩\n";

  const std::string marley_utils::marley_pic =
    "ZZ77I\?I7777\?\?+I\?\?+=====~=~~~:+=+7I$I$Z$$O"
    "Z77Z$O\?\?=\?\?$Z$Z=$=~=~,:::~~:=~===+=+=II\n"
    "Z$7IIIIII+I=++=+==+~=~~~~~===\?ZIZ$OOZZOOI8$ZZ"
    "O8Z8O$ZOO7Z+:,:::,::~:::::~+~:~=~=+\n"
    "$77I$II\?+I\?+++=+=+====+\?=~=I$7O7$OOZ$8OO8DO"
    "O8888O8ZOZZ8Z7I~:::::::::::~~:~==~==~\n"
    "I7I$I7++\?+++=+=++===~~++I\?OO$I$ZOZDO8O888D88"
    "8OO8888888ZOO$\?I::~,::,::::~:~~=~~==\n"
    "IIII+\?\?I\?=\?+=+=====~~=$$Z$$Z8IOZOZOD888D$="
    ":::,,,,,:Z88888ZZ$=\?:,~::::~~:~,:~~~:~\n"
    "77I$\?\?\?+\?+\?=+=~=~~~=7Z87ZZ$OIZZOO8888DD\?"
    "::::,,::::,,::O88O8ZZ=:=:::,::::::::=~:=\n"
    "II+\?+\?\?++++======++7Z8I7OO$$DO$8Z888DZ\?~~:"
    ",,,,,,,,,:,:~:$O8ZOZ\?=~,:,,:,:::::~~~:\n"
    "$+I\?+\?+=+===+=+\?+IZ8$DDDIZO8OOD8D88D87+:::="
    "7ZI7~,,,,:::::~8O8ZZ7+~::,:,::,:,:::~\n"
    "7+\?+\?+\?+======\?$OO$IOD8OD88O8OD7D8DD8I:~77"
    "$7ZZ$$O+::,:::~~=8O78OI\?~:,:,::,,,::~:\n"
    "I\?+++=\?+=~++$$ZZODZ88ZDDOZZDOZD8DDDD7~~Z\?=~"
    "==~Z88O+~,::::~~=88OO8$+~=:,,::,,,:,:\n"
    "+++++==++I$$ZO8$$O8DNDNDOZ8O88O8DDD$+\?$$8O$88"
    "8O8O$7:~~:~,,~~O8888OO$~,:,,:::,:,:\n"
    "+\?++====\?\?$O88Z88ODD$ZOZ$$OOOOODDD\?~~::~=$"
    "Z~Z7887I~::$Z$O\?::I888ZO8$\?~:,,::~::,,\n"
    "\?+++====\?$88DDZ$88OZDDD8$O88O88D8\?:::,:,,:+"
    "8OO+7$~,,,O8OOOZ==8888888O\?=~::::::,,\n"
    "++=====+\?OODZ$$DD8DDD8ZZ$$ODO8DD=::~,,:::,,,,"
    ":,,,,:,,:OI~:ZOZ88D8$8OO7+:,::::~::\n"
    "+++=+=\?\?8D87DZ8DDDDZ888O88OO88ZD\?\?\?:::~::"
    ":~=I::,,,:~,,88D:,7Z8D88O88O+=:,,,:~::,\n"
    "+=+=++78Z7OZZ8DDDD8DDD8$+78$8D8D7I~:~I\?II7$I\?"
    "ZI$DD8I\?+==~88II8DD8O88OZ~:,,:::,:,\n"
    "++===+ID$ODDDDNDDDDD7=$O\?$88NDD8I\?++I\?+77\?"
    "~~~:=\?IOZ$OO7O8787$8D88$88ZO=~::,:::::\n"
    "=+===+I$8D8=8NDDND=OZZD$\?7ODO8D87=\?++II+Z++\?"
    "I\?Z7$Z$I$D8~::~Z7DDD87O$OZ++,,,::::,\n"
    "===~+\?8DZ8N8DNDNNDD+I=77I7D8DDDOZ++==\?I8\?\?"
    "$8ZO$7=\?$\?Z$==::,:8DDD8888O$$~,,:::::,\n"
    "====$8OOODN\?DDNDNNDI=+ZIIZD8DDD$7+I+==IZ8Z~~\?"
    "~O8$7\?=ZI+$:::~DDDD8O88$87~,::::::,\n"
    "==+8$8O88N8ODNNNDNDND7Z=ZDD8NDNO+7+=\?+\?~Z\?,"
    ":,,+,+DDZ7Z+7~~,~DDDD888O8$$\?:,,,:::,\n"
    "+~=+ZO88DD8DNNDNDDND78DI\?O8DDDDZ\?7I=\?==~~ZD"
    "DDDO,:,=D87+I==~=DDDD8DDOO$I=::,::::,\n"
    "+=+\?8I8DDNO8NNNNDNNDDD\?IDON8DDD$7$I\?7=~=~~~"
    ":~8DDD.~888O7I\?\?8DDD88OD8O$~:::,,:,:,\n"
    "=~\?8$O8DDDZDNDNNDNDNDD78DDNNDDDZI+\?\?7~=+:=="
    "77~~:ZO,DOOZZ7$ZDDDD8D8D888$::::,::,,\n"
    "+=\?8IDDNDDODNNNNNNNN$D$DDDNDDDDO$=7=+=~\?+==\?"
    "+7Z7I7$I\?Z$$ODDDDD88DDDD88I\?,,::,:,:\n"
    "=\?$OI8NDND8DNNDNNNDD7788DNDD8DDNZ$I+\?~~~==+="
    "$$$II$I~ZIOO8DDDDD8D8D8DDO7=::,:,,:,\n"
    "I$8888DDNDDDNN8NNNNN7Z88NNDDNNDDD8ZO7I=~~:~\?="
    "77II\?\?7O$ZODDDDD8DDD8O88Z\?:,:::,::,\n"
    "$DDZ8DDNND8DNDNNNNND7$8ODN8DNND888D8$$ZI$$\?I~"
    "=+I\?+7Z7$ODDDNDD88888OODO:::,,,,:::\n"
    "8DDOD8DDNN$ND7NNNNND8I88NN8DDNDDDZDD$ZZIZ88ZO7"
    "\?==78$IZDDDDDDDOD888DD8Z,::,,,,::,\n"
    "DD8DDDDDDNNNO\?NNNNN\?878NNNZOZDDD8ZDNDD888O8D"
    "OD87ZOI\?O8DDDDDDDODO$D7OO~:::,:,::::\n"
    "8NZOO8D8NDNNN7D88NNIOID8ND88$NDD8O$DDDDDN88Z8D"
    "8DD8Z7D8DDDDDDDOD8\?8OO=::,:::,:,,,\n"
    "DOZ8ZD8DDDNDDDN$8ONDNI88NNZD7DND$$8DNDNND8DNOO"
    "$NODODDDNDDDNDDD887IOO=:::::::,,,:\n"
    "NID8ND8DDNN+DNN$DID8D=ODDNZD\?8DN8ODZDDDN8NDOD"
    "8DDD8DDDDDDDDNDDD8~IDOO~:::,:,:,:,:\n"
    "IDODDDNDNDNZNNNI77D8+IZDN888=O8DD8DDNDNNNN8ODD"
    "DDN8DDDDDDDDDZDZI77+OO$I~:::::,:,:\n"
    "7NODDDD8NNDDN\?I7$INO++ODDI8D+7O8$88NNNDDNNDDD"
    "D8D8O$DDDDDND88O$$887$\?OZZ:::,,::,,\n"
    "NNDODDNDDNDNO$7I\?+NZ\?7DD$88=+\?ZZ$OIDDNNDND8"
    "DODND8$=DDDNDD8OOI$O87$O8Z8\?+:~::,,::\n"
    "NNNDZD$DN+DNOZZ\?==D7=IDN$DNI=I$ZZZ7$Z$DDNN8D8"
    "O$8$::DDNN8D8D8I8O8Z88DD8I\?+:::,,::\n"
    "NNNZ$DDNNZ7NOII+==D8+ZO=8D\?++I$$Z$88O8OZ8OZOI"
    "Z8$::NNDNDOD88$$ZIIIZOD88I8O+:::,,,\n"
    "NNNZ8DDND8ON7ZI==~DOZ88OO7$I\?\?7+I$\?$Z8$$Z$Z"
    "ZZD+::ODNZOD8=7I$D8Z$~I8DDD888+~~,,,,\n"
    "NNN$DD8NDDNN+7I==~7\?+7ZZ8I7I\?7\?77I7+ZZZ$$7Z"
    "NI~::~88O8DD+O$$DODZZZZDDODOO+:$\?:~,,\n"
    "NNDD8DNN8ZDN8I\?~~=~Z\?$I8D$7I=\?ZZ77\?\?$II$Z"
    "DZ=~::8D8ODZINZO$ODDI78=8DD8D+\?Z$$Z\?::,\n"
    "NN\?\?NNNDD8DD$\?7=NN=OZD8ODZ\?I=7ZI$I7O$$OO$="
    "DDI:~\?7OOIDD8\?$DDDO\?Z\?OZ87I8~OZ888Z+~:\n"
    "NNIO8NND88NN8IIND87ZDDO8$I\?7$ZIII7\?\?Z8:=\?+"
    "OZD:IIDIODND7\?O$NDOOI7::O$\?77DZO$D8ZI~\n"
    "NNI7DDNNDONNO\?O$ONDNI\?DN$+\?ZIZ$7I778=++==++"
    "ID\?Z$O$ZDD\?\?7$DDO\?$$O~:8+~OO8$\?888OZ=\n"
    "NNO78NNN78NNDDZ7DN$8DZ78$I7I+$I7\?\?Z~+++=++\?"
    "II$OOIIZDDD87ODD\?ZO8Z:8Z$7ODDZ78D8O$\?\n"
    "NN=ZONNDDODD77$ZN7D8D=NOO$$I$7I77$8+\?++\?=+\?"
    "$7I\?DO7DDDD87ND8Z8$8~:=OO7Z$OI+8DDZZ\?\n"
    "NDOO8NND$ON87NDD\?7\?$D7N8Z7II8ZI+Z8~\?==+==I\?"
    "DZID$=ZDD8ZIZNODOZ8D~:~ZD$\?\?\?D7O8DOI\?\n"
    "NNIODNNNDNNI$NO\?+DDNDD8ZZ$78OO$ID~\?+\?=++=IZ"
    "$OO8OI$N$$Z8NDZZO$D::~=$OZI+=$$$+O8Z+\n"
    "NI\?O8NNNDD\?ODD8=$$887DOZZ$I7\?\?\?$~+\?=+~=\?"
    "\?~O\?ZZN\?I7\?OZ8NDOOIDZI:~==+$OI~~~ZI~78ZI\n"
    "NZ7DNNNNN8DNNOD\?78N8N8$\?7I\?\?++7==+\?\?~~=+"
    "+Z~OZD7Z=DO7NND7$78DO=~====O$+:~=OZ=\?$O+\n"
    "N8OZDNNNNZNMN7$\?8OO\?OZ\?7$7III7$+++=====++=8"
    "+$ZI7\?$7$8D$D\?ODD8O==~=IO8\?\?+=D$I=7I=\n"
    "ND88ODN8NNNNN7\?Z7ZN=D8I$\?\?\?I+7~\?\?=++===+"
    "Z8\?+I$\?O\?DD8+8\?DOO8DN7I+==+=Z8=Z~$Z7+\?7\?\n";

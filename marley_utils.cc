#include <algorithm>
#include <cctype>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
#include "marley_utils.hh"

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

// Random number generator that will be used when selecting gammas for
// cascade simulations. Seed it using the system time
unsigned marley_utils::seed = std::chrono::system_clock::now().time_since_epoch().count();
std::knuth_b marley_utils::rand_gen(seed);


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

  for (int i = 1; i < g + 2; i++) {
    x += p[i]/(z+std::complex<double>(i,0));
  }

  std::complex<double> t = z + (g + 0.5);

  return std::sqrt(2*pi) * std::pow(t, z + 0.5) * std::exp(-t) * x;
}

// This function is a modified version of a public-domain implementation of
// Brent's algorithm for minimizing a function. You can download the original
// source code from http://www.codeproject.com/Articles/30201/Optimizing-a-Function-of-One-Variable

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
    static const double c = 0.5*(3.0 - sqrt(5.0));
    static const double SQRT_DBL_EPSILON = sqrt(DBL_EPSILON);
    
    double& a = leftEnd; double& b = rightEnd; double& x = minLoc;

    v = w = x = a + c*(b - a); d = e = 0.0;
    fv = fw = fx = f(x);

    // Check stopping criteria
    while (fabs(x - m) > t2 - 0.5*(b - a))
    {
        m = 0.5*(a + b);
        tol = SQRT_DBL_EPSILON*fabs(x) + epsilon; t2 = 2.0*tol;
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

// We can maximize a function using the same technique by minimizing its opposite
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
// over the interval [a,b] using the composite trapezoidal
// rule over n subintervals.
// (see http://en.wikipedia.org/wiki/Numerical_integration)
double marley_utils::num_integrate(const std::function<double(double)> &f,
  double a, double b, int n)
{
  double integral = 0;
  for(int k = 1; k < n-1; k++) {
    integral += ((b - a)/n)*f(a + k*(b - a)/n);
  }
  integral += ((b - a)/n)*(f(a)/2 + f(b)/2);
  return integral;
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

  // Find both solutions of the quadratic equation while avoiding
  // an extra subtraction (which can potentially lead to catastrophic
  // loss of precision) between -b and the square root of the
  // discriminant.
  if (b > 0) {
    solMinus = -b - std::sqrt(b*b - c);
    solPlus = c/solMinus;
  }
  else {
    solPlus = -b + std::sqrt(b*b - c);
    solMinus = c/solPlus;
  }

}

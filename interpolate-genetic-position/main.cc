/*!
  \file main.cc
  \brief main entry/exit for software. interprets command line arguments,
  dispatches tasks, exits \copyright Released under the MIT License. Copyright
  2023 Lightning Auriga
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "interpolate-genetic-position/cargs.h"
#include "interpolate-genetic-position/interpolator.h"

namespace igp = interpolate_genetic_position;

/*!
  \brief main program implementation
  @param argc number of command line entries, including program name
  @param argv array of command line entries
 */
int main(int argc, const char** const argv) {
  // parse command line input
  igp::cargs ap(argc, argv);
  // if help is requested or no flags specified
  if (ap.help() || argc == 1) {
    // print a help message and exit
    ap.print_help(std::cout);
    return 0;
  }
  if (ap.version()) {
    // print version information and exit
    ap.print_version(std::cout);
    return 0;
  }

  mpf_set_default_prec(ap.get_mpf_precision());

  std::string input = ap.get_input_filename();
  std::string preset = ap.get_input_preset();
  std::string genetic_map = ap.get_recombination_map();
  std::string map_format = ap.get_map_format();
  std::string output = ap.get_output_filename();
  std::string output_format = ap.get_output_format();
  bool verbose = ap.verbose();
  bool output_morgans = ap.output_morgans();
  double step_interval = ap.get_region_step_interval();
  unsigned fixed_output_width = ap.get_fixed_output_width();
  if (input.empty() && genetic_map.empty()) {
    throw std::runtime_error("only one of -i and -g can be read from stdin");
  }
  igp::check_io_combinations(preset, output_format);
  if (igp::string_to_format_type(preset) != igp::BED &&
      fabs(step_interval) > DBL_EPSILON) {
    std::cerr << "warning: step interval parameter is only respected "
              << "with bedfile (range) input" << std::endl;
    step_interval = 0.0;
  }

  igp::interpolator ip;
  ip.interpolate(input, preset, genetic_map, map_format, output, output_format,
                 output_morgans, step_interval, fixed_output_width, verbose);

  if (verbose) {
    std::cout << "all done woo!" << std::endl;
  }
  return 0;
}

/*!
  \file cargs.cc
  \brief method implementation for command line argument parser class
  \copyright Released under the MIT License.
  Copyright 2023 Lightning Auriga
*/

#include "interpolate-genetic-position/cargs.h"

void interpolate_genetic_position::cargs::initialize_options() {
  _desc.add_options()("help,h", "emit this help message")(
      "verbose,v", "emit extremely verbose debug logs")(
      "input,i",
      boost::program_options::value<std::string>()->default_value(""),
      "name of input variant/region query file (default: read from stdin)")(
      "preset,p", boost::program_options::value<std::string>(),
      "format of input file (accepted values: bim, map, snp, bed)")(
      "genetic-map,g",
      boost::program_options::value<std::string>()->default_value(""),
      "name of input genetic recombination map (default: read from stdin)")(
      "map-format,m", boost::program_options::value<std::string>(),
      "format of input recombination map (accepted values: bolt, bedgraph, "
      "bigwig)")(
      "output,o",
      boost::program_options::value<std::string>()->default_value(""),
      "name of output file (default: write to stdout)");
}

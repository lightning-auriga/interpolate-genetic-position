/*!
  \file cargs.cc
  \brief method implementation for command line argument parser class
  \copyright Released under the MIT License.
  Copyright 2023 Lightning Auriga
*/

#include "interpolate-genetic-position/cargs.h"

namespace igp = interpolate_genetic_position;

igp::cargs::cargs() {
  throw std::domain_error("cargs: do not use default constructor");
}

igp::cargs::cargs(int argc, const char **const argv)
    : _desc("Recognized options") {
  initialize_options();
  boost::program_options::store(
      boost::program_options::parse_command_line(argc, argv, _desc), _vm);
  boost::program_options::notify(_vm);
}

igp::cargs::cargs(const cargs &obj) : _desc(obj._desc), _vm(obj._vm) {}

igp::cargs::~cargs() throw() {}

void igp::cargs::initialize_options() {
  _desc.add_options()("help,h", "emit this help message")(
      "verbose,v", "emit extremely verbose debug logs")(
      "version", "emit program version and exit")(
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
      "name of output file (default: write to stdout)")(
      "output-format,f",
      boost::program_options::value<std::string>()->default_value(""),
      "format of output file (accepted values: bim, map, snp, bed)")(
      "output-morgans",
      "emit output genetic position in morgans instead of centimorgans")(
      "region-step-interval",
      boost::program_options::value<double>()->default_value(0.0),
      "genetic distance increment to add at boundaries of bed interval; "
      "experimental, generally should be kept at default");
}

bool igp::cargs::help() const { return compute_flag("help"); }

bool igp::cargs::verbose() const { return compute_flag("verbose"); }

bool igp::cargs::version() const { return compute_flag("version"); }

std::string igp::cargs::get_input_filename() const {
  return compute_parameter<std::string>("input");
}

std::string igp::cargs::get_input_preset() const {
  std::string preset = compute_parameter<std::string>("preset");
  if (preset.compare("bim") && preset.compare("map") && preset.compare("bed") &&
      preset.compare("snp") && preset.compare("vcf")) {
    throw std::runtime_error("invalid input preset format: \"" + preset + "\"");
  }
  return preset;
}

std::string igp::cargs::get_recombination_map() const {
  return compute_parameter<std::string>("genetic-map");
}

std::string igp::cargs::get_map_format() const {
  std::string map_format = compute_parameter<std::string>("map-format");
  if (map_format.compare("bolt") && map_format.compare("bedgraph") &&
      map_format.compare("bigwig")) {
    throw std::runtime_error("invalid genetic map format: \"" + map_format +
                             "\"");
  }
  return map_format;
}

std::string igp::cargs::get_output_filename() const {
  return compute_parameter<std::string>("output");
}

std::string igp::cargs::get_output_format() const {
  std::string output_format = compute_parameter<std::string>("output-format");
  if (output_format.compare("bim") && output_format.compare("map") &&
      output_format.compare("bed") && output_format.compare("snp")) {
    throw std::runtime_error("invalid output format: \"" + output_format +
                             "\"");
  }
  return output_format;
}

bool igp::cargs::output_morgans() const {
  return compute_flag("output-morgans");
}

double igp::cargs::get_region_step_interval() const {
  return compute_parameter<double>("region-step-interval");
}

bool igp::cargs::compute_flag(const std::string &tag) const {
  return _vm.count(tag);
}

void igp::cargs::print_help(std::ostream &out) {
  if (!(out << _desc))
    throw std::domain_error("cargs::print_help: unable to write to stream");
}

void igp::cargs::print_version(std::ostream &out) {
  if (!(out << INTERPOLATE_GENETIC_POSITION_PACKAGE_STRING << std::endl))
    throw std::domain_error("cargs::print_version: unable to write to stream");
}

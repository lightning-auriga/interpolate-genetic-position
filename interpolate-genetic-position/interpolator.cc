/*!
 \file interpolator.cc
 \brief implementation for interpolator class
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/interpolator.h"

namespace igp = interpolate_genetic_position;

igp::interpolator::interpolator() {}
igp::interpolator::~interpolator() throw() {}
void igp::interpolator::interpolate(
    const std::string &input_filename, const std::string &preset,
    const std::string &genetic_map_filename, const std::string &map_format,
    const std::string &output_filename, const std::string &output_format,
    bool output_morgans, const double &step_interval,
    unsigned fixed_output_width, bool verbose) const {
  input_variant_file input_variant_interface;
  input_genetic_map_file genetic_map_interface;
  output_variant_file output_variant_interface;
  input_variant_interface.set_fallback_stream(&std::cin);
  output_variant_interface.output_morgans(output_morgans);
  output_variant_interface.set_fixed_width(fixed_output_width);
  genetic_map_interface.set_fallback_stream(&std::cin);
  genetic_map gm(&genetic_map_interface);
  format_type map_ft = string_to_format_type(map_format);
  gm.open(genetic_map_filename, map_ft);
  query_file qf(&input_variant_interface, &output_variant_interface);
  format_type query_ft = string_to_format_type(preset);
  format_type output_ft = string_to_format_type(output_format);
  qf.open(input_filename, query_ft);
  qf.initialize_output(output_filename, output_ft);
  qf.set_step_interval(step_interval);
  mpf_class gpos_interpolated;
  std::vector<query_result> results;
  while (qf.get()) {
    gm.query(qf.get_chr(), qf.get_pos1(), qf.get_pos2(), verbose, &results);
    qf.report(results);
  }
  qf.close();
}

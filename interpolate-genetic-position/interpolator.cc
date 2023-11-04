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
void igp::interpolator::interpolate(const std::string &input_filename,
                                    const std::string &preset,
                                    const std::string &genetic_map_filename,
                                    const std::string &map_format,
                                    const std::string &output_filename,
                                    bool verbose) const {
  input_variant_file input_variant_interface;
  input_genetic_map_file genetic_map_interface;
  output_variant_file output_variant_interface;
  genetic_map gm(&genetic_map_interface);
  format_type map_ft = string_to_format_type(map_format);
  gm.open(genetic_map_filename, map_ft);
  query_file qf(&input_variant_interface, &output_variant_interface);
  format_type query_ft = string_to_format_type(preset);
  qf.open(input_filename, query_ft);
  qf.initialize_output(output_filename, query_ft);
  mpf_class gpos_interpolated;
  while (qf.get()) {
    gm.query(qf.get_chr(), qf.get_pos(), verbose, &gpos_interpolated);
    qf.report(gpos_interpolated);
  }
  qf.close();
}

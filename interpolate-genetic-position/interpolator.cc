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
  genetic_map gm;
  format_type map_ft = string_to_format_type(map_format);
  gm.open(genetic_map_filename, map_ft);
  query_file qf;
  format_type query_ft = string_to_format_type(preset);
  qf.open(input_filename, query_ft);
  mpf_class gpos_interpolated;
  std::ofstream output;
  if (!output_filename.empty()) {
    output.open(output_filename.c_str());
    if (!output.is_open()) {
      throw std::runtime_error(
          "interpolator::interpolate: cannot write "
          "output file \"" +
          output_filename + "\"");
    }
  }
  while (qf.get()) {
    gm.query(qf.get_chr(), qf.get_pos(), verbose, &gpos_interpolated);
    qf.report(gpos_interpolated,
              output_filename.empty() ? &std::cout : &output);
  }
  if (output.is_open()) {
    output.close();
  }
}

/*!
 \file interpolator.h
 \brief primary controller for interpolation, refactored out
 for compatibility with testing.
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_INTERPOLATOR_H_
#define INTERPOLATE_GENETIC_POSITION_INTERPOLATOR_H_

#include <gmpxx.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "interpolate-genetic-position/genetic_map.h"
#include "interpolate-genetic-position/input_genetic_map_file.h"
#include "interpolate-genetic-position/input_variant_file.h"
#include "interpolate-genetic-position/query_file.h"
#include "interpolate-genetic-position/utilities.h"

namespace interpolate_genetic_position {
/*!
 * \class interpolator
 * \brief primary interpolation controller. needs to be its own
 * class outside of main for easy compatibility with unit testing.
 */
class interpolator {
 public:
  /*!
   * \brief default constructor
   */
  interpolator();
  /*!
   * \brief destructor
   */
  ~interpolator() throw();
  /*!
   * \brief run interpolation on input and create output
   * \param input_filename name of input variant query file
   * \param preset descriptor of input query file format
   * \param genetic_map name of input recombination map file
   * \param map_format descriptor of genetic map format
   * \param output_filename name of output results file
   * \param output_format descriptor of output results file
   * \param output_morgans whether genetic position should be output
   * in morgans, as opposed to the default centimorgans
   * \param step_interval fixed genetic distance to add to boundary
   * between successive query regions in an input bedfile. interpolation
   * within regions behaves linearly regardless. this functionality,
   * when set to 0, is effectively disabled. this parameter is ignored
   * for all input types other than bedfiles.
   * \param fixed_output_width number of digits to print after
   * decimal for output floating point values, or 0 for the number
   * of digits to be scaled dynamically based on the number of sigfigs
   * \param verbose whether to emit (extremely) verbose logging to std::cout
   * can be empty string, in which case output is std::cout
   */
  void interpolate(const std::string &input_filename, const std::string &preset,
                   const std::string &genetic_map_filename,
                   const std::string &map_format,
                   const std::string &output_filename,
                   const std::string &output_format, bool output_morgans,
                   const double &step_interval, unsigned fixed_output_width,
                   bool verbose) const;
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_INTERPOLATOR_H_

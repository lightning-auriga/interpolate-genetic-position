/*!
 \file cargs.h
 \brief command line argument handling
 \author Lightning Auriga
 \note requires boost::program_options library + headers
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Thanks to
 https://www.boost.org/doc/libs/1_70_0/doc/html/program_options/tutorial.html#id-1.3.32.4.3
 */

#ifndef INTERPOLATE_GENETIC_POSITION_CARGS_H_
#define INTERPOLATE_GENETIC_POSITION_CARGS_H_

#include <stdexcept>
#include <string>
#include <vector>

#include "boost/program_options.hpp"

namespace interpolate_genetic_position {
/*!
  \class cargs
  \brief command line argument parser using boost::program_options
 */
class cargs {
 public:
  /*!
    \brief constructor with program arguments
    @param argc number of arguments including program name
    @param argv string array containing actual arguments
   */
  cargs(int argc, const char **const argv);
  /*!
    \brief copy constructor
    @param obj existing cargs object
   */
  cargs(const cargs &obj);
  /*!
    \brief destructor
   */
  ~cargs() throw();

  /*!
    \brief set user help documentation and default values for parameters as
    needed

    Note the weird syntax with overloaded () operators. The lists are not
    separated by commas.
   */
  void initialize_options();

  /*!
    \brief determine whether the user has requested help documentation
    \return whether the user has requested help documentation

    This test is separate from whether the user has run the software with no
    flags
   */
  bool help() const;

  /*!
   * \brief determine whether the user has requested (extremely) verbose debug
   * logging \return whether the user has requested verbose logging
   */
  bool verbose() const;

  /*!
   * \brief get name of input variant/region query file
   * \return name of input variant/query region file
   *
   * this string can be empty, in which case input is pulled from cin
   */
  std::string get_input_filename() const;
  /*!
   * \brief get input variant/region file format
   * \return input variant/region file format
   *
   * the term "preset" is based on tabix's notation.
   * accepted values are: "bim", "map", "bed". input
   * format also controls output behavior.
   */
  std::string get_input_preset() const;
  /*!
   * \brief get name of input recombination map
   * \return name of input recombination map
   */
  std::string get_recombination_map() const;
  /*!
   * \brief get input recombination map format
   * \return input recombination map format
   *
   * this refers to the source of the genetic map file.
   * accepted values are:
   *   - "bolt"
   * (https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/)
   *   - "bedgraph"
   * (https://hgdownload.soe.ucsc.edu/gbdb/hg38/recombRate/recombAvg.bw)
   *     - note that this needs to be converted to bedgraph before being used
   * with this software
   */
  std::string get_map_format() const;
  /*!
   * \brief get name of output file
   * \return name of output file
   */
  std::string get_output_filename() const;

  /*!
    \brief determine whether genetic position should be output in morgans,
    instead of centimorgans.
    \return whether the user has requested that genetic position be output
    in morgans instead of centimorgans.
   */
  bool output_morgans() const;

  /*!
   * \brief get fixed genetic distance to add to the boundary between
   * successive interpolated regions, when using bedfile input.
   * \return fixed genetic distance to add.
   *
   * This parameter is designed to enable experimental injection of artificial
   * genetic breakpoints at fixed regions of the genome. Interpolation within
   * the regions remains the same, but the fixed distance is added between
   * regions, creating a step function-like distance. There are only extremely
   * experimental applications of this functionality; the default of 0 will
   * cause standard interpolation to be performed without introduciing
   * artificial breakpoints, and should be used in almost all cases.
   */
  double get_region_step_interval() const;
  /*!
    \brief find status of arbitrary flag
    @param tag name of flag
    \return whether the flag is set

    This is the underlying accessor function used by the custom flag
    accessors, and can be used for arbitrary flag additions if you don't want
    to type out the custom access functions or want to allow dynamic
    specification from a config file.
   */
  bool compute_flag(const std::string &tag) const;
  /*!
    \brief find value of arbitrary parameter
    @tparam value_type class to which the value should be cast
    @param tag name of parameter
    \return value of parameter if specified

    \warning throws exception if parameter was not set and had no default
   */
  template <class value_type>
  value_type compute_parameter(const std::string &tag) const {
    if (_vm.count(tag)) {
      return _vm[tag].as<value_type>();
    }
    throw std::domain_error("cargs: requested parameter \"" + tag + "\" unset");
  }

  /*!
    \brief report help documentation to arbitrary output stream
    @param out stream to which to write help documentation

    Parameter should probably be std::cout or std::cerr at your preference.
   */
  void print_help(std::ostream &out);

 private:
  /*!
    \brief default constructor
    \warning disabled
   */
  cargs();
  boost::program_options::options_description
      _desc;  //!< help documentation string
  boost::program_options::variables_map
      _vm;  //!< storage of parsed command line settings
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_CARGS_H_

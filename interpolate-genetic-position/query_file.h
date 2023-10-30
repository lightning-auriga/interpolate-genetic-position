/*!
 \file query_file.h
 \brief interface for variant or region query file
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_QUERY_FILE_H_
#define INTERPOLATE_GENETIC_POSITION_QUERY_FILE_H_

#include <gmpxx.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "interpolate-genetic-position/input_variant_file.h"
#include "interpolate-genetic-position/utilities.h"

namespace interpolate_genetic_position {
/*!
 * \class query_file
 * \brief interface class for reading data from a variant or region
 * query file. handles different input and output formats, and adjusts
 * query formats as needed. in the initial implementation, this
 * file must be sorted by chromosome and position.
 */
class query_file {
 public:
  /*!
   * \brief default constructor
   */
  query_file();
  /*!
   * \brief copy constructor
   *
   * Copy constructor is disabled due to internal
   * file handle.
   */
  query_file(const query_file &obj);
  /*!
   * \brief initialize object with pointer to existing interface object
   * \param ptr pointer to allocated interface object. this is managed upstream,
   * and will not be deleted by this object
   */
  explicit query_file(base_input_variant_file *ptr);
  /*!
   * \brief destructor
   */
  ~query_file() throw();
  /*!
   * \brief initialize file connection from string
   * \param filename name of input file. if an empty string,
   * the input is assumed to come from cin, and logic skips
   * opening the stream but otherwise proceeds as usual
   * \param ft descriptor of query file format
   */
  void open(const std::string &filename, format_type ft);
  /*!
   * \brief initialize file connection from char
   * \param filename name of input file
   * \param ft descriptor of query file format
   * \note included for thematic compatibility with std::ifstream
   */
  void open(const char *filename, format_type ft);
  /*!
   * \brief get the next variant entry from file and
   * store it in internal buffer
   * \return boolean indicating whether a new file entry
   * was successfully loaded. FALSE should indicate EOF.
   */
  bool get();
  /*!
   * \brief get chromosome of current query
   * \return chromosome of current query
   */
  const std::string &get_chr() const;
  /*!
   * \brief get position of current query
   * \return position of current query
   */
  const mpz_class &get_pos() const;
  /*!
   * \brief close any input connection
   */
  void close();
  /*!
   * \brief probe input connection to determine whether EOF has been reached
   * \return whether an open connection is EOF
   *
   * If no connection is established, then input is coming from cin, in which
   * case that stream's next character is tested for EOF.
   */
  bool eof();
  /*!
   * \brief given an interpolated genetic position and an output stream,
   * recapitulate input data with updated genetic position information
   * \param gpos_interpolated interpolated genetic position for this query
   * \param output pointer to output stream; opening/closing it is handled
   * upstream
   */
  void report(const mpf_class &gpos_interpolated, std::ostream *output) const;

 private:
  base_input_variant_file *_interface;  //!< input file handler
  format_type _ft;                      //!< stored format of input filestream
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_QUERY_FILE_H_

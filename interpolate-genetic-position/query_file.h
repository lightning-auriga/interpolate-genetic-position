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
#include "interpolate-genetic-position/output_variant_file.h"
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
   * \param inptr pointer to allocated input interface object. this is
   * managed upstream, and will not be deleted by this object
   * \param outptr pointer to allocated output interface object. this is
   * managed upstream, and will not be deleted by this object
   */
  explicit query_file(base_input_variant_file *inptr,
                      base_output_variant_file *outptr);
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
   * \brief initialize output file connection
   * \param filename name of output file. if an empty string,
   * the output is assumed to go to cout, and logic skips
   * opening the stream but otherwise proceeds as usual
   * \param ft descriptor of output file format
   */
  void initialize_output(const std::string &filename, format_type ft);
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
  const mpz_class &get_pos1() const;
  /*!
   * \brief get end position of current query
   * \return end position of current query or, if not
   * applicable for this query type, -1
   */
  const mpz_class &get_pos2() const;
  /*!
   * \brief get the fixed interval to add to output
   * \return fixed interval to add to output
   */
  const mpf_class &get_step_interval() const;
  /*!
   * \brief set the fixed interval to add to output
   * \param step fixed interval to add to output
   */
  void set_step_interval(const mpf_class &step);
  /*!
   * \brief get the index of this result among results
   * for this chromosome
   * \return index of this result
   */
  unsigned get_index_on_chromosome() const;
  /*!
   * \brief set the index of this result among results
   * for this chromosome
   * \param index index of this result
   */
  void set_index_on_chromosome(unsigned index);
  /*!
   * \brief get identifier of previous output result
   * \return identifier of previous output result
   */
  const std::string &get_previous_chromosome() const;
  /*!
   * \brief set identifier of previous output result
   * \param chr new value for chromosome
   */
  void set_previous_chromosome(const std::string &chr);
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
   * \brief given an interpolated result and an output stream,
   * report the result to an output stream
   * \param results a vector of results corresponding to an input query
   */
  void report(const std::vector<query_result> &results);

 private:
  base_input_variant_file *_interface;  //!< input file handler
  format_type _ft;                      //!< stored format of input filestream
  base_output_variant_file *_output;    //!< interface class to output
  unsigned _index_on_chromosome;     //!< how many queries have been returned on
                                     //!< this chromosome
  std::string _previous_chromosome;  //!< name of chromosome for prior output
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_QUERY_FILE_H_

/*!
 \file genetic_map.h
 \brief interface for reference genetic map file
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_GENETIC_MAP_H_
#define INTERPOLATE_GENETIC_POSITION_GENETIC_MAP_H_

#include <gmpxx.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "interpolate-genetic-position/utilities.h"

namespace interpolate_genetic_position {
/*!
 * \class genetic_map
 * \brief interface class for reading data from a genetic map file.
 * attempts to buffer data to prevent excessive memory allocation,
 * though in the initial implementation at least, this forces input
 * query variants to be appropriately sorted.
 */
class genetic_map {
 public:
  /*!
   * \brief default constructor
   */
  genetic_map();
  /*!
   * \brief copy constructor
   *
   * Copy constructor is disabled due to internal
   * file handle.
   */
  genetic_map(const genetic_map &obj);
  /*!
   * \brief destructor
   */
  ~genetic_map() throw();
  /*!
   * \brief initialize file connection from string
   * \param filename name of input file
   * \param ft descriptor of recombination rate file format
   */
  void open(const std::string &filename, format_type ft);
  /*!
   * \brief initialize file connection from char
   * \param filename name of input file
   * \param ft descriptor of recombination rate file format
   * \note included for thematic compatibility with std::ifstream
   */
  void open(const char *filename, format_type ft);
  /*!
   * \brief get the next genetic map entry from file and
   * store it in internal buffer
   * \return boolean indicating whether a new file entry
   * was successfully loaded. FALSE should indicate EOF.
   */
  bool get();
  /*!
   * \brief query the currently loaded data points to try
   * to interpolate genetic position. this action may or may
   * not cause the object to update its internal cache.
   * \param chr_query chromosome of query variant as string
   * \param pos_query physical position of query variant,
   * represented as mpz
   * \param verbose whether to emit (extremely) verbose logging to std::cout
   * \param gpos_interpolated pointer to mpf with interpolated
   * genetic position for the query
   */
  void query(const std::string &chr_query, const mpz_class &pos_query,
             bool verbose, mpf_class *gpos_interpolated);
  /*!
   * \brief close any input connection
   */
  void close();
  /*!
   * \brief probe input connection to determine whether EOF has been reached
   * \return whether an open connection is EOF
   *
   * If no connection is established... not really sure.
   */
  bool eof();
  /*!
   * \brief translate an assortment of chromosome representations
   * into simple integers for sort comparison.
   * \param chr chromosome code, as string, for conversion
   * \return input chromosome, as an integer on {[1,24], 26}
   *
   * This function should handle optional "chr" prefixes depending on input
   * build. Chromosomes 1-22 are labeled as normal; X is interpreted as 23; Y is
   * interpreted as 24; M or MT is interpreted as 26. XY/25, the PAR regions,
   * are not currently handled.
   */
  int chromosome_to_integer(const std::string &chr) const;
  /*!
   * \brief compare two chromosome representations to get a consistent
   * directionality of their comparison. Will strip prefixes and convert
   * most standard chromosome aliases.
   * \param chr1 first chromosome code to compare
   * \param chr2 second chromosome code to compare
   * \return the relationship of the first chromosome code to the second
   */
  direction chromosome_compare(const std::string &chr1,
                               const std::string &chr2) const;

 private:
  std::ifstream _input;   //!< file connection for plaintext input
  gzFile _gzinput;        //!< file connection for gzipped input
  char *_buffer;          //!< character buffer for gzinput line reading
  unsigned _buffer_size;  //!< size of gzinput buffer, in bytes
  format_type _ft;        //!< stored format of input filestream
  std::string _chr_old;   //!< chromosome of previous entry
  std::string _chr_new;   //!< chromosome of new entry
  mpz_class _pos_old;     //!< physical position of previous entry
  mpz_class _pos_new;     //!< physical position of new entry
  mpf_class _gpos_old;    //!< genetic position of previous entry
  mpf_class _gpos_new;    //!< genetic position of new entry
  mpf_class _rate_old;    //!< point recombination rate change of previous entry
  mpf_class _rate_new;    //!< point recombination rate change of new entry
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_GENETIC_MAP_H_

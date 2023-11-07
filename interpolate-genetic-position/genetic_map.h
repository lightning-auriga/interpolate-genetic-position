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

#include "interpolate-genetic-position/input_genetic_map_file.h"
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
   * \brief construct object with pointer to object
   * interfacing with genetic map file
   * \param ptr pointer to interface object. memory is
   * handled upstream and won't be freed by this class
   */
  explicit genetic_map(base_input_genetic_map_file *ptr);
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
   * \brief run possibly multiple queries over segments
   * of a larger query, aggregate the results, and return
   * them at once
   * \param chr_query chromosome of query variant as string
   * \param pos1_query physical position of query variant,
   * represented as mpz
   * \param pos2_query for applicable queries (e.g. bedfiles),
   * physical position of end of region; else, -1
   * \param verbose whether to emit (extremely) verbose logging to std::cout
   * \param results pointer to vector that should contain interpolated results
   */
  void query(const std::string &chr_query, const mpz_class &pos1_query,
             const mpz_class &pos2_query, bool verbose,
             std::vector<query_result> *results);
  /*!
   * \brief query the currently loaded data points to try
   * to interpolate genetic position. this action may or may
   * not cause the object to update its internal cache.
   * \param chr_query chromosome of query variant as string
   * \param pos1_query physical position of query variant,
   * represented as mpz
   * \param pos2_query for applicable queries (e.g. bedfiles),
   * physical position of end of region; else, -1
   * \param verbose whether to emit (extremely) verbose logging to std::cout
   * \param result pointer to object that should contain interpolated results
   */
  void query(const std::string &chr_query, const mpz_class &pos1_query,
             const mpz_class &pos2_query, bool verbose, query_result *results);
  /*!
   * \brief close any input connection
   */
  void close();

 private:
  base_input_genetic_map_file *_interface;  //!< pointer to interface object
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_GENETIC_MAP_H_

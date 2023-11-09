/*!
 \file output_variant_file.h
 \brief interface for writing output variant files
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_OUTPUT_VARIANT_FILE_H_
#define INTERPOLATE_GENETIC_POSITION_OUTPUT_VARIANT_FILE_H_

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
 * \class base_output_variant_file
 * \brief virtual base class for output variant file types;
 * split out in this manner to facilitate mocks.
 */
class base_output_variant_file {
 public:
  /*!
   * \brief basic constructor
   */
  base_output_variant_file();
  /*!
   * \brief pure virtual destructor
   */
  virtual ~base_output_variant_file() throw() = 0;
  /*!
   * \brief open file connection
   * \param filename name of input file to open
   * \param ft format of output file
   */
  virtual void open(const std::string &filename, format_type ft) = 0;
  /*!
   * \brief close file connection
   */
  virtual void close() = 0;
  /*!
   * \brief report output data to appropriate target
   */
  virtual void write(const std::string &chr, const mpz_class &pos1,
                     const mpz_class &pos2, const std::string &id,
                     const mpf_class &gpos, const std::string &a1,
                     const std::string &a2) = 0;
  /*!
   * \brief get descriptor of format of output file
   * \return output file format
   */
  virtual format_type get_format() const = 0;
  /*!
   * \brief set whether output unit should be morgans
   * \param use_morgans whether output unit should be morgans
   */
  virtual void output_morgans(bool use_morgans) = 0;
  /*!
   * \brief determine whether output genetic position unit
   * should be morgans
   * \return whether output genetic position unit
   * should be morgans
   */
  virtual bool output_morgans() const = 0;
};

/*!
 * \class output_variant_file
 * \brief the base class for these file interfaces needs to be purely virtual
 * for unit testing, but many of the derived classes' methods share
 * identical implementations.
 */
class output_variant_file : public base_output_variant_file {
 public:
  /*!
   * \brief basic constructor
   */
  output_variant_file();
  /*!
   * \brief destructor
   */
  ~output_variant_file() throw();
  /*!
   * \brief open file connection
   * \param filename name of input file to open
   * \param ft format of output file
   */
  void open(const std::string &filename, format_type ft);
  /*!
   * \brief close file connection
   */
  void close();
  /*!
   * \brief report output data to appropriate target
   */
  void write(const std::string &chr, const mpz_class &pos1,
             const mpz_class &pos2, const std::string &id,
             const mpf_class &gpos, const std::string &a1,
             const std::string &a2);
  /*!
   * \brief get descriptor of format of output file
   * \return output file format
   */
  format_type get_format() const;
  /*!
   * \brief set whether output unit should be morgans
   * \param use_morgans whether output unit should be morgans
   */
  void output_morgans(bool use_morgans);
  /*!
   * \brief determine whether output genetic position unit
   * should be morgans
   * \return whether output genetic position unit
   * should be morgans
   */
  bool output_morgans() const;

 private:
  std::ofstream _output;  //!< output uncompressed file stream
  format_type _ft;        //!< format of output file
  bool _output_morgans;   //!< whether to emit genetic position as morgans
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_OUTPUT_VARIANT_FILE_H_

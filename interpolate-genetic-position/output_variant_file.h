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
                     const mpf_class &gpos, const mpf_class &rate,
                     const std::string &a1, const std::string &a2) = 0;
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
  /*!
   * \brief get the chromosome code of the most recently
   * emitted result
   * \return the chromosome code of the most recently
   * emitted result
   */
  virtual std::string get_last_chr() const = 0;
  /*!
   * \brief set the chromosome code of the most recently
   * emitted result
   * \param chr the new chromosome code of the most recently
   * emitted result
   */
  virtual void set_last_chr(const std::string &chr) = 0;
  /*!
   * \brief get the start position of the most recently
   * emitted result
   * \return the start position of the most recently
   * emitted result
   */
  virtual mpz_class get_last_pos1() const = 0;
  /*!
   * \brief set the start position of the most recently
   * emitted result
   * \param pos1 the new start position of the most recently
   * emitted result
   */
  virtual void set_last_pos1(const mpz_class &pos1) = 0;
  /*!
   * \brief get the end position of the most recently
   * emitted result
   * \return the end position of the most recently
   * emitted result
   */
  virtual mpz_class get_last_pos2() const = 0;
  /*!
   * \brief set the end position of the most recently
   * emitted result
   * \param pos1 the new end position of the most recently
   * emitted result
   */
  virtual void set_last_pos2(const mpz_class &pos2) = 0;
  /*!
   * \brief get the genetic position of the most recently
   * emitted result
   * \return the genetic position of the most recently
   * emitted result
   */
  virtual mpf_class get_last_gpos() const = 0;
  /*!
   * \brief set the genetic position of the most recently
   * emitted result
   * \param gpos the new genetic position of the most recently
   * emitted result
   */
  virtual void set_last_gpos(const mpf_class &gpos) = 0;
  /*!
   * \brief get the rate of the most recently
   * emitted result
   * \return the rate of the most recently
   * emitted result
   */
  virtual mpf_class get_last_rate() const = 0;
  /*!
   * \brief set the rate of the most recently
   * emitted result
   * \param gpos the new rate of the most recently
   * emitted result
   */
  virtual void set_last_rate(const mpf_class &rate) = 0;
  /*!
   * \brief get the experimental step interval at region boundaries
   * \return the experimental step interval at region boundaries
   */
  virtual const mpf_class &get_step_interval() const = 0;
  /*!
   * \brief set the experimental step interval at region boundaries
   * \param step_interval the experimental step interval
   */
  virtual void set_step_interval(const mpf_class &step_interval) = 0;
  /*!
   * \brief get the index of this result among results
   * for this chromosome
   * \return index of this result
   */
  virtual unsigned get_index_on_chromosome() const = 0;
  /*!
   * \brief set the index of this result among results
   * for this chromosome
   * \param index index of this result
   */
  virtual void set_index_on_chromosome(unsigned index) = 0;
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
             const mpf_class &gpos, const mpf_class &rate,
             const std::string &a1, const std::string &a2);
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
  /*!
   * \brief get the chromosome code of the most recently
   * emitted result
   * \return the chromosome code of the most recently
   * emitted result
   */
  std::string get_last_chr() const;
  /*!
   * \brief set the chromosome code of the most recently
   * emitted result
   * \param chr the new chromosome code of the most recently
   * emitted result
   */
  void set_last_chr(const std::string &chr);
  /*!
   * \brief get the start position of the most recently
   * emitted result
   * \return the start position of the most recently
   * emitted result
   */
  mpz_class get_last_pos1() const;
  /*!
   * \brief set the start position of the most recently
   * emitted result
   * \param pos1 the new start position of the most recently
   * emitted result
   */
  void set_last_pos1(const mpz_class &pos1);
  /*!
   * \brief get the end position of the most recently
   * emitted result
   * \return the end position of the most recently
   * emitted result
   */
  mpz_class get_last_pos2() const;
  /*!
   * \brief set the end position of the most recently
   * emitted result
   * \param pos2 the new end position of the most recently
   * emitted result
   */
  void set_last_pos2(const mpz_class &pos2);
  /*!
   * \brief get the genetic position of the most recently
   * emitted result
   * \return the genetic position of the most recently
   * emitted result
   */
  mpf_class get_last_gpos() const;
  /*!
   * \brief set the genetic position of the most recently
   * emitted result
   * \param gpos the new genetic position of the most recently
   * emitted result
   */
  void set_last_gpos(const mpf_class &gpos);
  /*!
   * \brief get the rate of the most recently
   * emitted result
   * \return the rate of the most recently
   * emitted result
   */
  mpf_class get_last_rate() const;
  /*!
   * \brief set the rate of the most recently
   * emitted result
   * \param gpos the new rate of the most recently
   * emitted result
   */
  void set_last_rate(const mpf_class &rate);
  /*!
   * \brief get the experimental step interval at region boundaries
   * \return the experimental step interval at region boundaries
   */
  const mpf_class &get_step_interval() const;
  /*!
   * \brief set the experimental step interval at region boundaries
   * \param step_interval the experimental step interval
   */
  void set_step_interval(const mpf_class &step_interval);
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

 private:
  std::ofstream _output;     //!< output uncompressed file stream
  format_type _ft;           //!< format of output file
  bool _output_morgans;      //!< whether to emit genetic position as morgans
  std::string _last_chr;     //!< chromosome of most recently emitted result
  mpz_class _last_pos1;      //!< start position of most recent result
  mpz_class _last_pos2;      //!< end position of most recent result
  mpf_class _last_gpos;      //!< genetic position of most recent result
  mpf_class _last_rate;      //!< rate of most recent result
  mpf_class _step_interval;  //!< gpos increment to edge of bed queries
  unsigned _index_on_chromosome;  //!< how many queries have been returned on
                                  //!< this chromosome
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_OUTPUT_VARIANT_FILE_H_

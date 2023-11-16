/*!
 \file input_genetic_map.h
 \brief interface for reading input genetic maps
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_INPUT_GENETIC_MAP_FILE_H_
#define INTERPOLATE_GENETIC_POSITION_INPUT_GENETIC_MAP_FILE_H_

#include <gmpxx.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "interpolate-genetic-position/bigwig_reader.h"
#include "interpolate-genetic-position/utilities.h"

namespace interpolate_genetic_position {
/*!
 * \class base_input_genetic_map_file
 * \brief virtual base class for input genetic map files;
 * split out in this manner to facilitate mocks.
 */
class base_input_genetic_map_file {
 public:
  /*!
   * \brief basic constructor
   */
  base_input_genetic_map_file();
  /*!
   * \brief pure virtual destructor
   */
  virtual ~base_input_genetic_map_file() throw() = 0;
  /*!
   * \brief open file connection
   * \param filename name of input file to open
   * \param ft format of input genetic map
   */
  virtual void open(const std::string &filename, format_type ft) = 0;
  /*!
   * \brief set the fallback stream for data
   * \param ptr pointer to the fallback stream
   *
   * This is intended to be cin, but is exposed for testing
   */
  virtual void set_fallback_stream(std::istream *ptr) = 0;
  /*!
   * \brief get the fallback stream for data
   * \return pointer to the fallback stream
   *
   * This is intended to be cin, but is exposed for testing
   */
  virtual std::istream *get_fallback_stream() const = 0;
  /*!
   * \brief get the next genetic map entry from file and
   * store it in internal buffer
   * \return boolean indicating whether a new file entry
   * was successfully loaded. FALSE should indicate EOF.
   */
  virtual bool get() = 0;
  /*!
   * \brief close file connection
   */
  virtual void close() = 0;
  /*!
   * \brief test input connection for EOF
   * \return whether input connection has encountered (or will encounter
   * as next character) EOF
   */
  virtual bool eof() = 0;
  /*!
   * \brief get chromosome of lower boundary of cached range
   * \return chromosome of lower boundary of cached range
   */
  virtual std::string get_chr_lower_bound() const = 0;
  /*!
   * \brief get chromosome of upper boundary of cached range
   * \return chromosome of upper boundary of cached range
   */
  virtual std::string get_chr_upper_bound() const = 0;
  /*!
   * \brief get start physical position of lower boundary of cached range
   * \return start physical position of lower boundary of cached range
   */
  virtual mpz_class get_startpos_lower_bound() const = 0;
  /*!
   * \brief get start physical position of upper boundary of cached range
   * \return start physical position of upper boundary of cached range
   */
  virtual mpz_class get_startpos_upper_bound() const = 0;
  /*!
   * \brief get end physical position of lower boundary of cached range
   * \return end physical position of lower boundary of cached range
   */
  virtual mpz_class get_endpos_lower_bound() const = 0;
  /*!
   * \brief get end physical position of upper boundary of cached range
   * \return end physical position of upper boundary of cached range
   */
  virtual mpz_class get_endpos_upper_bound() const = 0;
  /*!
   * \brief get genetic position of lower boundary of cached range
   * \return genetic position of lower boundary of cached range
   */
  virtual mpf_class get_gpos_lower_bound() const = 0;
  /*!
   * \brief get genetic position of upper boundary of cached range
   * \return genetic position of upper boundary of cached range
   */
  virtual mpf_class get_gpos_upper_bound() const = 0;
  /*!
   * \brief get rate of change of genetic distance at lower boundary of cached
   * range \return rate of change of genetic distance at lower boundary of
   * cached range
   */
  virtual mpf_class get_rate_lower_bound() const = 0;
  /*!
   * \brief get rate of change of genetic distance at upper boundary of cached
   * range \return rate of change of genetic distance at upper boundary of
   * cached range
   */
  virtual mpf_class get_rate_upper_bound() const = 0;
};

/*!
 * \class input_genetic_map_file
 */
class input_genetic_map_file : public base_input_genetic_map_file {
 public:
  /*!
   * \brief basic constructor
   */
  input_genetic_map_file();
  /*!
   * \brief destructor
   */
  ~input_genetic_map_file() throw();
  /*!
   * \brief open file connection
   * \param filename name of input file to open
   * \param ft format of input genetic map
   */
  void open(const std::string &filename, format_type ft);
  /*!
   * \brief set the fallback stream for data
   * \param ptr pointer to the fallback stream
   *
   * This is intended to be cin, but is exposed for testing
   */
  void set_fallback_stream(std::istream *ptr);
  /*!
   * \brief get the fallback stream for data
   * \return pointer to the fallback stream
   *
   * This is intended to be cin, but is exposed for testing
   */
  std::istream *get_fallback_stream() const;
  /*!
   * \brief get the next genetic map entry from file and
   * store it in internal buffer
   * \return boolean indicating whether a new file entry
   * was successfully loaded. FALSE should indicate EOF.
   */
  bool get();
  /*!
   * \brief close file connection
   */
  void close();
  /*!
   * \brief test input connection for EOF
   * \return whether input connection has encountered (or will encounter
   * as next character) EOF
   */
  bool eof();
  /*!
   * \brief get chromosome of lower boundary of cached range
   * \return chromosome of lower boundary of cached range
   */
  std::string get_chr_lower_bound() const;
  /*!
   * \brief get chromosome of upper boundary of cached range
   * \return chromosome of upper boundary of cached range
   */
  std::string get_chr_upper_bound() const;
  /*!
   * \brief get start physical position of lower boundary of cached range
   * \return start physical position of lower boundary of cached range
   */
  mpz_class get_startpos_lower_bound() const;
  /*!
   * \brief get start physical position of upper boundary of cached range
   * \return start physical position of upper boundary of cached range
   */
  mpz_class get_startpos_upper_bound() const;
  /*!
   * \brief get end physical position of lower boundary of cached range
   * \return end physical position of lower boundary of cached range
   */
  mpz_class get_endpos_lower_bound() const;
  /*!
   * \brief get end physical position of upper boundary of cached range
   * \return end physical position of upper boundary of cached range
   */
  mpz_class get_endpos_upper_bound() const;
  /*!
   * \brief get genetic position of lower boundary of cached range
   * \return genetic position of lower boundary of cached range
   */
  mpf_class get_gpos_lower_bound() const;
  /*!
   * \brief get genetic position of upper boundary of cached range
   * \return genetic position of upper boundary of cached range
   */
  mpf_class get_gpos_upper_bound() const;
  /*!
   * \brief get rate of change of genetic distance at lower boundary of cached
   * range \return rate of change of genetic distance at lower boundary of
   * cached range
   */
  mpf_class get_rate_lower_bound() const;
  /*!
   * \brief get rate of change of genetic distance at upper boundary of cached
   * range \return rate of change of genetic distance at upper boundary of
   * cached range
   */
  mpf_class get_rate_upper_bound() const;

 private:
  std::ifstream _input;          //!< file connection for plaintext input
  gzFile _gzinput;               //!< file connection for gzipped input
  std::istream *_fallback;       //!< pointer to fallback stream connection
  bigwig_reader _bwinput;        //!< file connection for bigwigs
  char *_buffer;                 //!< character buffer for gzinput line reading
  unsigned _buffer_size;         //!< size of gzinput buffer, in bytes
  format_type _ft;               //!< stored format of input filestream
  std::string _chr_lower_bound;  //!< chromosome of previous entry
  std::string _chr_upper_bound;  //!< chromosome of new entry
  mpz_class
      _startpos_lower_bound;  //!< start physical position of previous entry
  mpz_class _startpos_upper_bound;  //!< start physical position of new entry
  mpz_class _endpos_lower_bound;    //!< end physical position of previous entry
  mpz_class _endpos_upper_bound;    //!< end physical position of new entry
  mpf_class _gpos_lower_bound;      //!< genetic position of previous entry
  mpf_class _gpos_upper_bound;      //!< genetic position of new entry
  mpf_class
      _rate_lower_bound;  //!< point recombination rate change of previous entry
  mpf_class
      _rate_upper_bound;  //!< point recombination rate change of new entry
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_INPUT_GENETIC_MAP_FILE_H_

/*!
 \file input_file.h
 \brief interface for reading input variant files
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_INPUT_VARIANT_FILE_H_
#define INTERPOLATE_GENETIC_POSITION_INPUT_VARIANT_FILE_H_

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
 * \class base_input_variant_file
 * \brief virtual base class for input variant file types;
 * split out in this manner to facilitate mocks.
 */
class base_input_variant_file {
 public:
  /*!
   * \brief basic constructor
   */
  base_input_variant_file();
  /*!
   * \brief pure virtual destructor
   */
  virtual ~base_input_variant_file() throw() = 0;
  /*!
   * \brief open file connection
   * \param filename name of input file to open
   */
  virtual void open(const std::string &filename) = 0;
  /*!
   * \brief close file connection
   */
  virtual void close() = 0;
  /*!
   * \brief set parameters controlling interpretation
   * of input lines
   * \param chr_index base 0 index of chromosome in line tokens
   * \param pos_index base 0 index of physical position in line tokens
   * \param gpos_index base 0 index of genetic position in line tokens
   * \param base0 whether physical position is base 0
   * \param n_tokens expected number of tokens in line
   */
  virtual void set_format_parameters(unsigned chr_index, unsigned pos_index,
                                     unsigned gpos_index, bool base0,
                                     unsigned n_tokens) = 0;
  /*!
   * \brief get a single line from an input stream
   * \param line pointer to allocated string that will contain line,
   * or equivalent, depending on input type
   * \return whether access operation was successful
   */
  virtual bool get_input_line(std::string *line) = 0;
  /*!
   * \brief get next marker's metadata from input connection
   * \return whether a variant was successfully accessed from file
   */
  virtual bool get_variant() = 0;
  /*!
   * \brief get chromosome of currently loaded marker
   * \return chromosome of currently loaded marker
   */
  virtual const std::string &get_chr() const = 0;
  /*!
   * \brief get physical position of currently loaded marker
   * \return physical position of currently loaded marker
   */
  virtual const mpz_class &get_pos() const = 0;
  /*!
   * \brief get vector containing tokenized representation of current line
   * \return vector containing tokenized representation of current line
   */
  virtual const std::vector<std::string> &get_line_contents() const = 0;
  /*!
   * \brief test input connection for EOF
   * \return whether input connection has encountered (or will encounter
   * as next character) EOF
   */
  virtual bool eof() = 0;
};

/*!
 * \class input_variant_file
 * \brief the base class for these file interfaces needs to be purely virtual
 * for unit testing, but many of the derived classes' methods share
 * identical implementations.
 */
class input_variant_file : public base_input_variant_file {
 public:
  /*!
   * \brief basic constructor
   */
  input_variant_file();
  /*!
   * \brief destructor
   */
  ~input_variant_file() throw();
  /*!
   * \brief open file connection
   * \param filename name of input file to open
   */
  void open(const std::string &filename);
  /*!
   * \brief close file connection
   */
  void close();
  /*!
   * \brief set parameters controlling interpretation
   * of input lines
   * \param chr_index base 0 index of chromosome in line tokens
   * \param pos_index base 0 index of physical position in line tokens
   * \param gpos_index base 0 index of genetic position in line tokens
   * \param base0 whether physical position is base 0
   * \param n_tokens expected number of tokens in line
   */
  void set_format_parameters(unsigned chr_index, unsigned pos_index,
                             unsigned gpos_index, bool base0,
                             unsigned n_tokens);
  /*!
   * \brief get a single line from an input stream
   * \param line pointer to allocated string that will contain line,
   * or equivalent, depending on input type
   * \return whether access operation was successful
   */
  bool get_input_line(std::string *line);
  /*!
   * \brief get next marker's metadata from input connection
   * \return whether a variant was successfully accessed from file
   */
  bool get_variant();
  /*!
   * \brief get chromosome of currently loaded marker
   * \return chromosome of currently loaded marker
   */
  const std::string &get_chr() const;
  /*!
   * \brief get physical position of currently loaded marker
   * \return physical position of currently loaded marker
   */
  const mpz_class &get_pos() const;
  /*!
   * \brief get vector containing tokenized representation of current line
   * \return vector containing tokenized representation of current line
   */
  const std::vector<std::string> &get_line_contents() const;
  /*!
   * \brief test input connection for EOF
   * \return whether input connection has encountered (or will encounter
   * as next character) EOF
   */
  bool eof();

 private:
  std::ifstream _input;                     //!< input uncompressed file stream
  gzFile _gzinput;                          //!< input gzipped file pointer
  char *_buffer;                            //!< character buffer for zlib reads
  std::vector<std::string> _line_contents;  //!< tokenized input line
  unsigned _buffer_size;                    //!< size of allocated buffer
  unsigned _chr_index;                      //!< index of chromosome in line
  unsigned _pos_index;   //!< index of physical position in line
  unsigned _gpos_index;  //!< index of genetic position in line
  std::string _chr;      //!< chromosome of current marker
  mpz_class _pos;        //!< physical position of current marker
  bool _base0;           //!< whether physical position is base 0
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_INPUT_VARIANT_FILE_H_

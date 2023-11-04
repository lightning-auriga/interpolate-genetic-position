/*!
 \file bigwig_reader.h
 \brief interface library for reading directly from bigwigs.
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_BIGWIG_READER_H_
#define INTERPOLATE_GENETIC_POSITION_BIGWIG_READER_H_

#include <bigWig.h>
#include <gmpxx.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "interpolate-genetic-position/utilities.h"

namespace interpolate_genetic_position {
/*!
 * \class bigwig_reader
 * \brief C++ style interface on top of libBigWig for streaming
 * reading of bigwig genetic recombination data.
 */
class bigwig_reader {
 public:
  /*!
   * \brief default contructor
   *
   * This does *not* allocate internal memory.
   */
  bigwig_reader();
  /*!
   * \brief default destructor
   *
   * This calls close.
   */
  ~bigwig_reader() throw();
  /*!
   * \brief initialize libBigWig internals and open file pointer
   * \param filename name of bigwig file to open
   *
   * Reportedly, this can be a remote.
   */
  void open(const std::string &filename);
  /*!
   * \brief free internal memory, close file, destruct library allocation
   */
  void close();
  /*!
   * \brief get next entry from bigwig file. how exactly the progression
   * between chromosomes will be handled is not yet clear.
   * \param chr pointer to string for result chromosome
   * \param pos1 pointer to arbitrary precision integer for start position
   * \param pos2 pointer to arbitrary precision integer for end position
   * \param rate pointer to arbitrary precision float for rate (cm/mb)
   * \return whether access operation was successful
   */
  bool get(std::string *chr, mpz_class *pos1, mpz_class *pos2, mpf_class *rate);
  /*!
   * \brief get chromosome of currently loaded overlapping intervals
   * \return chromosome of currently loaded overlapping intervals
   *
   * Note that if no intervals are currently loaded, this will return
   * an empty string.
   */
  const std::string &get_loaded_chr() const;
  /*!
   * \brief determine whether there is an active file connection
   * \return whether there is an active file connection
   */
  bool is_open() const;
  /*!
   * \brief determine whether entire bigwig has been iterated
   * \return whether valid chromosomes have been exhausted
   */
  bool eof() const;
  /*!
   * \brief load the next chromosome in the sorted sequence
   * \return whether the next chromosome could be loaded and the
   * current file is not EOF
   */
  bool load_next_chr();

 protected:
  /*!
   * \brief attempt to load all region data for a specified chromosome
   * \param chr requested chromosome
   * \return whether load operation returned anything
   */
  bool load_chr(const std::string &chr);

 private:
  bigWigFile_t *_input;                  //!< file handle to bigwig file
  bwOverlappingIntervals_t *_intervals;  //!< loaded bedgraph regions
  std::string _chr;          //!< chromosome range currently stored in intervals
  unsigned _interval_index;  //!< current accessor location in interval
};
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_BIGWIG_READER_H_

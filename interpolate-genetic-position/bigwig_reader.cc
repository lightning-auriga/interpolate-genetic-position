/*!
 \file bigwig_reader.cc
 \brief implementation for bigwig reader class
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/bigwig_reader.h"

namespace igp = interpolate_genetic_position;

igp::bigwig_reader::bigwig_reader()
    : _input(NULL), _intervals(NULL), _chr(""), _interval_index(0) {}

igp::bigwig_reader::~bigwig_reader() throw() {
  close();
  bwCleanup();
}

void igp::bigwig_reader::open(const std::string &filename) {
  if (!bwIsBigWig(filename.c_str(), NULL)) {
    throw std::runtime_error("bigwig_reader: input file is not bigwig: \"" +
                             filename + "\"");
  }
  if (bwInit(1 << 17)) {
    throw std::runtime_error("bigwig_reader: unable to bwInit");
  }
  _input = bwOpen(filename.c_str(), NULL, "r");
  if (!_input) {
    throw std::runtime_error("bigwig_reader: cannot open \"" + filename + "\"");
  }
  load_chr(find_minimum_chromosome());
}

std::string igp::bigwig_reader::find_minimum_chromosome() const {
  std::string bwchr = "";
  int bigwig_chrint = 0, min_bigwig_chrint = 777777;
  if (_input && _input->cl) {
    for (int64_t i = 0; i < _input->cl->nKeys; ++i) {
      bwchr = std::string(_input->cl->chrom[i]);
      if (chromosome_to_integer(bwchr, &bigwig_chrint)) {
        if (bigwig_chrint < min_bigwig_chrint) {
          min_bigwig_chrint = bigwig_chrint;
        }
      }
    }
    if (min_bigwig_chrint == 777777) {
      throw std::runtime_error("bigwig_reader: min chromosome not recognized");
    }
    return integer_to_chromosome(min_bigwig_chrint);
  }
  throw std::runtime_error(
      "bigwig_reader: find_minimum_chromosome was called without "
      "an open file handle, which is an undefined operation.");
}

void igp::bigwig_reader::close() {
  if (_input) {
    bwClose(_input);
    _input = NULL;
    bwCleanup();
  }
  if (_intervals) {
    bwDestroyOverlappingIntervals(_intervals);
    _intervals = NULL;
  }
  _chr = "";
  _interval_index = 0;
}

const std::string &igp::bigwig_reader::get_loaded_chr() const { return _chr; }

bool igp::bigwig_reader::load_chr(const std::string &chr) {
  _interval_index = 0;
  if (_intervals) {
    bwDestroyOverlappingIntervals(_intervals);
    _intervals = 0;
  }
  _chr = interpret_chr(chr);
  _intervals = bwGetOverlappingIntervals(_input, _chr.c_str(), 0, 1000000000);
  // failed load is denoted by NULL return pointer
  return _intervals != NULL;
}

bool igp::bigwig_reader::load_next_chr() {
  if (eof()) return false;
  return load_chr(next_chromosome(get_loaded_chr()));
}

bool igp::bigwig_reader::is_open() const { return _input != NULL; }

bool igp::bigwig_reader::eof() const {
  return _chr == "chrM" && (!_intervals || _interval_index == _intervals->l);
}

bool igp::bigwig_reader::get(std::string *chr, mpz_class *pos1, mpz_class *pos2,
                             mpf_class *rate) {
  if (!_intervals || eof()) return false;
  if (_interval_index == _intervals->l) {
    while (!eof()) {
      if (load_next_chr()) break;
    }
  }
  if (eof()) return false;
  *chr = get_loaded_chr();
  *pos1 = _intervals->start[_interval_index];
  *pos2 = _intervals->end[_interval_index];
  *rate = _intervals->value[_interval_index];
  ++_interval_index;
  return true;
}

std::string igp::bigwig_reader::interpret_chr(const std::string &chr) const {
  /*
   * From libBigWig's docs, it seems that a bigwig contains a chromosome
   * index in its leading metadata, and that index is available in the
   * loaded file object.
   *
   * The class anticipates that this will be called a very small number
   * of times relative to the total number of file operations. As such,
   * this is gonna be sloppy and not record its work for the moment.
   */
  std::string bwchr = "";
  int query_chrint = 0, bigwig_chrint = 0;
  if (!chromosome_to_integer(chr, &query_chrint)) {
    throw std::runtime_error("unrecognized query chromosome: \"" + chr + "\"");
  }
  if (_input && _input->cl) {
    for (int64_t i = 0; i < _input->cl->nKeys; ++i) {
      bwchr = std::string(_input->cl->chrom[i]);
      if (chromosome_to_integer(bwchr, &bigwig_chrint)) {
        if (query_chrint == bigwig_chrint) {
          return bwchr;
        }
      }
    }
    throw std::runtime_error("bigwig_reader: query chromosome \"" + chr +
                             "\" was not detected in any recognized form "
                             "in the genetic map.");
  }
  throw std::runtime_error(
      "bigwig_reader: interpret_chr was called without "
      "an open file handle, which is an undefined operation.");
}

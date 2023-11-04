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
  if (bwInit(1 << 17)) {
    throw std::runtime_error("bigwig_reader: unable to bwInit");
  }
  _input = bwOpen(filename.c_str(), NULL, "r");
  if (!_input) {
    throw std::runtime_error("bigwig_reader: cannot open \"" + filename + "\"");
  }
  load_chr("chr1");
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
  _chr = chr;
  _intervals = bwGetOverlappingIntervals(_input, chr.c_str(), 0, 1000000000);
  // failed load is denoted by NULL return pointer
  return _intervals != NULL;
}

bool igp::bigwig_reader::load_next_chr() {
  if (eof()) return false;
  return load_chr(next_chromosome(get_loaded_chr()));
}

bool igp::bigwig_reader::eof() const {
  return _chr == "chrM" && (!_intervals || _interval_index == _intervals->l);
}

bool igp::bigwig_reader::get(std::string *chr, mpz_class *pos1, mpz_class *pos2,
                             mpf_class *rate) {
  if (!_intervals) return false;
  if (_interval_index == _intervals->l) {
    if (!load_next_chr()) return false;
  }
  *chr = get_loaded_chr();
  *pos1 = _intervals->start[_interval_index];
  *pos2 = _intervals->end[_interval_index];
  *rate = _intervals->value[_interval_index];
  return true;
}

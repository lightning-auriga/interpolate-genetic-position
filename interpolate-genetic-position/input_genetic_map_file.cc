/*!
 \file input_genetic_map_file.cc
 \brief implementation for input genetic map file classes
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/input_genetic_map_file.h"

namespace igp = interpolate_genetic_position;

igp::base_input_genetic_map_file::base_input_genetic_map_file() {}
igp::base_input_genetic_map_file::~base_input_genetic_map_file() throw() {}

igp::input_genetic_map_file::input_genetic_map_file() {
  _gzinput = NULL;
  _buffer = 0;
  _buffer_size = 1000;
  _ft = UNKNOWN;
  _chr_lower_bound = "";
  _pos_lower_bound = 0;
  _gpos_lower_bound = 0.0;
  _rate_lower_bound = 0.0;
  _chr_upper_bound = "";
  _pos_upper_bound = 0;
  _gpos_upper_bound = 0.0;
  _rate_upper_bound = 0.0;

  _buffer = new char[_buffer_size];
}
igp::input_genetic_map_file::~input_genetic_map_file() throw() { close(); }
void igp::input_genetic_map_file::open(const std::string &filename,
                                       format_type ft) {
  _ft = ft;
  // function must chomp the expected header line of bolt format files
  if (filename.rfind(".gz") == filename.size() - 3) {
    _gzinput = gzopen(filename.c_str(), "rb");
    if (!_gzinput) {
      throw std::runtime_error("input_genetic_map_file: cannot read file \"" +
                               filename + "\"");
    }
    if (_ft == BOLT) {
      gzgets(_gzinput, _buffer, _buffer_size);
    }
  } else {
    std::string line = "";
    _input.open(filename.c_str());
    if (!_input.is_open()) {
      throw std::runtime_error("input_genetic_map_file: cannot read file \"" +
                               filename + "\"");
    }
    if (_ft == BOLT) {
      getline(_input, line);
    }
  }
  // Load the first two values, such that a valid range is available at the
  // beginning of iteration.
  get();
  get();
}
bool igp::input_genetic_map_file::get() {
  _chr_lower_bound = _chr_upper_bound;
  _pos_lower_bound = _pos_upper_bound;
  _gpos_lower_bound = _gpos_upper_bound;
  _rate_lower_bound = _rate_upper_bound;
  std::string line = "", pos1 = "", pos2 = "", gpos = "", rate = "";
  if (_input.is_open()) {
    if (_input.peek() == EOF) {
      return false;
    }
    getline(_input, line);
  } else {
    if (gzgets(_gzinput, _buffer, _buffer_size) == Z_NULL) {
      return false;
    }
    line = std::string(_buffer);
  }
  std::istringstream strm1(line);
  if (_ft == BOLT) {
    // BOLT format includes a first column chromosome indicator,
    // without a "chr" prefix independent of genome build.
    // It includes 1-22 and X encoded as 23.
    if (!(strm1 >> _chr_upper_bound >> pos1 >> rate >> gpos)) {
      throw std::runtime_error(
          "input_genetic_map_file::get: cannot parse BOLT ratefile line \"" +
          line + "\"");
    }
    _pos_upper_bound = pos1;
    _gpos_upper_bound = gpos;
    _rate_upper_bound = rate;
  } else if (_ft == UCSC) {
    // UCSC format includes rate but not genetic position itself. Assuming the
    // first position in a rate file is 0 genetic position, as is conventional,
    // we need to track the accumulated genetic position each time we get a new
    // entry from the file.
    if (!(strm1 >> _chr_upper_bound >> pos1 >> pos2 >> rate)) {
      throw std::runtime_error(
          "input_genetic_map_file::get: cannot parse UCSC bedfile line \"" +
          line + "\"");
    }
    _pos_upper_bound = pos1;
    _rate_upper_bound = rate;
    if (_chr_upper_bound == _chr_lower_bound) {
      _gpos_upper_bound =
          _gpos_lower_bound + _rate_lower_bound *
                                  (_pos_upper_bound - _pos_lower_bound) /
                                  mpf_class(1000000.0);
    } else {
      _gpos_upper_bound = mpf_class("0.0");
    }
  } else if (_ft == UNKNOWN) {
    throw std::runtime_error("input_genetic_map_file::get: format is unset");
  } else {
    throw std::runtime_error(
        "input_genetic_map_file::get: unrecognized file format "
        "(check format_type enum)");
  }
  return true;
}
void igp::input_genetic_map_file::close() {
  if (_input.is_open()) {
    _input.close();
  }
  if (_gzinput) {
    gzclose(_gzinput);
    _gzinput = 0;
  }
  if (_buffer) {
    delete[] _buffer;
    _buffer = 0;
  }
}
bool igp::input_genetic_map_file::eof() {
  if (_input.is_open()) {
    return _input.peek() == EOF;
  }
  if (_gzinput) {
    return gzeof(_gzinput);
  }
  return false;
}
std::string igp::input_genetic_map_file::get_chr_lower_bound() const {
  return _chr_lower_bound;
}
std::string igp::input_genetic_map_file::get_chr_upper_bound() const {
  return _chr_upper_bound;
}
mpz_class igp::input_genetic_map_file::get_pos_lower_bound() const {
  return _pos_lower_bound;
}
mpz_class igp::input_genetic_map_file::get_pos_upper_bound() const {
  return _pos_upper_bound;
}
mpf_class igp::input_genetic_map_file::get_gpos_lower_bound() const {
  return _gpos_lower_bound;
}
mpf_class igp::input_genetic_map_file::get_gpos_upper_bound() const {
  return _gpos_upper_bound;
}
mpf_class igp::input_genetic_map_file::get_rate_lower_bound() const {
  return _rate_lower_bound;
}
mpf_class igp::input_genetic_map_file::get_rate_upper_bound() const {
  return _rate_upper_bound;
}

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
  _startpos_lower_bound = 0;
  _endpos_lower_bound = -1;
  _gpos_lower_bound = 0.0;
  _rate_lower_bound = 0.0;
  _chr_upper_bound = "";
  _startpos_upper_bound = 0;
  _endpos_upper_bound = -1;
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
  } else if (ft == BIGWIG) {
    _bwinput.open(filename);
  } else if (!filename.empty()) {
    std::string line = "";
    _input.open(filename.c_str());
    if (!_input.is_open()) {
      throw std::runtime_error("input_genetic_map_file: cannot read file \"" +
                               filename + "\"");
    }
    if (_ft == BOLT) {
      getline(_input, line);
    }
  } else {  // file is streamed from cin
    if (_ft == BOLT) {
      std::string line = "";
      getline(*get_fallback_stream(), line);
    }
  }
  // Load the first two values, such that a valid range is available at the
  // beginning of iteration.
  get();
  get();
}
void igp::input_genetic_map_file::set_fallback_stream(std::istream *ptr) {
  _fallback = ptr;
}
std::istream *igp::input_genetic_map_file::get_fallback_stream() const {
  if (!_fallback) {
    throw std::runtime_error(
        "igmf::get_fallback_stream: called on NULL pointer");
  }
  return _fallback;
}
bool igp::input_genetic_map_file::get() {
  _chr_lower_bound = _chr_upper_bound;
  _startpos_lower_bound = _startpos_upper_bound;
  _endpos_lower_bound = _endpos_upper_bound;
  _gpos_lower_bound = _gpos_upper_bound;
  _rate_lower_bound = _rate_upper_bound;
  std::string line = "", pos1 = "", gpos = "", rate = "";
  if (_input.is_open()) {
    if (_input.peek() == EOF) {
      return false;
    }
    getline(_input, line);
  } else if (_bwinput.is_open()) {
    mpz_class pos2;
    if (!_bwinput.get(&_chr_upper_bound, &_startpos_upper_bound,
                      &_endpos_upper_bound, &_rate_upper_bound)) {
      return false;
    }
    _startpos_upper_bound = _startpos_upper_bound + 1;
    _endpos_upper_bound = _endpos_upper_bound + 1;
    // as with bedgraph below, requires interpolation
    if (_chr_upper_bound == _chr_lower_bound) {
      _gpos_upper_bound = _gpos_lower_bound +
                          _rate_lower_bound *
                              (_startpos_upper_bound - _startpos_lower_bound) /
                              mpf_class(1000000.0);
    } else {
      _gpos_upper_bound = mpf_class("0.0");
    }
    // return immediately here, as string parsing is handled automatically
    // upstream
    return true;
  } else if (_gzinput != NULL) {
    if (gzgets(_gzinput, _buffer, _buffer_size) == Z_NULL) {
      return false;
    }
    line = std::string(_buffer);
  } else {
    if (get_fallback_stream()->peek() == EOF) {
      return false;
    }
    getline(*get_fallback_stream(), line);
  }
  std::istringstream strm1(line);
  if (_ft == BOLT) {
    // BOLT format includes a first column chromosome indicator,
    // without a "chr" prefix independent of genome build.
    // It includes 1-22 and X encoded as 23.
    if (!(strm1 >> _chr_upper_bound >> _startpos_upper_bound >>
          _rate_upper_bound >> _gpos_upper_bound)) {
      throw std::runtime_error(
          "input_genetic_map_file::get: cannot parse BOLT ratefile line \"" +
          line + "\"");
    }
  } else if (_ft == BEDGRAPH) {
    // bedgraph tracks from UCSC include rate but not genetic position itself.
    // Assuming the first position in a rate file is 0 genetic position, as is
    // conventional, we need to track the accumulated genetic position each
    // time we get a new entry from the file.
    if (!(strm1 >> _chr_upper_bound >> _startpos_upper_bound >>
          _endpos_upper_bound >> _rate_upper_bound)) {
      throw std::runtime_error(
          "input_genetic_map_file::get: cannot parse bedgraph line \"" + line +
          "\"");
    }
    _startpos_upper_bound = _startpos_upper_bound + 1;
    _endpos_upper_bound = _endpos_upper_bound + 1;
    if (_chr_upper_bound == _chr_lower_bound) {
      _gpos_upper_bound = _gpos_lower_bound +
                          _rate_lower_bound *
                              (_startpos_upper_bound - _startpos_lower_bound) /
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
  if (_bwinput.is_open()) {
    return _bwinput.eof();
  }
  return get_fallback_stream()->peek() == EOF;
}
std::string igp::input_genetic_map_file::get_chr_lower_bound() const {
  return _chr_lower_bound;
}
std::string igp::input_genetic_map_file::get_chr_upper_bound() const {
  return _chr_upper_bound;
}
mpz_class igp::input_genetic_map_file::get_startpos_lower_bound() const {
  return _startpos_lower_bound;
}
mpz_class igp::input_genetic_map_file::get_startpos_upper_bound() const {
  return _startpos_upper_bound;
}
mpz_class igp::input_genetic_map_file::get_endpos_lower_bound() const {
  return _endpos_lower_bound;
}
mpz_class igp::input_genetic_map_file::get_endpos_upper_bound() const {
  return _endpos_upper_bound;
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

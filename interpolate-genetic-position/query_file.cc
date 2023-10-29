/*!
 \file query_file.cc
 \brief implementation for query file
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/query_file.h"

namespace igp = interpolate_genetic_position;

igp::query_file::query_file() {
  _gzinput = NULL;
  _buffer = 0;
  _buffer_size = 1000;
  _ft = UNKNOWN;
  _chr_index = 0;
  _pos_index = 0;
  _base0 = false;
  _chr = "";
  _pos = 0;

  _buffer = new char[_buffer_size];
}
igp::query_file::query_file(const query_file &obj) {
  throw std::runtime_error(
      "query_file: copy constructor operation is invalid for this class");
}
igp::query_file::~query_file() throw() { close(); }
void igp::query_file::open(const std::string &filename, format_type ft) {
  _ft = ft;
  if (!filename.empty()) {
    if (filename.rfind(".gz") == filename.size() - 3) {
      _gzinput = gzopen(filename.c_str(), "rb");
      if (!_gzinput) {
        throw std::runtime_error("query_file: cannot read file \"" + filename +
                                 "\"");
      }
    } else {
      _input.open(filename.c_str());
      if (!_input.is_open()) {
        throw std::runtime_error("query_file: cannot read file \"" + filename +
                                 "\"");
      }
    }
  }  // otherwise, input comes from cin
  if (_ft == BIM) {
    _chr_index = 0;
    _pos_index = 3;
    _gpos_index = 2;
    _base0 = false;
    _line_contents.resize(6);
  } else if (_ft == MAP) {
    _chr_index = 0;
    _pos_index = 3;
    _gpos_index = 2;
    _base0 = false;
    _line_contents.resize(4);
  } else if (_ft == BED) {
    _chr_index = 0;
    _pos_index = 1;
    _gpos_index = 4;
    _base0 = true;
    _line_contents.resize(4);
  }
}
void igp::query_file::open(const char *filename, format_type ft) {
  open(std::string(filename), ft);
}
bool igp::query_file::get() {
  std::string line = "", catcher = "";
  if (eof()) {
    return false;
  }
  if (_input.is_open()) {
    getline(_input, line);
  } else if (_gzinput) {
    if (gzgets(_gzinput, _buffer, _buffer_size) == Z_NULL) {
      return false;
    }
    line = std::string(_buffer);
  } else {
    getline(std::cin, line);
  }
  std::istringstream strm1(line);
  for (unsigned i = 0; i < _line_contents.size(); ++i) {
    if (!(strm1 >> _line_contents.at(i))) {
      throw std::runtime_error("query_file: insufficient line tokens");
    }
  }
  _chr = _line_contents.at(_chr_index);
  _pos = _line_contents.at(_pos_index);
  if (_base0) {
    _pos = _pos + 1;
  }
  return true;
}
const std::string &igp::query_file::get_chr() const { return _chr; }
const mpz_class &igp::query_file::get_pos() const { return _pos; }
void igp::query_file::close() {
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
bool igp::query_file::eof() {
  if (_input.is_open()) {
    return _input.peek() == EOF;
  }
  if (_gzinput) {
    return gzeof(_gzinput);
  }
  return std::cin.peek() == EOF;
}
void igp::query_file::report(const mpf_class &gpos_interpolated,
                             std::ostream *output) const {
  for (unsigned i = 0; i < _line_contents.size(); ++i) {
    if (i) {
      if (!(*output << '\t')) {
        throw std::runtime_error("query_file::report: cannot write to stream");
      }
    }
    if ((_ft == BIM || _ft == MAP) && i == 2) {
      if (!(*output << gpos_interpolated)) {
        throw std::runtime_error("query_file::report: cannot write to stream");
      }
    } else {
      if (!(*output << _line_contents.at(i))) {
        throw std::runtime_error("query_file::report: cannot write to stream");
      }
    }
  }
  if (_ft == BED) {
    if (!(*output << '\t' << gpos_interpolated)) {
      throw std::runtime_error("query_file::report: cannot write to stream");
    }
  }
  if (!(*output << '\n')) {
    throw std::runtime_error("query_file::report: cannot write to stream");
  }
}

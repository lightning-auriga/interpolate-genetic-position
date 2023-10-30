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
  throw std::runtime_error("query_file: default constructor disabled");
}
igp::query_file::query_file(const query_file &obj) {
  throw std::runtime_error(
      "query_file: copy constructor operation is invalid for this class");
}
igp::query_file::query_file(base_input_variant_file *ptr)
    : _interface(ptr), _ft(UNKNOWN) {}
igp::query_file::~query_file() throw() {}
void igp::query_file::open(const std::string &filename, format_type ft) {
  _ft = ft;
  _interface->open(filename);
  // handle file formats
  if (_ft == BIM) {
    _interface->set_format_parameters(0, 3, 2, false, 6);
  } else if (_ft == MAP) {
    _interface->set_format_parameters(0, 3, 2, false, 4);
  } else if (_ft == BED) {
    _interface->set_format_parameters(0, 1, 4, true, 4);
  }
}
void igp::query_file::open(const char *filename, format_type ft) {
  open(std::string(filename), ft);
}
bool igp::query_file::get() { return _interface->get_variant(); }
const std::string &igp::query_file::get_chr() const {
  return _interface->get_chr();
}
const mpz_class &igp::query_file::get_pos() const {
  return _interface->get_pos();
}
void igp::query_file::close() { _interface->close(); }
bool igp::query_file::eof() { return _interface->eof(); }
void igp::query_file::report(const mpf_class &gpos_interpolated,
                             std::ostream *output) const {
  for (unsigned i = 0; i < _interface->get_line_contents().size(); ++i) {
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
      if (!(*output << _interface->get_line_contents().at(i))) {
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

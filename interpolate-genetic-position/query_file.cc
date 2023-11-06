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
igp::query_file::query_file(base_input_variant_file *inptr,
                            base_output_variant_file *outptr)
    : _interface(inptr), _ft(UNKNOWN), _output(outptr) {}
igp::query_file::~query_file() throw() {}
void igp::query_file::open(const std::string &filename, format_type ft) {
  _ft = ft;
  _interface->open(filename);
  // handle file formats
  if (_ft == BIM) {
    _interface->set_format_parameters(0, 3, -1, 2, false, 6);
  } else if (_ft == MAP) {
    _interface->set_format_parameters(0, 3, -1, 2, false, 4);
  } else if (_ft == BED) {
    _interface->set_format_parameters(0, 1, 2, 4, true, 4);
  } else if (_ft == SNP) {
    _interface->set_format_parameters(1, 3, -1, 2, false, 4);
  }
}
void igp::query_file::initialize_output(const std::string &filename,
                                        format_type ft) {
  _output->open(filename, ft);
}
bool igp::query_file::get() { return _interface->get_variant(); }
const std::string &igp::query_file::get_chr() const {
  return _interface->get_chr();
}
const mpz_class &igp::query_file::get_pos1() const {
  return _interface->get_pos1();
}
const mpz_class &igp::query_file::get_pos2() const {
  return _interface->get_pos2();
}
void igp::query_file::close() {
  _interface->close();
  _output->close();
}
bool igp::query_file::eof() { return _interface->eof(); }
void igp::query_file::report(const mpf_class &gpos_interpolated) const {
  std::string chr = "", id = "", a1 = "", a2 = "";
  mpz_class pos1 = -1, pos2 = -1;
  if (_ft == BIM || _ft == MAP) {
    chr = _interface->get_line_contents().at(0);
    pos1 = mpz_class(_interface->get_line_contents().at(3));
    id = _interface->get_line_contents().at(1);
  }
  if (_ft == BIM) {
    a1 = _interface->get_line_contents().at(4);
    a2 = _interface->get_line_contents().at(5);
  }
  if (_ft == BED) {
    chr = _interface->get_line_contents().at(0);
    pos1 = mpz_class(_interface->get_line_contents().at(1));
    pos2 = mpz_class(_interface->get_line_contents().at(2));
  }
  _output->write(chr, pos1, pos2, id, gpos_interpolated, a1, a2);
}

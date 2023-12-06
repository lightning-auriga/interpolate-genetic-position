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
    : _interface(inptr),
      _ft(UNKNOWN),
      _output(outptr),
      _index_on_chromosome(0) {}
igp::query_file::~query_file() throw() {}
void igp::query_file::open(const std::string &filename, format_type ft) {
  _ft = ft;
  _interface->open(filename);
  // handle file formats
  if (_ft == BIM) {
    _interface->set_format_parameters(0, 3, -1, 2, -1, false, 6);
  } else if (_ft == MAP) {
    _interface->set_format_parameters(0, 3, -1, 2, -1, false, 4);
  } else if (_ft == BED) {
    _interface->set_format_parameters(0, 1, 2, 4, 3, true, 4);
  } else if (_ft == SNP) {
    _interface->set_format_parameters(1, 3, -1, 2, -1, false, 4);
  } else if (_ft == VCF) {
    // these are ignored due to use of htslib
    _interface->set_format_parameters(0, 1, -1, -1, -1, false, 9);
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
void igp::query_file::report(const std::vector<query_result> &results) {
  std::string id = "", a1 = "", a2 = "";
  if (get_previous_chromosome().compare(results.begin()->get_chr())) {
    set_index_on_chromosome(0);
    set_previous_chromosome(results.begin()->get_chr());
  }
  for (std::vector<query_result>::const_iterator iter = results.begin();
       iter != results.end(); ++iter) {
    if (_ft == BIM || _ft == MAP) {
      id = _interface->get_line_contents().at(1);
    }
    if (_ft == BIM) {
      a1 = _interface->get_line_contents().at(4);
      a2 = _interface->get_line_contents().at(5);
    }
    if (_ft == VCF) {
      id = _interface->get_varid();
      a1 = _interface->get_a1();
      a2 = _interface->get_a2();
    }
    _output->write(
        iter->get_chr(), iter->get_startpos(), iter->get_endpos(), id,
        iter->get_gpos() +
            (_ft == BED ? get_step_interval() * get_index_on_chromosome()
                        : mpf_class("0.0")),
        iter->get_rate(), a1, a2);
  }
  set_index_on_chromosome(get_index_on_chromosome() + 1);
}
const mpf_class &igp::query_file::get_step_interval() const {
  if (_output) {
    return _output->get_step_interval();
  }
  throw std::runtime_error(
      "query_file::get_step_interval: called on unset output");
}
void igp::query_file::set_step_interval(const mpf_class &step) {
  if (_output) {
    _output->set_step_interval(step);
    return;
  }
  throw std::runtime_error(
      "query_file::set_step_interval: called on unset output");
}
unsigned igp::query_file::get_index_on_chromosome() const {
  return _index_on_chromosome;
}
void igp::query_file::set_index_on_chromosome(unsigned index) {
  _index_on_chromosome = index;
}
const std::string &igp::query_file::get_previous_chromosome() const {
  return _previous_chromosome;
}
void igp::query_file::set_previous_chromosome(const std::string &chr) {
  _previous_chromosome = chr;
}

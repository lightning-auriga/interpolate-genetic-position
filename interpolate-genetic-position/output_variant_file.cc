/*!
 \file output_variant_file.cc
 \brief implementation for output variant file base class
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/output_variant_file.h"

namespace igp = interpolate_genetic_position;

igp::base_output_variant_file::base_output_variant_file() {}
igp::base_output_variant_file::~base_output_variant_file() throw() {}

igp::output_variant_file::output_variant_file()
    : igp::base_output_variant_file::base_output_variant_file(),
      _ft(UNKNOWN),
      _output_morgans(false) {}

igp::output_variant_file::~output_variant_file() throw() { close(); }

void igp::output_variant_file::open(const std::string &filename,
                                    format_type ft) {
  // set output format
  _ft = ft;
  // this needs to be updated to catch vcfs
  if (filename.rfind(".gz") == filename.size() - 3) {
    throw std::runtime_error("output gzipped files not yet supported");
  } else if (!filename.empty()) {
    _output.open(filename.c_str());
    if (!_output.is_open()) {
      throw std::runtime_error("output_variant_file: cannot open file \"" +
                               filename + "\"");
    }
  }  // otherwise, filename is empty, and the object will write to std::cout
}

void igp::output_variant_file::close() {
  if (_output.is_open()) {
    _output.close();
    _output.clear();
  }
}

void igp::output_variant_file::write(
    const std::string &chr, const mpz_class &pos1, const mpz_class &pos2,
    const std::string &id, const mpf_class &gpos, const std::string &a1,
    const std::string &a2) {
  // the idea is: format an output line, then emit it to appropriate target
  std::ostringstream out;
  format_type ft = get_format();
  mpf_class output_gpos = output_morgans() ? gpos / mpf_class("100") : gpos;
  if (ft == BIM || ft == MAP) {
    out << chr << '\t' << id << '\t' << output_gpos << '\t'
        << (cmp(pos2, 0) < 0 ? pos1 + 1 : pos1);
    if (ft == BIM) {
      out << '\t' << a1 << '\t' << a2;
    }
  } else if (ft == SNP) {
    out << id << '\t' << chr << '\t' << output_gpos << '\t'
        << (cmp(pos2, 0) < 0 ? pos1 + 1 : pos1);
  } else if (ft == BED) {
    out << chr << '\t';
    if (cmp(pos2, 0) < 0) {
      out << (pos1 - 1) << '\t' << pos1 << '\t';
    } else {
      out << pos1 << '\t' << pos2 << '\t';
    }
    out << output_gpos;
  } else {
    throw std::runtime_error(
        "output_variant_file::write: format not supported");
  }

  if (_output.is_open()) {
    if (!(_output << out.str() << '\n')) {
      throw std::runtime_error(
          "output_variant_file::write: cannot write to file");
    }
  } else {
    std::cout << out.str() << '\n';
  }
}

igp::format_type igp::output_variant_file::get_format() const { return _ft; }

bool igp::output_variant_file::output_morgans() const {
  return _output_morgans;
}

void igp::output_variant_file::output_morgans(bool use_morgans) {
  _output_morgans = use_morgans;
}

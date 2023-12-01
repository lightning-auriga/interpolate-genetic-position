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
      _output_morgans(false),
      _last_chr(""),
      _last_pos1(0),
      _last_pos2(0),
      _last_gpos(0.0),
      _last_rate(0.0),
      _step_interval(0.0) {}

igp::output_variant_file::~output_variant_file() throw() { close(); }

void igp::output_variant_file::open(const std::string &filename,
                                    format_type ft) {
  std::string bolt_header =
      "chr\tposition\tCOMBINED_rate(cM/Mb)\tGenetic_Map(cM)\n";
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
    if (_ft == BOLT) {
      _output << bolt_header;
    }
  } else {  // otherwise, filename is empty, and the object will write to
            // std::cout
    if (_ft == BOLT) {
      std::cout << bolt_header;
    }
  }
}

void igp::output_variant_file::close() {
  if (_output.is_open()) {
    // only for BOLT output: make sure the end of the last chromosome has
    // a placeholder entry with 0 rate
    if (get_format() == BOLT && cmp(get_last_rate(), 0) != 0) {
      std::ostringstream out;
      out << get_last_chr() << '\t' << get_last_pos2() << '\t' << "0\t"
          << (get_last_gpos() +
              get_last_rate() * (get_last_pos2() - get_last_pos1()) /
                  mpf_class(1000000.0) +
              get_step_interval());
      if (_output.is_open()) {
        if (!(_output << out.str() << '\n')) {
          throw std::runtime_error(
              "output_variant_file::write: cannot write to file");
        }
      } else {
        std::cout << out.str() << '\n';
      }
    }
    _output.close();
    _output.clear();
  }
}

void igp::output_variant_file::write(
    const std::string &chr, const mpz_class &pos1, const mpz_class &pos2,
    const std::string &id, const mpf_class &gpos, const mpf_class &rate,
    const std::string &a1, const std::string &a2) {
  // the idea is: format an output line, then emit it to appropriate target
  std::ostringstream out;
  format_type ft = get_format();
  mpf_class output_gpos = output_morgans() ? gpos / mpf_class("100") : gpos;
  mpf_class step_interval = get_step_interval();
  // track when a result is on a different chromosome than the previous ones
  if (get_last_chr().compare(chr)) {
    // for bolt output only, emit dummy results at the end of a chromosome
    if (ft == BOLT && pos2 > 0 && !get_last_chr().empty()) {
      out << get_last_chr() << '\t' << get_last_pos2() << '\t' << "0\t"
          << (get_last_gpos() +
              get_last_rate() * (get_last_pos2() - get_last_pos1()) /
                  mpf_class(1000000.0) +
              step_interval)
          << '\n';
    }
    set_last_chr(chr);
  }
  set_last_pos1(pos1);
  set_last_pos2(pos2);
  set_last_gpos(output_gpos);
  set_last_rate(rate);
  if (ft == BIM || ft == MAP) {
    out << chr << '\t' << id << '\t' << output_gpos << '\t' << pos1;
    if (ft == BIM) {
      out << '\t' << a1 << '\t' << a2;
    }
  } else if (ft == SNP) {
    out << id << '\t' << chr << '\t' << output_gpos << '\t' << pos1;
  } else if (ft == BOLT) {
    out << chr << '\t' << pos1 << '\t' << rate << '\t' << output_gpos;
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

std::string igp::output_variant_file::get_last_chr() const { return _last_chr; }

void igp::output_variant_file::set_last_chr(const std::string &chr) {
  _last_chr = chr;
}

mpf_class igp::output_variant_file::get_last_gpos() const { return _last_gpos; }

void igp::output_variant_file::set_last_gpos(const mpf_class &gpos) {
  _last_gpos = gpos;
}

mpf_class igp::output_variant_file::get_last_rate() const { return _last_rate; }

void igp::output_variant_file::set_last_rate(const mpf_class &rate) {
  _last_rate = rate;
}

mpz_class igp::output_variant_file::get_last_pos1() const { return _last_pos1; }

void igp::output_variant_file::set_last_pos1(const mpz_class &pos1) {
  _last_pos1 = pos1;
}

mpz_class igp::output_variant_file::get_last_pos2() const { return _last_pos2; }

void igp::output_variant_file::set_last_pos2(const mpz_class &pos2) {
  _last_pos2 = pos2;
}

const mpf_class &igp::output_variant_file::get_step_interval() const {
  return _step_interval;
}

void igp::output_variant_file::set_step_interval(
    const mpf_class &step_interval) {
  _step_interval = step_interval;
}

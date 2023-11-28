/*!
 \file utilities.cc
 \brief implementation for global utility functions
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/utilities.h"

namespace igp = interpolate_genetic_position;

igp::format_type igp::string_to_format_type(const std::string &name) {
  if (!name.compare("bolt")) return BOLT;
  if (!name.compare("bedgraph")) return BEDGRAPH;
  if (!name.compare("bigwig")) return BIGWIG;
  if (!name.compare("bim")) return BIM;
  if (!name.compare("map")) return MAP;
  if (!name.compare("bed")) return BED;
  if (!name.compare("snp")) return SNP;
  if (!name.compare("vcf")) return VCF;
  throw std::runtime_error(
      "string_to_format_type: unrecognized type "
      "descriptor: \"" +
      name + "\"");
}
bool igp::chromosome_to_integer(const std::string &chr, int *res) {
  std::string stripped_chr = chr;
  if (stripped_chr.find("chr") == 0) {
    stripped_chr = stripped_chr.substr(3);
  }
  if (!stripped_chr.compare("X")) {
    *res = 23;
  } else if (!stripped_chr.compare("Y")) {
    *res = 24;
  } else if (!stripped_chr.compare("MT") || !stripped_chr.compare("M")) {
    *res = 26;
  } else {
    std::istringstream strm1(stripped_chr);
    if (!(strm1 >> *res)) {
      return false;
    }
  }
  return true;
}
std::string igp::integer_to_chromosome(int chr) {
  if (chr >= 1 && chr <= 22) {
    return "chr" + std::to_string(chr);
  } else if (chr == 23) {
    return "chrX";
  } else if (chr == 24) {
    return "chrY";
  } else if (chr == 26) {
    return "chrM";
  } else {
    throw std::runtime_error("integer_to_chromosome: unknown chromosome");
  }
}
igp::direction igp::chromosome_compare(const std::string &chr1,
                                       const std::string &chr2) {
  int int1 = 0;
  chromosome_to_integer(chr1, &int1);
  int int2 = 0;
  chromosome_to_integer(chr2, &int2);
  if (int1 == int2) return EQUAL;
  return int1 < int2 ? LESS_THAN : GREATER_THAN;
}
std::string igp::next_chromosome(const std::string &current_chr) {
  int chrint = 0;
  chromosome_to_integer(current_chr, &chrint);
  chrint = chrint == 24 ? 26 : chrint + 1;
  return integer_to_chromosome(chrint);
}

igp::query_result::query_result()
    : _chr(""), _startpos(0), _endpos(0), _gpos(0.0), _rate(0.0) {}

igp::query_result::query_result(const query_result &obj)
    : _chr(obj._chr),
      _startpos(obj._startpos),
      _endpos(obj._endpos),
      _gpos(obj._gpos),
      _rate(obj._rate) {}

igp::query_result::~query_result() throw() {}

void igp::query_result::set_chr(const std::string &chr) { _chr = chr; }

const std::string &igp::query_result::get_chr() const { return _chr; }

void igp::query_result::set_startpos(const mpz_class &pos) { _startpos = pos; }

const mpz_class &igp::query_result::get_startpos() const { return _startpos; }

void igp::query_result::set_endpos(const mpz_class &pos) { _endpos = pos; }

const mpz_class &igp::query_result::get_endpos() const { return _endpos; }

void igp::query_result::set_gpos(const mpf_class &gpos) { _gpos = gpos; }

const mpf_class &igp::query_result::get_gpos() const { return _gpos; }

void igp::query_result::set_rate(const mpf_class &rate) { _rate = rate; }

const mpf_class &igp::query_result::get_rate() const { return _rate; }

void igp::check_io_combinations(const std::string &informat_str,
                                const std::string &outformat_str) {
  format_type informat = string_to_format_type(informat_str);
  format_type outformat = string_to_format_type(outformat_str);
  if (informat == SNP || informat == BIM || informat == VCF) {
    if (outformat != MAP && outformat != SNP && outformat != BIM) {
      throw std::domain_error("for input format " + informat_str +
                              ", valid output formats are: bim, map, snp");
    }
  } else if (informat == MAP) {
    if (outformat != MAP) {
      throw std::domain_error("for input format " + informat_str +
                              ", valid output format is map");
    }
  } else if (informat == BED) {
    if (outformat != BED) {
      throw std::domain_error("for input format " + informat_str +
                              ", valid output format is bed");
    }
  } else {
    throw std::domain_error("unrecognized input format: " + informat_str);
  }
}

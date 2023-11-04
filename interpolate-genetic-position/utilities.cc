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
  if (!name.compare("bim")) return BIM;
  if (!name.compare("map")) return MAP;
  if (!name.compare("bed")) return BED;
  throw std::runtime_error(
      "string_to_format_type: unrecognized type "
      "descriptor: \"" +
      name + "\"");
}
int igp::chromosome_to_integer(const std::string &chr) {
  int rep = 0;
  std::string stripped_chr = chr;
  if (stripped_chr.find("chr") == 0) {
    stripped_chr = stripped_chr.substr(3);
  }
  if (!stripped_chr.compare("X")) {
    rep = 23;
  } else if (!stripped_chr.compare("Y")) {
    rep = 24;
  } else if (!stripped_chr.compare("MT") || !stripped_chr.compare("M")) {
    rep = 26;
  } else {
    std::istringstream strm1(stripped_chr);
    if (!(strm1 >> rep)) {
      throw std::runtime_error(
          "chromosome_to_integer: unhandled chromosome code \"" + chr + "\"");
    }
  }
  return rep;
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
  unsigned int1 = chromosome_to_integer(chr1);
  unsigned int2 = chromosome_to_integer(chr2);
  if (int1 == int2) return EQUAL;
  return int1 < int2 ? LESS_THAN : GREATER_THAN;
}
std::string igp::next_chromosome(const std::string &current_chr) {
  int chrint = chromosome_to_integer(current_chr);
  chrint = chrint == 24 ? 26 : chrint + 1;
  return integer_to_chromosome(chrint);
}

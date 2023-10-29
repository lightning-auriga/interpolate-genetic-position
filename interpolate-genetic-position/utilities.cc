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
  if (!name.compare("ucsc")) return UCSC;
  if (!name.compare("bim")) return BIM;
  if (!name.compare("map")) return MAP;
  if (!name.compare("bed")) return BED;
  throw std::runtime_error(
      "string_to_format_type: unrecognized type "
      "descriptor: \"" +
      name + "\"");
}

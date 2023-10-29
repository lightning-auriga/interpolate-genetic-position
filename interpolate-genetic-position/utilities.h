/*!
 \file utilities.h
 \brief general utility functions and definitions
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_UTILITIES_H_
#define INTERPOLATE_GENETIC_POSITION_UTILITIES_H_

namespace interpolate_genetic_position {
typedef enum { UNKNOWN, BOLT, UCSC, BIM, MAP, BED, VCF } format_type;
typedef enum { LESS_THAN, EQUAL, GREATER_THAN } direction;
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_UTILITIES_H_

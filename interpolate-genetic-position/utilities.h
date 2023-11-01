/*!
 \file utilities.h
 \brief general utility functions and definitions
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef INTERPOLATE_GENETIC_POSITION_UTILITIES_H_
#define INTERPOLATE_GENETIC_POSITION_UTILITIES_H_

#include <sstream>
#include <stdexcept>
#include <string>

namespace interpolate_genetic_position {
typedef enum { UNKNOWN, BOLT, UCSC, BIM, MAP, BED, VCF } format_type;
typedef enum { LESS_THAN, EQUAL, GREATER_THAN } direction;
format_type string_to_format_type(const std::string &name);
/*!
 * \brief translate an assortment of chromosome representations
 * into simple integers for sort comparison.
 * \param chr chromosome code, as string, for conversion
 * \return input chromosome, as an integer on {[1,24], 26}
 *
 * This function should handle optional "chr" prefixes depending on input
 * build. Chromosomes 1-22 are labeled as normal; X is interpreted as 23; Y is
 * interpreted as 24; M or MT is interpreted as 26. XY/25, the PAR regions,
 * are not currently handled.
 */
int chromosome_to_integer(const std::string &chr);
/*!
 * \brief compare two chromosome representations to get a consistent
 * directionality of their comparison. Will strip prefixes and convert
 * most standard chromosome aliases.
 * \param chr1 first chromosome code to compare
 * \param chr2 second chromosome code to compare
 * \return the relationship of the first chromosome code to the second
 */
direction chromosome_compare(const std::string &chr1, const std::string &chr2);
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_UTILITIES_H_

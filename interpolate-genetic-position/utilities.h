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
typedef enum { UNKNOWN, BOLT, BEDGRAPH, BIM, MAP, BED, VCF } format_type;
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
 * \brief translate an integer chromosome representation
 * into a simple "chr[1-22XYM]" string.
 * \param chr chromosome as integer in {1-24,26}
 * \return input chromosome, as a string "chr[1-22XYM]"
 *
 * I'm not yet clear on whether the calling logic,
 * interacting with libBigWig, needs to deal with alternate
 * string representations for old genome builds.
 */
std::string integer_to_chromosome(int chr);
/*!
 * \brief compare two chromosome representations to get a consistent
 * directionality of their comparison. Will strip prefixes and convert
 * most standard chromosome aliases.
 * \param chr1 first chromosome code to compare
 * \param chr2 second chromosome code to compare
 * \return the relationship of the first chromosome code to the second
 */
direction chromosome_compare(const std::string &chr1, const std::string &chr2);
/*!
 * \brief represent a chromosome in a way that bigwigs are ok with
 * \param chr a chromosome represented as a string
 * \return the input chromosome but in a representation bigwigs like
 */
std::string make_chr_bigwig_friendly(const std::string &chr);
/*!
 * \brief get next chromosome according to some sort criterion
 * \param current_chr current chromosome
 * \return the chromosome that comes after the current one
 *
 * This logic is required to control the behavior of cycling through
 * a bigwig without other guidance.
 */
std::string next_chromosome(const std::string &current_chr);
}  // namespace interpolate_genetic_position

#endif  // INTERPOLATE_GENETIC_POSITION_UTILITIES_H_

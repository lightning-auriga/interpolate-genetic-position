/*!
 \file utilities_test.cc
 \brief test of global utility functions.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/utilities.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyNumber;
using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

TEST(utilitiesTest, ucscNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("ucsc"), igp::UCSC);
}

TEST(utilitiesTest, boltNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("bolt"), igp::BOLT);
}

TEST(utilitiesTest, mapNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("map"), igp::MAP);
}

TEST(utilitiesTest, bimNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("bim"), igp::BIM);
}

TEST(utilitiesTest, bedNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("bed"), igp::BED);
}

TEST(utilitiesTest, chromosomeToIntegerAutosome) {
  EXPECT_EQ(igp::chromosome_to_integer("2"), 2);
}

TEST(utilitiesTest, chromosomeToIntegerXchr) {
  EXPECT_EQ(igp::chromosome_to_integer("X"), 23);
}

TEST(utilitiesTest, chromosomeToIntegerYchr) {
  EXPECT_EQ(igp::chromosome_to_integer("Y"), 24);
}

TEST(utilitiesTest, chromosomeToIntegerMitochrondrion) {
  EXPECT_EQ(igp::chromosome_to_integer("M"), 26);
}

TEST(utilitiesTest, chromosomeToIntegerMitochondrion_mt) {
  EXPECT_EQ(igp::chromosome_to_integer("MT"), 26);
}

TEST(utilitiesTest, chromosomeToIntegerAutosomeChrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chr3"), 3);
}

TEST(utilitiesTest, chromosomeToIntegerXchrChrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chrX"), 23);
}

TEST(utilitiesTest, chromosomeToIntegerYchrChrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chrY"), 24);
}

TEST(utilitiesTest, chromosomeToIntegerMitochondrionChrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chrM"), 26);
}

TEST(utilitiesTest, errorOnInvalidChromosome) {
  EXPECT_THROW(igp::chromosome_to_integer("XY"), std::runtime_error);
}

TEST(utilitiesTest, chromosomeCompareLessthan) {
  EXPECT_EQ(igp::chromosome_compare("chr4", "chrX"), igp::LESS_THAN);
}

TEST(utilitiesTest, chromosomeCompareGreaterthan) {
  EXPECT_EQ(igp::chromosome_compare("chrM", "chrY"), igp::GREATER_THAN);
}

TEST(utilitiesTest, chromosomeCompareEqualto) {
  EXPECT_EQ(igp::chromosome_compare("chrX", "X"), igp::EQUAL);
}

TEST(utilitiesTest, nextChromosomeStandardAutosome) {
  EXPECT_EQ(igp::next_chromosome("chr4"), "chr5");
}

TEST(utilitiesTest, nextChromosome22) {
  EXPECT_EQ(igp::next_chromosome("chr22"), "chrX");
}

TEST(utilitiesTest, nextChromosomeX) {
  EXPECT_EQ(igp::next_chromosome("chrX"), "chrY");
}

TEST(utilitiesTest, nextChromosomeY) {
  EXPECT_EQ(igp::next_chromosome("chrY"), "chrM");
}

TEST(utilitiesTest, nextChromosomeM) {
  EXPECT_THROW(igp::next_chromosome("chrM"), std::runtime_error);
}

TEST(utilitiesTest, integerToChromosomeStandardAutosome) {
  EXPECT_EQ(igp::integer_to_chromosome(10), "chr10");
}

TEST(utilitiesTest, integerToChromosome23) {
  EXPECT_EQ(igp::integer_to_chromosome(23), "chrX");
}

TEST(utilitiesTest, integerToChromosome24) {
  EXPECT_EQ(igp::integer_to_chromosome(24), "chrY");
}

TEST(utilitiesTest, integerToChromosome25) {
  EXPECT_THROW(igp::integer_to_chromosome(25), std::runtime_error);
}

TEST(utilitiesTest, integerToChromosome26) {
  EXPECT_EQ(igp::integer_to_chromosome(26), "chrM");
}

TEST(utilitiesTest, integerToChromosome27) {
  EXPECT_THROW(igp::integer_to_chromosome(27), std::runtime_error);
}

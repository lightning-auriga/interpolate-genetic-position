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

TEST(utilitiesTest, ucsc_name_conversion) {
  EXPECT_EQ(igp::string_to_format_type("ucsc"), igp::UCSC);
}

TEST(utilitiesTest, bolt_name_conversion) {
  EXPECT_EQ(igp::string_to_format_type("bolt"), igp::BOLT);
}

TEST(utilitiesTest, map_name_conversion) {
  EXPECT_EQ(igp::string_to_format_type("map"), igp::MAP);
}

TEST(utilitiesTest, bim_name_conversion) {
  EXPECT_EQ(igp::string_to_format_type("bim"), igp::BIM);
}

TEST(utilitiesTest, bed_name_conversion) {
  EXPECT_EQ(igp::string_to_format_type("bed"), igp::BED);
}

TEST(utilitiesTest, chromosome_to_integer_autosome) {
  EXPECT_EQ(igp::chromosome_to_integer("2"), 2);
}

TEST(utilitiesTest, chromosome_to_integer_xchr) {
  EXPECT_EQ(igp::chromosome_to_integer("X"), 23);
}

TEST(utilitiesTest, chromosome_to_integer_ychr) {
  EXPECT_EQ(igp::chromosome_to_integer("Y"), 24);
}

TEST(utilitiesTest, chromosome_to_integer_mitochrondrion) {
  EXPECT_EQ(igp::chromosome_to_integer("M"), 26);
}

TEST(utilitiesTest, chromosome_to_integer_mitochondrion_mt) {
  EXPECT_EQ(igp::chromosome_to_integer("MT"), 26);
}

TEST(utilitiesTest, chromosome_to_integer_autosome_chrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chr3"), 3);
}

TEST(utilitiesTest, chromosome_to_integer_xchr_chrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chrX"), 23);
}

TEST(utilitiesTest, chromosome_to_integer_ychr_chrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chrY"), 24);
}

TEST(utilitiesTest, chromosome_to_integer_mitochondrion_chrprefix) {
  EXPECT_EQ(igp::chromosome_to_integer("chrM"), 26);
}

TEST(utilitiesTest, error_on_invalid_chromosome) {
  EXPECT_THROW(igp::chromosome_to_integer("XY"), std::runtime_error);
}

TEST(utilitiesTest, chromosome_compare_lessthan) {
  EXPECT_EQ(igp::chromosome_compare("chr4", "chrX"), igp::LESS_THAN);
}

TEST(utilitiesTest, chromosome_compare_greaterthan) {
  EXPECT_EQ(igp::chromosome_compare("chrM", "chrY"), igp::GREATER_THAN);
}

TEST(utilitiesTest, chromosome_compare_equalto) {
  EXPECT_EQ(igp::chromosome_compare("chrX", "X"), igp::EQUAL);
}

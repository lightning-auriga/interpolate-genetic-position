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

TEST(utilitiesTest, bedgraphNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("bedgraph"), igp::BEDGRAPH);
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

TEST(utilitiesTest, snpNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("snp"), igp::SNP);
}

TEST(utilitiesTest, vcfNameConversion) {
  EXPECT_EQ(igp::string_to_format_type("vcf"), igp::VCF);
}

TEST(utilitiesTest, chromosomeToIntegerAutosome) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("2", &chrint));
  EXPECT_EQ(chrint, 2);
}

TEST(utilitiesTest, chromosomeToIntegerXchr) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("X", &chrint));
  EXPECT_EQ(chrint, 23);
}

TEST(utilitiesTest, chromosomeToIntegerYchr) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("Y", &chrint));
  EXPECT_EQ(chrint, 24);
}

TEST(utilitiesTest, chromosomeToIntegerMitochrondrion) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("M", &chrint));
  EXPECT_EQ(chrint, 26);
}

TEST(utilitiesTest, chromosomeToIntegerMitochondrion_mt) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("MT", &chrint));
  EXPECT_EQ(chrint, 26);
}

TEST(utilitiesTest, chromosomeToIntegerAutosomeChrprefix) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("chr3", &chrint));
  EXPECT_EQ(chrint, 3);
}

TEST(utilitiesTest, chromosomeToIntegerXchrChrprefix) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("chrX", &chrint));
  EXPECT_EQ(chrint, 23);
}

TEST(utilitiesTest, chromosomeToIntegerYchrChrprefix) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("chrY", &chrint));
  EXPECT_EQ(chrint, 24);
}

TEST(utilitiesTest, chromosomeToIntegerMitochondrionChrprefix) {
  int chrint = 0;
  EXPECT_TRUE(igp::chromosome_to_integer("chrM", &chrint));
  EXPECT_EQ(chrint, 26);
}

TEST(utilitiesTest, errorOnInvalidChromosome) {
  int chrint = 0;
  EXPECT_FALSE(igp::chromosome_to_integer("XY", &chrint));
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

TEST(utilitiesTest, detectInvalidFormatCombinations) {
  EXPECT_THROW(igp::check_io_combinations("bed", "bim"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("bed", "bolt"));
  EXPECT_THROW(igp::check_io_combinations("bed", "map"), std::domain_error);
  EXPECT_THROW(igp::check_io_combinations("bed", "snp"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("bim", "bim"));
  EXPECT_THROW(igp::check_io_combinations("bim", "bolt"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("bim", "map"));
  EXPECT_NO_THROW(igp::check_io_combinations("bim", "snp"));
  EXPECT_THROW(igp::check_io_combinations("map", "bim"), std::domain_error);
  EXPECT_THROW(igp::check_io_combinations("map", "bolt"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("map", "map"));
  EXPECT_THROW(igp::check_io_combinations("map", "snp"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("snp", "bim"));
  EXPECT_THROW(igp::check_io_combinations("snp", "bolt"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("snp", "map"));
  EXPECT_NO_THROW(igp::check_io_combinations("snp", "snp"));
  EXPECT_NO_THROW(igp::check_io_combinations("vcf", "bim"));
  EXPECT_THROW(igp::check_io_combinations("vcf", "bolt"), std::domain_error);
  EXPECT_NO_THROW(igp::check_io_combinations("vcf", "map"));
  EXPECT_NO_THROW(igp::check_io_combinations("vcf", "snp"));
}

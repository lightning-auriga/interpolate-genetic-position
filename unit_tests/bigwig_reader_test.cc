/*!
 \file bigwig_reader_test.cc
 \brief test of bigwig_reader class.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/bigwig_reader.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/bigwig_reader_test.h"

using ::testing::AnyNumber;
using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

bigwigReaderTest::bigwigReaderTest()
    : testing::Test(), _mpf_error_tolerance(FLT_EPSILON) {}

bigwigReaderTest::~bigwigReaderTest() {}

TEST_F(bigwigReaderTest, flagsMissingInputFile) {
  igp::bigwig_reader bw;
  EXPECT_THROW(bw.open("dummyfile.txt"), std::runtime_error);
}

TEST_F(bigwigReaderTest, flagsInvalidInputFile) {
  igp::bigwig_reader bw;
  EXPECT_THROW(bw.open("README.md"), std::runtime_error);
}

TEST_F(bigwigReaderTest, canProcessTestBigwig) {
  igp::bigwig_reader bw;
  std::string chr = "";
  mpz_class pos1, pos2;
  mpf_class rate;
  EXPECT_NO_THROW(bw.open("unit_tests/test.bw"));
  EXPECT_EQ(bw.get_loaded_chr(), "chr1");
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr1");
  EXPECT_EQ(pos1, mpz_class(1431813));
  EXPECT_EQ(pos2, mpz_class(1515567));
  EXPECT_EQ(cmp(abs(rate - mpf_class("0.0324492")), _mpf_error_tolerance), -1);
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr1");
  EXPECT_EQ(pos1, mpz_class(1515567));
  EXPECT_EQ(pos2, mpz_class(1530002));
  EXPECT_EQ(cmp(abs(rate - mpf_class("0.189597")), _mpf_error_tolerance), -1);
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr2");
  EXPECT_EQ(pos1, mpz_class(522921));
  EXPECT_EQ(pos2, mpz_class(535043));
  EXPECT_EQ(cmp(abs(rate - mpf_class("0.16944")), _mpf_error_tolerance), -1);
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr2");
  EXPECT_EQ(pos1, mpz_class(535043));
  EXPECT_EQ(pos2, mpz_class(548121));
  EXPECT_EQ(cmp(abs(rate - mpf_class("0.0914548")), _mpf_error_tolerance), -1);
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr10");
  EXPECT_EQ(pos1, mpz_class(589917));
  EXPECT_EQ(pos2, mpz_class(616253));
  EXPECT_EQ(rate, mpf_class("0.0"));
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr10");
  EXPECT_EQ(pos1, mpz_class(616253));
  EXPECT_EQ(pos2, mpz_class(634488));
  EXPECT_EQ(cmp(abs(rate - mpf_class("3.01017e-39")), _mpf_error_tolerance),
            -1);
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr22");
  EXPECT_EQ(pos1, mpz_class(17076254));
  EXPECT_EQ(pos2, mpz_class(17081023));
  EXPECT_EQ(rate, mpf_class("0.0"));
  EXPECT_FALSE(bw.eof());
  EXPECT_TRUE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_EQ(chr, "chr22");
  EXPECT_EQ(pos1, mpz_class(17081023));
  EXPECT_EQ(pos2, mpz_class(17083846));
  EXPECT_EQ(cmp(abs(rate - mpf_class("0.137252")), _mpf_error_tolerance), -1);
  EXPECT_FALSE(bw.eof());
  EXPECT_FALSE(bw.get(&chr, &pos1, &pos2, &rate));
  EXPECT_TRUE(bw.eof());
  EXPECT_NO_THROW(bw.close());
}

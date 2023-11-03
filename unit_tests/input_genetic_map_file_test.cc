/*!
 \file input_genetic_map_file_test.cc
 \brief test of interface between genetic_map and recombination map files.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "unit_tests/input_genetic_map_file_test.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyNumber;
using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

igp::mock_input_genetic_map_file::mock_input_genetic_map_file()
    : igp::base_input_genetic_map_file() {}
igp::mock_input_genetic_map_file::~mock_input_genetic_map_file() throw() {}

inputGeneticMapFileTest::inputGeneticMapFileTest()
    : testing::Test(),
      _tmpfile(boost::filesystem::unique_path().native()),
      _mpf_error_tolerance("0.0000000000001") {}

inputGeneticMapFileTest::~inputGeneticMapFileTest() {
  if (boost::filesystem::exists(_tmpfile)) {
    boost::filesystem::remove(_tmpfile);
  }
}

TEST_F(inputGeneticMapFileTest, geneticMapUnderstandsBedgraph) {
  // write example bedgraph recombination data
  std::ofstream output;
  output.open(_tmpfile.c_str());
  ASSERT_TRUE(output.is_open());
  ASSERT_TRUE(
      output << "chr1\t0\t1\t0.012\nchr1\t1\t2\t0.111\nchr4\t2\t3\t0.213"
             << std::endl);
  output.close();

  igp::input_genetic_map_file mapfile;
  EXPECT_NO_THROW(mapfile.open(_tmpfile, igp::UCSC));
  EXPECT_FALSE(mapfile.eof());
  EXPECT_EQ(mapfile.get_chr_lower_bound(), "chr1");
  EXPECT_EQ(mapfile.get_chr_upper_bound(), "chr1");
  EXPECT_EQ(mapfile.get_pos_lower_bound(), mpz_class(0));
  EXPECT_EQ(mapfile.get_pos_upper_bound(), mpz_class(1));
  EXPECT_EQ(mapfile.get_gpos_lower_bound(), mpf_class("0.0"));
  EXPECT_EQ(cmp(abs(mapfile.get_gpos_upper_bound() - mpf_class("0.000000012")),
                _mpf_error_tolerance),
            -1);
  EXPECT_EQ(mapfile.get_rate_lower_bound(), mpf_class("0.012"));
  EXPECT_EQ(mapfile.get_rate_upper_bound(), mpf_class("0.111"));
  EXPECT_TRUE(mapfile.get());
  EXPECT_EQ(mapfile.get_chr_lower_bound(), "chr1");
  EXPECT_EQ(mapfile.get_chr_upper_bound(), "chr4");
  EXPECT_EQ(mapfile.get_pos_lower_bound(), mpz_class(1));
  EXPECT_EQ(mapfile.get_pos_upper_bound(), mpz_class(2));
  EXPECT_EQ(cmp(abs(mapfile.get_gpos_lower_bound() - mpf_class("0.000000012")),
                _mpf_error_tolerance),
            -1);
  EXPECT_EQ(mapfile.get_gpos_upper_bound(), mpf_class("0.0"));
  EXPECT_EQ(mapfile.get_rate_lower_bound(), mpf_class("0.111"));
  EXPECT_EQ(mapfile.get_rate_upper_bound(), mpf_class("0.213"));
  EXPECT_TRUE(mapfile.eof());
  EXPECT_FALSE(mapfile.get());
  EXPECT_NO_THROW(mapfile.close());
}

TEST_F(inputGeneticMapFileTest, geneticMapUnderstandsBolt) {
  // write example bolt format recombination data
  std::ofstream output;
  output.open(_tmpfile.c_str());
  ASSERT_TRUE(output.is_open());
  ASSERT_TRUE(output << "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
                     << "chr1 1 0.012 0\n"
                     << "chr1 2 0.111 0.000000012\n"
                     << "chr4 3 0.213 0" << std::endl);
  output.close();

  igp::input_genetic_map_file mapfile;
  EXPECT_NO_THROW(mapfile.open(_tmpfile, igp::BOLT));
  EXPECT_FALSE(mapfile.eof());
  EXPECT_EQ(mapfile.get_chr_lower_bound(), "chr1");
  EXPECT_EQ(mapfile.get_chr_upper_bound(), "chr1");
  EXPECT_EQ(mapfile.get_pos_lower_bound(), mpz_class(1));
  EXPECT_EQ(mapfile.get_pos_upper_bound(), mpz_class(2));
  EXPECT_EQ(mapfile.get_gpos_lower_bound(), mpf_class("0.0"));
  EXPECT_EQ(cmp(abs(mapfile.get_gpos_upper_bound() - mpf_class("0.000000012")),
                _mpf_error_tolerance),
            -1);
  EXPECT_EQ(mapfile.get_rate_lower_bound(), mpf_class("0.012"));
  EXPECT_EQ(mapfile.get_rate_upper_bound(), mpf_class("0.111"));
  EXPECT_TRUE(mapfile.get());
  EXPECT_EQ(mapfile.get_chr_lower_bound(), "chr1");
  EXPECT_EQ(mapfile.get_chr_upper_bound(), "chr4");
  EXPECT_EQ(mapfile.get_pos_lower_bound(), mpz_class(2));
  EXPECT_EQ(mapfile.get_pos_upper_bound(), mpz_class(3));
  EXPECT_EQ(mapfile.get_gpos_lower_bound(), mpf_class("0.000000012"));
  EXPECT_EQ(mapfile.get_gpos_upper_bound(), mpf_class("0.0"));
  EXPECT_EQ(mapfile.get_rate_lower_bound(), mpf_class("0.111"));
  EXPECT_EQ(mapfile.get_rate_upper_bound(), mpf_class("0.213"));
  EXPECT_TRUE(mapfile.eof());
  EXPECT_FALSE(mapfile.get());
  EXPECT_NO_THROW(mapfile.close());
}

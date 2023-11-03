/*!
 \file input_genetic_map_file_test.cc
 \brief test of interface between genetic_map and recombination map files.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "unit_tests/input_genetic_map_file_test.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "interpolate-genetic-position/genetic_map.h"

using ::testing::AnyNumber;
using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;
mpf_class mpf_error_tolerance("0.000000001");

igp::mock_input_genetic_map_file::mock_input_genetic_map_file()
    : igp::base_input_genetic_map_file() {}
igp::mock_input_genetic_map_file::~mock_input_genetic_map_file() throw() {}

TEST(input_genetic_map_fileTest, open_command_executed) {
  igp::mock_input_genetic_map_file mockfile;
  std::string filename = "test.ucsc.bedgraph";
  EXPECT_CALL(mockfile, open(filename, igp::UCSC)).Times(1).WillOnce(Return());
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  igp::genetic_map gm(&mockfile);
  gm.open(filename, igp::UCSC);
}

TEST(input_genetic_map_fileTest, basic_constructor_disabled) {
  EXPECT_THROW(igp::genetic_map gm, std::runtime_error);
}

TEST(input_genetic_map_fileTest, copy_constructor_disabled) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_THROW(igp::genetic_map gm2(gm), std::runtime_error);
}

TEST(input_genetic_map_fileTest, query_beyond_end_of_chromosome) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(23.4)));
  mpf_class gpos_interpolated;
  std::string query_chr = "1";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  EXPECT_EQ(gpos_interpolated, mpf_class(23.4));
}

TEST(input_genetic_map_fileTest, query_beyond_end_of_file) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(2).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("23"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("23"));
  mpf_class gpos_interpolated;
  std::string query_chr = "M";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  EXPECT_EQ(gpos_interpolated, mpf_class(0.0));
}

TEST(input_genetic_map_fileTest, query_interior_uncovered_chromosome) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("3"));
  mpf_class gpos_interpolated;
  std::string query_chr = "2";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  EXPECT_EQ(gpos_interpolated, mpf_class(0.0));
}

TEST(input_genetic_map_fileTest, query_beginning_of_chromosome) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_pos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_pos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  mpf_class gpos_interpolated;
  std::string query_chr = "1";
  mpz_class query_pos = 50000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  EXPECT_EQ(gpos_interpolated, mpf_class(0.0));
}

TEST(input_genetic_map_fileTest, query_standard_interpolation) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_pos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_pos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(0.001)));
  EXPECT_CALL(mockfile, get_rate_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(0.05)));
  mpf_class gpos_interpolated;
  std::string query_chr = "1";
  mpz_class query_pos = 150000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  mpf_class gpos_expected("0.0035");
  EXPECT_EQ(cmp(abs(gpos_interpolated - gpos_expected), mpf_error_tolerance),
            -1);
}

TEST(input_genetic_map_fileTest, query_overlap_lower_bound) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_pos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_pos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("0.001")));
  mpf_class gpos_interpolated;
  std::string query_chr = "1";
  mpz_class query_pos = 100000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  EXPECT_EQ(gpos_interpolated, mpf_class("0.001"));
}

TEST(input_genetic_map_fileTest, query_overlap_upper_bound) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_pos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_pos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_gpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("0.002")));
  mpf_class gpos_interpolated;
  std::string query_chr = "1";
  mpz_class query_pos = 200000;
  bool verbose = false;
  gm.query(query_chr, query_pos, verbose, &gpos_interpolated);
  EXPECT_EQ(gpos_interpolated, mpf_class("0.002"));
}

TEST(input_genetic_map_fileTest, query_detects_unsorted_chromosomes) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("4"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("3"));
  mpf_class gpos_interpolated;
  std::string query_chr = "2";
  mpz_class query_pos = 100000;
  bool verbose = false;
  EXPECT_THROW(gm.query(query_chr, query_pos, verbose, &gpos_interpolated),
               std::runtime_error);
}

TEST(input_genetic_map_fileTest, query_detects_unsorted_positions) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_pos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_pos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(100000)));
  mpf_class gpos_interpolated;
  std::string query_chr = "1";
  mpz_class query_pos = 5000;
  bool verbose = false;
  EXPECT_THROW(gm.query(query_chr, query_pos, verbose, &gpos_interpolated),
               std::runtime_error);
}

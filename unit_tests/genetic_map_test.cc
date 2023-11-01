/*!
 \file genetic_map_test.cc
 \brief test of operations of genetic_map.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/genetic_map.h"

#include "gtest/gtest.h"

using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

TEST(genetic_mapTest, query_chromosome_less_than_window) {}

TEST(genetic_mapTest, query_chromosome_greater_than_window) {}

TEST(genetic_mapTest, query_beyond_lower_bound_before_upper_chromosome) {}

TEST(genetic_mapTest, query_chromosome_has_no_map_internal) {}

TEST(genetic_mapTest, query_chromosome_has_no_map_before_range) {}

TEST(genetic_mapTest, query_chromosome_has_no_map_after_range) {}

TEST(input_variant_fileTest, can_initialize_bimfile) {
  igp::mock_input_variant_file mockfile;
  std::string bimfilename = "test.bim";
  EXPECT_CALL(mockfile, set_format_parameters(0, 3, 2, false, 6))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(bimfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(bimfilename, igp::BIM);
}

TEST(input_variant_fileTest, can_initialize_mapfile) {
  igp::mock_input_variant_file mockfile;
  std::string mapfilename = "test.map";
  EXPECT_CALL(mockfile, set_format_parameters(0, 3, 2, false, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(mapfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(mapfilename, igp::MAP);
}

TEST(input_variant_fileTest, can_initialize_bedfile) {
  igp::mock_input_variant_file mockfile;
  std::string bedfilename = "test.bed";
  EXPECT_CALL(mockfile, set_format_parameters(0, 1, 4, true, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(bedfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(bedfilename, igp::BED);
}

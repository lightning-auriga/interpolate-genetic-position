/*!
 \file input_variant_file_test.cc
 \brief test of interface between query_file and input variant files.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "unit_tests/input_variant_file_test.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "interpolate-genetic-position/query_file.h"

using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

igp::mock_input_variant_file::mock_input_variant_file()
    : igp::base_input_variant_file() {}
igp::mock_input_variant_file::~mock_input_variant_file() throw() {}

TEST(inputVariantFileTest, can_initialize_bimfile) {
  igp::mock_input_variant_file mockfile;
  std::string bimfilename = "test.bim";
  EXPECT_CALL(mockfile, set_format_parameters(0, 3, 2, false, 6))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(bimfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(bimfilename, igp::BIM);
}

TEST(inputVariantFileTest, can_initialize_mapfile) {
  igp::mock_input_variant_file mockfile;
  std::string mapfilename = "test.map";
  EXPECT_CALL(mockfile, set_format_parameters(0, 3, 2, false, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(mapfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(mapfilename, igp::MAP);
}

TEST(inputVariantFileTest, can_initialize_bedfile) {
  igp::mock_input_variant_file mockfile;
  std::string bedfilename = "test.bed";
  EXPECT_CALL(mockfile, set_format_parameters(0, 1, 4, true, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(bedfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(bedfilename, igp::BED);
}

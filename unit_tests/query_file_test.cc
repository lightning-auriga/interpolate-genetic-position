/*!
 \file query_file_test.cc
 \brief test of variant query file operations.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/query_file.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/input_variant_file_test.h"
#include "unit_tests/output_variant_file_test.h"

using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

TEST(queryFileTest, canInitializeBimfile) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string bimfilename = "test.bim";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 3, 2, false, 6))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(bimfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(bimfilename, igp::BIM);
}

TEST(queryFileTest, canInitializeMapfile) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string mapfilename = "test.map";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 3, 2, false, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(mapfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(mapfilename, igp::MAP);
}

TEST(queryFileTest, canInitializeBedfile) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string bedfilename = "test.bed";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 1, 4, true, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(bedfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(bedfilename, igp::BED);
}

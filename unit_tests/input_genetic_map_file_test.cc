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

using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

igp::mock_input_genetic_map_file::mock_input_genetic_map_file()
    : igp::base_input_genetic_map_file() {}
igp::mock_input_genetic_map_file::~mock_input_genetic_map_file() throw() {}

/*
TEST(input_genetic_map_fileTest, can_initialize_bimfile) {
  igp::mock_input_genetic_map_file mockfile;
  std::string bimfilename = "test.bim";
  EXPECT_CALL(mockfile, set_format_parameters(0, 3, 2, false, 6))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, open(bimfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mockfile);
  qf.open(bimfilename, igp::BIM);
}
*/

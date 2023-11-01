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

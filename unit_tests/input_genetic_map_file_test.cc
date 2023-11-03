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

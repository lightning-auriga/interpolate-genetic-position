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

TEST(input_variant_fileTest, can_load_bimfile) {}

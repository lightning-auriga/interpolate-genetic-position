/*!
 \file output_variant_file_test.cc
 \brief test of interface between query_file and output variant files.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "unit_tests/output_variant_file_test.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

igp::mock_output_variant_file::mock_output_variant_file()
    : igp::base_output_variant_file() {}
igp::mock_output_variant_file::~mock_output_variant_file() throw() {}

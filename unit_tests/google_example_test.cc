/*!
 \file google_example_test.cc
 \brief example google unit test.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Docs at <https://github.com/google/googletest>
 */

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AtLeast;
using ::testing::Return;

TEST(basicTest, can_do_numbers) { EXPECT_EQ(1, 1); }

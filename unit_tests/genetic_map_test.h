/*!
 \file genetic_map_test.h
 \brief declaration of test suite for genetic map class
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Powered by gmock from
 <https://github.com/google/googletest/tree/main/googlemock>
 */

#ifndef UNIT_TESTS_GENETIC_MAP_TEST_H_
#define UNIT_TESTS_GENETIC_MAP_TEST_H_

#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "interpolate-genetic-position/genetic_map.h"

class geneticMapTest : public testing::Test {
 protected:
  geneticMapTest();
  ~geneticMapTest() throw();
  mpf_class _mpf_error_tolerance;
};
#endif  // UNIT_TESTS_GENETIC_MAP_TEST_H_

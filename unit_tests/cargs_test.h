/*!
 \file cargs_test.h
 \brief declaration of test suite for command line interface
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef UNIT_TESTS_CARGS_TEST_H_
#define UNIT_TESTS_CARGS_TEST_H_

#include <cfloat>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "interpolate-genetic-position/cargs.h"

class cargsTest : public testing::Test {
 protected:
  cargsTest();
  ~cargsTest() throw();
  void populate(const std::string &str, std::vector<std::string> *vec,
                const char ***arr);
  std::vector<std::string> _argvec1;
  const char **_argv1;
};
#endif  // UNIT_TESTS_CARGS_TEST_H_

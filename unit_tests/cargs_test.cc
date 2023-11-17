/*!
 \file cargs_test.cc
 \brief test of command line interface.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/cargs.h"

#include "gtest/gtest.h"
#include "unit_tests/cargs_test.h"

namespace igp = interpolate_genetic_position;

cargsTest::cargsTest() : testing::Test(), _argv1(NULL) {
  std::string test1 = "progname -h";
  populate(test1, &_argvec1, &_argv1);
}

cargsTest::~cargsTest() {
  if (_argv1) {
    delete[] _argv1;
  }
}

void cargsTest::populate(const std::string &str, std::vector<std::string> *vec,
                         const char ***arr) {
  if (!vec || !arr) {
    throw std::runtime_error("populate: null pointer(s)");
  }
  vec->clear();
  std::istringstream strm1(str);
  std::string token = "";
  while (strm1 >> token) {
    vec->push_back(token);
  }
  if (*arr) {
    delete[] *arr;
    *arr = 0;
  }
  *arr = new const char *[vec->size()];
  for (unsigned i = 0; i < vec->size(); ++i) {
    (*arr)[i] = vec->at(i).c_str();
  }
}

TEST_F(cargsTest, basicInitialization) {
  igp::cargs ap(_argvec1.size(), _argv1);
}

TEST_F(cargsTest, copyConstructor) {
  igp::cargs ap1(_argvec1.size(), _argv1);
  igp::cargs ap2(ap1);
}

TEST_F(cargsTest, printHelp) {
  std::ostringstream o;
  igp::cargs ap(_argvec1.size(), _argv1);
  ap.print_help(o);
  EXPECT_EQ(o.str().find("Recognized options:"), 0UL);
}

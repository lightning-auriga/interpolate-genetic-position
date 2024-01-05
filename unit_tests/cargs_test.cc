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

cargsTest::cargsTest()
    : testing::Test(),
      _argv1(NULL),
      _argv2(NULL),
      _argv3(NULL),
      _argv4(NULL),
      _argv5(NULL) {
  std::string test1 = "progname -h";
  populate(test1, &_argvec1, &_argv1);
  std::string test2 =
      "progname -i fn1 -p bed -g fn2 -m bolt -o fn3 --output-morgans "
      "--region-step-interval 1.23 -f map";
  populate(test2, &_argvec2, &_argv2);
  std::string test3 = "progname -i fn1 -p bim -g fn2 -m bedgraph -o fn3 -v";
  populate(test3, &_argvec3, &_argv3);
  std::string test4 = "progname -i fn1 -p bedgraph -g fn2 -m map -o fn3 -f vcf";
  populate(test4, &_argvec4, &_argv4);
  std::string test5 = "progname --version --precision 68";
  populate(test5, &_argvec5, &_argv5);
}

cargsTest::~cargsTest() {
  if (_argv1) {
    delete[] _argv1;
  }
  if (_argv2) {
    delete[] _argv2;
  }
  if (_argv3) {
    delete[] _argv3;
  }
  if (_argv4) {
    delete[] _argv4;
  }
  if (_argv5) {
    delete[] _argv5;
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

TEST_F(cargsTest, printVersion) {
  std::ostringstream o;
  igp::cargs ap(_argvec5.size(), _argv5);
  ap.print_version(o);
  EXPECT_EQ(o.str().find("interpolate-genetic-position"), 0UL);
}

TEST_F(cargsTest, detectHelp) {
  igp::cargs ap1(_argvec1.size(), _argv1);
  EXPECT_TRUE(ap1.help());
  igp::cargs ap2(_argvec2.size(), _argv2);
  EXPECT_FALSE(ap2.help());
}

TEST_F(cargsTest, detectVerbose) {
  igp::cargs ap1(_argvec3.size(), _argv3);
  EXPECT_TRUE(ap1.verbose());
  igp::cargs ap2(_argvec2.size(), _argv2);
  EXPECT_FALSE(ap2.verbose());
}

TEST_F(cargsTest, detectVersion) {
  igp::cargs ap1(_argvec5.size(), _argv5);
  EXPECT_TRUE(ap1.version());
  igp::cargs ap2(_argvec1.size(), _argv1);
  EXPECT_FALSE(ap2.version());
}

TEST_F(cargsTest, basicAccessors) {
  igp::cargs ap(_argvec2.size(), _argv2);
  EXPECT_EQ(ap.get_input_filename(), "fn1");
  EXPECT_EQ(ap.get_input_preset(), "bed");
  EXPECT_EQ(ap.get_recombination_map(), "fn2");
  EXPECT_EQ(ap.get_map_format(), "bolt");
  EXPECT_EQ(ap.get_output_filename(), "fn3");
  EXPECT_EQ(ap.get_output_format(), "map");
}

TEST_F(cargsTest, detectOutputMorgans) {
  igp::cargs ap1(_argvec2.size(), _argv2);
  EXPECT_TRUE(ap1.output_morgans());
  igp::cargs ap2(_argvec3.size(), _argv3);
  EXPECT_FALSE(ap2.output_morgans());
}

TEST_F(cargsTest, detectStepInterval) {
  igp::cargs ap1(_argvec2.size(), _argv2);
  EXPECT_EQ(ap1.get_region_step_interval(), 1.23);
  igp::cargs ap2(_argvec3.size(), _argv3);
  EXPECT_EQ(ap2.get_region_step_interval(), 0.0);
}

TEST_F(cargsTest, detectImproperFormats) {
  igp::cargs ap(_argvec4.size(), _argv4);
  EXPECT_THROW(ap.get_input_preset(), std::runtime_error);
  EXPECT_THROW(ap.get_map_format(), std::runtime_error);
  EXPECT_THROW(ap.get_output_format(), std::runtime_error);
}

TEST_F(cargsTest, detectPrecision) {
  igp::cargs ap(_argvec5.size(), _argv5);
  EXPECT_EQ(ap.get_mpf_precision(), 68u);
}

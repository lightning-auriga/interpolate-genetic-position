/*!
 \file integration_test.h
 \brief mock class for  interfacing with input recombination map iles
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Powered by gmock from
 <https://github.com/google/googletest/tree/main/googlemock>
 */

#ifndef INTEGRATION_TESTS_INTEGRATION_TEST_H_
#define INTEGRATION_TESTS_INTEGRATION_TEST_H_

#include <zlib.h>

#include <iostream>
#include <string>
#include <vector>

#include "boost/filesystem.hpp"
#include "gtest/gtest.h"
#include "interpolate-genetic-position/interpolator.h"

class integrationTest : public testing::Test {
 protected:
  integrationTest();
  ~integrationTest() throw();
  std::string create_plaintext_file(const std::string &filename,
                                    const std::string &content) const;
  std::string create_compressed_file(const std::string &filename,
                                     const std::string &content) const;
  std::string create_bigwig(const std::string &filename,
                            const std::string &content) const;
  std::string load_plaintext_file(const std::string &filename) const;
  std::string load_compressed_file(const std::string &filename) const;
  std::string get_bedfile_content() const;
  std::string get_bim_content() const;
  std::string get_map_content() const;
  std::string get_bolt_content() const;
  std::string get_bedgraph_content() const;
  std::string get_bigwig_content() const;
  const std::string _in_query_tmpfile;
  const std::string _in_gmap_tmpfile;
  const std::string _out_tmpfile;
};
#endif  // INTEGRATION_TESTS_INTEGRATION_TEST_H_

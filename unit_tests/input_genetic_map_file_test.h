/*!
 \file input_genetic_map_file_test.h
 \brief mock class for interfacing with input recombination map iles
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Powered by gmock from
 <https://github.com/google/googletest/tree/main/googlemock>
 */

#ifndef UNIT_TESTS_INPUT_GENETIC_MAP_FILE_TEST_H_
#define UNIT_TESTS_INPUT_GENETIC_MAP_FILE_TEST_H_

#include <string>
#include <vector>

#include "boost/filesystem.hpp"
#include "gmock/gmock.h"
#include "interpolate-genetic-position/input_genetic_map_file.h"

namespace interpolate_genetic_position {
/*!
  \class mock_input_genetic_map_file
  \brief mock interface for testing input query file behavior
 */
class mock_input_genetic_map_file : public base_input_genetic_map_file {
 public:
  // basic constructors/destructors as usual
  mock_input_genetic_map_file();
  ~mock_input_genetic_map_file() throw();
  // mocks
  MOCK_METHOD(void, open, (const std::string &filename, format_type ft),
              (override));
  MOCK_METHOD(bool, get, (), (override));
  MOCK_METHOD(void, close, (), (override));
  MOCK_METHOD(bool, eof, (), (override));
  MOCK_METHOD(std::string, get_chr_lower_bound, (), (const, override));
  MOCK_METHOD(std::string, get_chr_upper_bound, (), (const, override));
  MOCK_METHOD(mpz_class, get_pos_lower_bound, (), (const, override));
  MOCK_METHOD(mpz_class, get_pos_upper_bound, (), (const, override));
  MOCK_METHOD(mpf_class, get_gpos_lower_bound, (), (const, override));
  MOCK_METHOD(mpf_class, get_gpos_upper_bound, (), (const, override));
  MOCK_METHOD(mpf_class, get_rate_lower_bound, (), (const, override));
  MOCK_METHOD(mpf_class, get_rate_upper_bound, (), (const, override));
};
}  // namespace interpolate_genetic_position

class inputGeneticMapFileTest : public testing::Test {
 protected:
  inputGeneticMapFileTest();
  ~inputGeneticMapFileTest() throw();
  const std::string _tmpfile;
  mpf_class _mpf_error_tolerance;
};
#endif  // UNIT_TESTS_INPUT_GENETIC_MAP_FILE_TEST_H_

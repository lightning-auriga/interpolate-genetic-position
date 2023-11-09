/*!
 \file input_variant_file_test.h
 \brief mock class for interfacing with input variant query files
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Powered by gmock from
 <https://github.com/google/googletest/tree/main/googlemock>
 */

#ifndef UNIT_TESTS_INPUT_VARIANT_FILE_TEST_H_
#define UNIT_TESTS_INPUT_VARIANT_FILE_TEST_H_

#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "interpolate-genetic-position/input_variant_file.h"

namespace interpolate_genetic_position {
/*!
  \class mock_input_variant_file
  \brief mock interface for testing input query file behavior
 */
class mock_input_variant_file : public base_input_variant_file {
 public:
  // basic constructors/destructors as usual
  mock_input_variant_file();
  ~mock_input_variant_file() throw();
  // mocks
  MOCK_METHOD(void, open, (const std::string &filename), (override));
  MOCK_METHOD(void, close, (), (override));
  MOCK_METHOD(void, set_format_parameters,
              (unsigned chr_index, unsigned pos1_index, int pos2_index,
               unsigned gpos_index, bool base0, unsigned n_tokens),
              (override));
  MOCK_METHOD(bool, get_input_line, (std::string * line), (override));
  MOCK_METHOD(bool, get_variant, (), (override));
  MOCK_METHOD(const std::string &, get_chr, (), (const, override));
  MOCK_METHOD(const mpz_class &, get_pos1, (), (const, override));
  MOCK_METHOD(const mpz_class &, get_pos2, (), (const, override));
  MOCK_METHOD(const std::vector<std::string> &, get_line_contents, (),
              (const, override));
  MOCK_METHOD(bool, eof, (), (override));
};
}  // namespace interpolate_genetic_position
#endif  // UNIT_TESTS_INPUT_VARIANT_FILE_TEST_H_

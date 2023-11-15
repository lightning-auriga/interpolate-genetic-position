/*!
 \file output_variant_file_test.h
 \brief mock class for interfacing with output variant query files
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Powered by gmock from
 <https://github.com/google/googletest/tree/main/googlemock>
 */

#ifndef UNIT_TESTS_OUTPUT_VARIANT_FILE_TEST_H_
#define UNIT_TESTS_OUTPUT_VARIANT_FILE_TEST_H_

#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "interpolate-genetic-position/output_variant_file.h"

namespace interpolate_genetic_position {
/*!
  \class mock_output_variant_file
  \brief mock interface for testing input query file behavior
 */
class mock_output_variant_file : public base_output_variant_file {
 public:
  // basic constructors/destructors as usual
  mock_output_variant_file();
  ~mock_output_variant_file() throw();
  // mocks
  MOCK_METHOD(void, open, (const std::string &filename, format_type ft),
              (override));
  MOCK_METHOD(void, close, (), (override));
  MOCK_METHOD(void, write,
              (const std::string &chr, const mpz_class &pos1,
               const mpz_class &pos2, const std::string &id,
               const mpf_class &gpos, const mpf_class &rate,
               const std::string &a1, const std::string &a2,
               const double &step_interval),
              (override));
  MOCK_METHOD(format_type, get_format, (), (const, override));
  MOCK_METHOD(void, output_morgans, (bool use_morgans), (override));
  MOCK_METHOD(bool, output_morgans, (), (const, override));
  MOCK_METHOD(std::string, get_last_chr, (), (const, override));
  MOCK_METHOD(void, set_last_chr, (const std::string &), (override));
  MOCK_METHOD(mpf_class, get_last_gpos, (), (const, override));
  MOCK_METHOD(void, set_last_gpos, (const mpf_class &), (override));
  MOCK_METHOD(mpf_class, get_last_rate, (), (const, override));
  MOCK_METHOD(void, set_last_rate, (const mpf_class &), (override));
  MOCK_METHOD(mpz_class, get_last_pos1, (), (const, override));
  MOCK_METHOD(void, set_last_pos1, (const mpz_class &), (override));
  MOCK_METHOD(mpz_class, get_last_pos2, (), (const, override));
  MOCK_METHOD(void, set_last_pos2, (const mpz_class &), (override));
};
}  // namespace interpolate_genetic_position
#endif  // UNIT_TESTS_OUTPUT_VARIANT_FILE_TEST_H_

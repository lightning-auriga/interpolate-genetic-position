/*!
 \file bedfile_queries_test.cc
 \brief tests specific to bedfile query input.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/genetic_map_test.h"
#include "unit_tests/input_genetic_map_file_test.h"

using ::testing::AnyNumber;
using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

TEST_F(geneticMapTest, bedfileQuery1) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(500000));
  EXPECT_CALL(mockfile, get_endpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(1500000));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(1500000));
  EXPECT_CALL(mockfile, get_endpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(2500000));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(1.5)));
  EXPECT_CALL(mockfile, get_rate_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(3.32)));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_startpos = 1000000;
  mpz_class query_endpos = 2000000;
  bool verbose = false;
  gm.query(query_chr, query_startpos, query_endpos, verbose, &result);
  EXPECT_EQ(result.get_chr(), "1");
  EXPECT_EQ(result.get_startpos(), mpz_class(1000000));
  EXPECT_EQ(result.get_endpos(), mpz_class(1500000));
  EXPECT_EQ(result.get_rate(), mpf_class(3.32));
  EXPECT_EQ(cmp(abs(result.get_gpos() - mpf_class(3.16)), _mpf_error_tolerance),
            -1);
}

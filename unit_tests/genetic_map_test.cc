/*!
 \file genetic_map_test.cc
 \brief test of operations of genetic_map.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "unit_tests/genetic_map_test.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/input_genetic_map_file_test.h"

using ::testing::AnyNumber;
using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

geneticMapTest::geneticMapTest()
    : testing::Test(), _mpf_error_tolerance("0.0000000000001") {}

geneticMapTest::~geneticMapTest() {}

TEST_F(geneticMapTest, openCommandExecuted) {
  igp::mock_input_genetic_map_file mockfile;
  std::string filename = "test.ucsc.bedgraph";
  EXPECT_CALL(mockfile, open(filename, igp::BEDGRAPH))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  igp::genetic_map gm(&mockfile);
  gm.open(filename, igp::BEDGRAPH);
}

TEST_F(geneticMapTest, basicConstructorDisabled) {
  EXPECT_THROW(igp::genetic_map gm, std::runtime_error);
}

TEST_F(geneticMapTest, copyConstructorDisabled) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_THROW(igp::genetic_map gm2(gm), std::runtime_error);
}

TEST_F(geneticMapTest, queryBeyondEndOfChromosomePointEstimates) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(900000)));
  EXPECT_CALL(mockfile, get_endpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(-1)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(23.4)));
  EXPECT_CALL(mockfile, get_rate_lower_bound()).Times(1).WillOnce(Return(1.0));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class(23.4));
}

TEST_F(geneticMapTest, queryBeyondEndOfChromosomeRegionEstimates) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(900000)));
  EXPECT_CALL(mockfile, get_endpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(950000)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(23.4)));
  EXPECT_CALL(mockfile, get_rate_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(1.23)));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  mpf_class gpos_expected("23.4615");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
}

TEST_F(geneticMapTest, queryBeyondEndOfFile) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(2).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("23"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("23"));
  igp::query_result result;
  std::string query_chr = "M";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class(0.0));
}

TEST_F(geneticMapTest, queryInteriorUncoveredChromosome) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(1).WillOnce(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("3"));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 1000000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class(0.0));
}

TEST_F(geneticMapTest, queryBeginningOfChromosome) {
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
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 50000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class(0.0));
}

TEST_F(geneticMapTest, queryStandardInterpolation) {
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
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(0.001)));
  EXPECT_CALL(mockfile, get_rate_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class(0.05)));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 150000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  mpf_class gpos_expected("0.0035");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
}

TEST_F(geneticMapTest, queryOverlapLowerBound) {
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
      .WillRepeatedly(Return(mpz_class(100000)));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("0.001")));
  EXPECT_CALL(mockfile, get_rate_lower_bound()).Times(1).WillOnce(Return(1.0));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 100000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class("0.001"));
}

TEST_F(geneticMapTest, queryOverlapUpperBound) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(2).WillRepeatedly(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillOnce(Return(mpz_class(100000)))
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillOnce(Return(mpz_class(200000)))
      .WillRepeatedly(Return(mpz_class(300000)));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("0.002")));
  EXPECT_CALL(mockfile, get_rate_lower_bound()).Times(1).WillOnce(Return(1.0));
  EXPECT_CALL(mockfile, get()).Times(1).WillOnce(Return(true));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 200000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class("0.002"));
}

TEST_F(geneticMapTest, queryDetectsUnsortedChromosomes) {
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("4"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("3"));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 100000;
  bool verbose = false;
  EXPECT_THROW(gm.query(query_chr, query_pos, -1, verbose, &result),
               std::runtime_error);
}

TEST_F(geneticMapTest, queryDetectsUnsortedPositions) {
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
      .WillRepeatedly(Return(mpz_class(200000)));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpz_class(100000)));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 5000;
  bool verbose = false;
  EXPECT_THROW(gm.query(query_chr, query_pos, -1, verbose, &result),
               std::runtime_error);
}

TEST_F(geneticMapTest, queryIncrementGeneticMapWhenExceedsPositionRange) {
  // if the query and genetic map are on the same chromosome but the query
  // position exceeds the upper chromosome, the map should be incremented
  // exactly once.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(2).WillRepeatedly(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillOnce(Return(1000))
      .WillRepeatedly(Return(2000));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillOnce(Return(2000))
      .WillRepeatedly(Return(3000));
  EXPECT_CALL(mockfile, get_gpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("0.245")));
  EXPECT_CALL(mockfile, get_rate_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("1.23")));
  EXPECT_CALL(mockfile, get()).Times(1).WillOnce(Return(true));
  igp::query_result result;
  std::string query_chr = "1";
  mpz_class query_pos = 2500;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  mpf_class gpos_expected("0.245615");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
}

TEST_F(geneticMapTest, queryIncrementGeneticMapWhenExceedsLowerChromosome) {
  // if the query exceeds the lower chromosome but not the upper,
  // the range should be incremented exactly once.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(2).WillRepeatedly(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillOnce(Return("1"))
      .WillOnce(Return("1"))
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(100000));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get()).Times(AnyNumber()).WillRepeatedly(Return(true));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 50000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class("0.0"));
}

TEST_F(geneticMapTest, queryIncrementGeneticMapWhenExceedsUpperChromosome) {
  // if the query exceeds the upper chromosome,
  // the range should be incremented exactly twice.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(2).WillRepeatedly(Return(false));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillOnce(Return("1"))
      .WillOnce(Return("1"))
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillOnce(Return("1"))
      .WillOnce(Return("1"))
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(100000));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get()).Times(AnyNumber()).WillRepeatedly(Return(true));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 50000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_gpos(), mpf_class("0.0"));
}

TEST_F(geneticMapTest, queryChromMatchEndOfBoltStyleMap) {
  // if the query is exactly the same as the final entry
  // in the genetic map, it should return as usual
  // when there's a direct collision; but in the implementation,
  // the logic is separate and needs separate testing.
  // In the instance that the genetic map is bolt-style, without
  // an end position per estimate, the rate is guaranteed to be 0
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(AnyNumber()).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get_endpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(-1));
  EXPECT_CALL(mockfile, get_gpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("234.5")));
  EXPECT_CALL(mockfile, get_rate_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("0.0")));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 200000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_startpos(), mpz_class(200000));
  EXPECT_EQ(result.get_endpos(), mpz_class(-1));
  EXPECT_EQ(result.get_gpos(), mpf_class("234.5"));
  EXPECT_EQ(result.get_rate(), mpf_class("0.0"));
}

TEST_F(geneticMapTest, querySameChromBeyondEndOfBedgraphStyleMap) {
  // if the query is exactly the same as the final entry
  // in the genetic map, it should return as usual
  // when there's a direct collision; but in the implementation,
  // the logic is separate and needs separate testing.
  // In the instance that the genetic map is bedgraph/bigwig-style,
  // with an end position per estimate, the rate is probably nonzero
  // and some amount of interpolation will be required.
  // For queries beyond that end position, the effective rate
  // must be set to 0.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(AnyNumber()).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get_endpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(300000));
  EXPECT_CALL(mockfile, get_gpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("234.5")));
  EXPECT_CALL(mockfile, get_rate_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("20.1")));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 400000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_startpos(), mpz_class(400000));
  EXPECT_EQ(result.get_endpos(), mpz_class(-1));
  mpf_class gpos_expected("236.51");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
  EXPECT_EQ(result.get_rate(), mpf_class("0.0"));
}

TEST_F(geneticMapTest, querySameChromWithinEndOfBedgraphStyleMap) {
  // if the query is exactly the same as the final entry
  // in the genetic map, it should return as usual
  // when there's a direct collision; but in the implementation,
  // the logic is separate and needs separate testing.
  // In the instance that the genetic map is bedgraph/bigwig-style,
  // with an end position per estimate, the rate is probably nonzero
  // and some amount of interpolation will be required.
  // For queries within the final range, the effective rate
  // is whatever is set for the final range.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(AnyNumber()).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get_endpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(300000));
  EXPECT_CALL(mockfile, get_gpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("234.5")));
  EXPECT_CALL(mockfile, get_rate_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("20.1")));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_pos = 250000;
  bool verbose = false;
  gm.query(query_chr, query_pos, -1, verbose, &result);
  EXPECT_EQ(result.get_startpos(), mpz_class(250000));
  EXPECT_EQ(result.get_endpos(), mpz_class(-1));
  mpf_class gpos_expected("235.505");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
  EXPECT_EQ(result.get_rate(), mpf_class("20.1"));
}

TEST_F(geneticMapTest,
       querySameChromWithinEndOfBedgraphStyleMapContainedRegion) {
  // if the query is exactly the same as the final entry
  // in the genetic map, it should return as usual
  // when there's a direct collision; but in the implementation,
  // the logic is separate and needs separate testing.
  // In the instance that the genetic map is bedgraph/bigwig-style,
  // with an end position per estimate, the rate is probably nonzero
  // and some amount of interpolation will be required.
  // For queries within the final range, the effective rate
  // is whatever is set for the final range.
  // Query is region contained within final map interval.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(AnyNumber()).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get_endpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(300000));
  EXPECT_CALL(mockfile, get_gpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("234.5")));
  EXPECT_CALL(mockfile, get_rate_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("20.1")));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_startpos = 250000;
  mpz_class query_endpos = 275000;
  bool verbose = false;
  gm.query(query_chr, query_startpos, query_endpos, verbose, &result);
  EXPECT_EQ(result.get_startpos(), mpz_class(250000));
  EXPECT_EQ(result.get_endpos(), mpz_class(275000));
  mpf_class gpos_expected("235.505");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
  EXPECT_EQ(result.get_rate(), mpf_class("20.1"));
}

TEST_F(geneticMapTest,
       querySameChromOverlappingEndOfBedgraphStyleMapContainedRegion) {
  // if the query is exactly the same as the final entry
  // in the genetic map, it should return as usual
  // when there's a direct collision; but in the implementation,
  // the logic is separate and needs separate testing.
  // In the instance that the genetic map is bedgraph/bigwig-style,
  // with an end position per estimate, the rate is probably nonzero
  // and some amount of interpolation will be required.
  // For queries within the final range, the effective rate
  // is whatever is set for the final range.
  // Query is region overlapping endpoint of final map interval.
  igp::mock_input_genetic_map_file mockfile;
  igp::genetic_map gm(&mockfile);
  EXPECT_CALL(mockfile, close()).Times(AnyNumber());
  EXPECT_CALL(mockfile, eof()).Times(AnyNumber()).WillRepeatedly(Return(true));
  EXPECT_CALL(mockfile, get_chr_lower_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("1"));
  EXPECT_CALL(mockfile, get_chr_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return("2"));
  EXPECT_CALL(mockfile, get_startpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(200000));
  EXPECT_CALL(mockfile, get_endpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(300000));
  EXPECT_CALL(mockfile, get_gpos_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("234.5")));
  EXPECT_CALL(mockfile, get_rate_upper_bound())
      .Times(AnyNumber())
      .WillRepeatedly(Return(mpf_class("20.1")));
  igp::query_result result;
  std::string query_chr = "2";
  mpz_class query_startpos = 250000;
  mpz_class query_endpos = 350000;
  bool verbose = false;
  gm.query(query_chr, query_startpos, query_endpos, verbose, &result);
  EXPECT_EQ(result.get_startpos(), mpz_class(250000));
  EXPECT_EQ(result.get_endpos(), mpz_class(300000));
  mpf_class gpos_expected("235.505");
  EXPECT_EQ(cmp(abs(result.get_gpos() - gpos_expected), _mpf_error_tolerance),
            -1);
  EXPECT_EQ(result.get_rate(), mpf_class("20.1"));
}

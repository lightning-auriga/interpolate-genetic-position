/*!
 \file query_file_test.cc
 \brief test of variant query file operations.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/query_file.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "unit_tests/input_variant_file_test.h"
#include "unit_tests/output_variant_file_test.h"

using ::testing::AtLeast;
using ::testing::Return;

namespace igp = interpolate_genetic_position;

TEST(queryFileTest, canInitializeBimfile) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string bimfilename = "test.bim";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 3, -1, 2, false, 6))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(bimfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(bimfilename, igp::BIM);
}

TEST(queryFileTest, canInitializeMapfile) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string mapfilename = "test.map";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 3, -1, 2, false, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(mapfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(mapfilename, igp::MAP);
}

TEST(queryFileTest, canInitializeBedfile) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string bedfilename = "test.bed";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 1, 2, 4, true, 4))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(bedfilename)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(bedfilename, igp::BED);
}

TEST(queryFileTest, canInitializeVcf) {
  igp::mock_input_variant_file mock_infile;
  igp::mock_output_variant_file mock_outfile;
  std::string vcfname = "test.vcf";
  EXPECT_CALL(mock_infile, set_format_parameters(0, 1, -1, -1, false, 9))
      .Times(1)
      .WillOnce(Return());
  EXPECT_CALL(mock_infile, open(vcfname)).Times(1).WillOnce(Return());
  igp::query_file qf(&mock_infile, &mock_outfile);
  qf.open(vcfname, igp::VCF);
}

TEST(queryFileTest, canDetectInvalidVcf) {
  igp::input_variant_file infile;
  igp::output_variant_file outfile;
  std::string fakename = "test.vcf";
  igp::query_file qf(&infile, &outfile);
  EXPECT_THROW(qf.open(fakename, igp::VCF), std::runtime_error);
}

TEST(queryFileTest, canReadFromCin) {
  std::string map_content =
      "1 rs1 0 100\n"
      "1 rs2 0 200\n"
      "2 rs3 0 300\n";
  std::istringstream strm1(map_content);
  igp::input_variant_file infile;
  infile.set_fallback_stream(&strm1);
  igp::mock_output_variant_file mock_outfile;
  igp::query_file qf(&infile, &mock_outfile);
  qf.open("", igp::MAP);
  EXPECT_TRUE(qf.get());
  EXPECT_EQ(qf.get_chr(), "1");
  EXPECT_EQ(qf.get_pos1(), mpz_class(100));
  EXPECT_EQ(qf.get_pos2(), mpz_class(-1));
  EXPECT_TRUE(qf.get());
  EXPECT_EQ(qf.get_chr(), "1");
  EXPECT_EQ(qf.get_pos1(), mpz_class(200));
  EXPECT_EQ(qf.get_pos2(), mpz_class(-1));
  EXPECT_TRUE(qf.get());
  EXPECT_EQ(qf.get_chr(), "2");
  EXPECT_EQ(qf.get_pos1(), mpz_class(300));
  EXPECT_EQ(qf.get_pos2(), mpz_class(-1));
  EXPECT_FALSE(qf.get());
}

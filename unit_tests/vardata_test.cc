/*!
 \file vardata_test.cc
 \brief test of variant metadata storage class.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "gtest/gtest.h"
#include "unit_tests/input_variant_file_test.h"

namespace igp = interpolate_genetic_position;

TEST(vardataTest, defaultConstructor) { igp::vardata vd; }

TEST(vardataTest, copyConstructor) {
  igp::vardata vd1;
  vd1.set_chr("chr1");
  vd1.set_pos1(100);
  vd1.set_pos2(200);
  vd1.set_a1("A");
  vd1.set_a2("T");
  vd1.set_varid("rs1");
  vd1.set_annotation("annot");
  igp::vardata vd2(vd1);
  EXPECT_EQ(vd1.get_chr(), vd2.get_chr());
  EXPECT_EQ(vd1.get_pos1(), vd2.get_pos1());
  EXPECT_EQ(vd1.get_pos2(), vd2.get_pos2());
  EXPECT_EQ(vd1.get_a1(), vd2.get_a1());
  EXPECT_EQ(vd1.get_a2(), vd2.get_a2());
  EXPECT_EQ(vd1.get_varid(), vd2.get_varid());
  EXPECT_EQ(vd1.get_annotation(), vd2.get_annotation());
}

TEST(vardataTest, assignmentOperator) {
  igp::vardata vd1;
  vd1.set_chr("chr1");
  vd1.set_pos1(100);
  vd1.set_pos2(200);
  vd1.set_a1("A");
  vd1.set_a2("T");
  vd1.set_varid("rs1");
  vd1.set_annotation("annot");
  igp::vardata vd2 = vd1;
  EXPECT_EQ(vd1.get_chr(), vd2.get_chr());
  EXPECT_EQ(vd1.get_pos1(), vd2.get_pos1());
  EXPECT_EQ(vd1.get_pos2(), vd2.get_pos2());
  EXPECT_EQ(vd1.get_a1(), vd2.get_a1());
  EXPECT_EQ(vd1.get_a2(), vd2.get_a2());
  EXPECT_EQ(vd1.get_varid(), vd2.get_varid());
  EXPECT_EQ(vd1.get_annotation(), vd2.get_annotation());
}

TEST(vardataTest, testGetSet) {
  igp::vardata vd;
  vd.set_chr("chr1");
  vd.set_pos1(100);
  vd.set_pos2(200);
  vd.set_a1("A");
  vd.set_a2("T");
  vd.set_varid("rs1");
  vd.set_annotation("annot");
  EXPECT_EQ(vd.get_chr(), "chr1");
  EXPECT_EQ(vd.get_pos1(), mpz_class(100));
  EXPECT_EQ(vd.get_pos2(), mpz_class(200));
  EXPECT_EQ(vd.get_a1(), "A");
  EXPECT_EQ(vd.get_a2(), "T");
  EXPECT_EQ(vd.get_varid(), "rs1");
  EXPECT_EQ(vd.get_annotation(), "annot");
}

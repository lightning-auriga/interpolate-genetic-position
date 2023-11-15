/*!
 \file integration_test.cc
 \brief test of full interpolation passes
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "integration_tests/integration_test.h"

#include "gtest/gtest.h"

namespace igp = interpolate_genetic_position;

integrationTest::integrationTest()
    : testing::Test(),
      _in_query_tmpfile(boost::filesystem::unique_path().native()),
      _in_gmap_tmpfile(boost::filesystem::unique_path().native()),
      _out_tmpfile(boost::filesystem::unique_path().native()) {}

integrationTest::~integrationTest() {
  if (boost::filesystem::exists(_in_query_tmpfile)) {
    boost::filesystem::remove(_in_query_tmpfile);
  }
  if (boost::filesystem::exists(_in_gmap_tmpfile)) {
    boost::filesystem::remove(_in_gmap_tmpfile);
  }
  if (boost::filesystem::exists(_out_tmpfile)) {
    boost::filesystem::remove(_out_tmpfile);
  }
}

std::string integrationTest::create_plaintext_file(
    const std::string &filename, const std::string &content) const {
  std::ofstream output;
  output.open(filename.c_str());
  if (!output.is_open()) {
    throw std::runtime_error(
        "create_plaintext_file: cannot open write connection");
  }
  if (!(output << content)) {
    throw std::runtime_error("create_plaintext_file: cannot write to file");
  }
  output.close();
  return content;
}

std::string integrationTest::create_compressed_file(
    const std::string &filename, const std::string &content) const {
  gzFile output = NULL;
  try {
    output = gzopen(filename.c_str(), "wb");
    if (!output) {
      throw std::runtime_error(
          "create_compressed_file: cannot open write connection");
    }
    if (gzwrite(output, content.c_str(), content.size()) !=
        static_cast<int>(content.size())) {
      throw std::runtime_error("create_compressed_file: cannot write to file");
    }
    gzclose(output);
    output = NULL;
  } catch (...) {
    if (output) {
      gzclose(output);
      output = NULL;
    }
    throw;
  }
  return content;
}

std::string integrationTest::create_bigwig(const std::string &filename,
                                           const std::string &content) const {
  // requires direct interaction with libBigWig
  return "";
}

std::string integrationTest::load_plaintext_file(
    const std::string &filename) const {
  std::ifstream input;
  std::string res = "", line = "";
  input.open(filename.c_str());
  if (!input.is_open()) {
    throw std::runtime_error("load_plaintext_file: unable to open file");
  }
  while (input.peek() != EOF) {
    getline(input, line);
    res += line + "\n";
  }
  input.close();
  return res;
}

std::string integrationTest::load_compressed_file(
    const std::string &filename) const {
  gzFile input = NULL;
  std::string res = "";
  char *buffer = NULL;
  int buffer_size = 2000;
  try {
    buffer = new char[buffer_size];
    input = gzopen(filename.c_str(), "rb");
    if (!input) {
      throw std::runtime_error("load_compressed_file: unable to open file");
    }
    while (gzgets(input, buffer, buffer_size) != Z_NULL) {
      res += std::string(buffer);
    }
    gzclose(input);
    input = 0;
    delete[] buffer;
    buffer = 0;
    return res;
  } catch (...) {
    if (input) {
      gzclose(input);
    }
    if (buffer) {
      delete[] buffer;
    }
    throw;
  }
}

std::string integrationTest::get_bedfile_content() const {
  return "chr1 499999 1499999 0\n"
         "chr1 1499999 2499999 0\n"
         "chr3 999999 1999999 0\n";
}

std::string integrationTest::get_bim_content() const { return ""; }

std::string integrationTest::get_map_content() const { return ""; }

std::string integrationTest::get_bolt_content() const {
  return "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
         "1 1000000 0.1 0\n"
         "1 2000000 0.2 0.1\n"
         "1 3000000 0.25 0.3\n"
         "2 1000000 0.1 0\n"
         "2 3000000 0.25 0.3\n"
         "2 5000000 0.55 0.8\n";
}

std::string integrationTest::get_bedgraph_content() const { return ""; }

std::string integrationTest::get_bigwig_content() const { return ""; }

TEST_F(integrationTest, bedfileInputBoltOutput) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_bedfile_content());
  std::string input_gmap =
      create_plaintext_file(_in_gmap_tmpfile, get_bolt_content());
  std::string expected_output =
      "chr1\t499999\t999999\t0\n"
      "chr1\t999999\t1499999\t0.1\n"
      "chr1\t1499999\t1999999\t0.1\n"
      "chr1\t1999999\t2499999\t0.2\n"
      "chr3\t999999\t1999999\t0\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "bed", _in_gmap_tmpfile, "bolt",
                 _out_tmpfile, false, 0.0, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

TEST_F(integrationTest, bimfileInputBimfileOutput) {}

TEST_F(integrationTest, mapfileInputMapfileOutput) {}

TEST_F(integrationTest, bedgraphGeneticMap) {}

TEST_F(integrationTest, bigwigGeneticMap) {}

TEST_F(integrationTest, boltGeneticMap) {}

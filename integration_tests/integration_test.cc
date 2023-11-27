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

void integrationTest::write_bigwig_content(const std::string &filename) const {
  bigWigFile_t *fp = NULL;
  const char *chroms[] = {"chr1", "chr2"};
  const char *chromsUse[] = {"chr1", "chr1", "chr1", "chr2", "chr2", "chr2"};
  uint32_t chrLens[] = {248956422, 242193529};
  uint32_t starts[] = {999999, 1999999, 2999999, 999999, 2999999, 4999999};
  uint32_t ends[] = {1999999, 2999999, 3999999, 2999999, 4999999, 5999999};
  float values[] = {0.1f, 0.2f, 0.25f, 0.1f, 0.25f, 0.55f};
  try {
    if (bwInit(1 << 17) != 0) {
      throw std::runtime_error("write_bigwig_content: bwInit failed");
    }
    fp = bwOpen(filename.c_str(), NULL, "w");
    if (!fp) {
      throw std::runtime_error("write_bigwig_content: bwOpen failed");
    }
    if (bwCreateHdr(fp, 10)) {
      throw std::runtime_error("write_bigwig_content: bwCreateHdr failed");
    }
    fp->cl = bwCreateChromList(chroms, chrLens, 2);
    if (!fp->cl) {
      throw std::runtime_error(
          "write_bigwig_content: bwCreateChromList failed");
    }
    if (bwWriteHdr(fp)) {
      throw std::runtime_error("write_bigwig_content: bwWriteHdr failed");
    }
    if (bwAddIntervals(fp, chromsUse, starts, ends, values, 3)) {
      throw std::runtime_error(
          "write_bigwig_content: bwAddIntervals (1) failed");
    }
    if (bwAddIntervals(fp, chromsUse + 3, starts + 3, ends + 3, values + 3,
                       3)) {
      throw std::runtime_error(
          "write_bigwig_content: bwAddIntervals (2) failed");
    }
    bwClose(fp);
    bwCleanup();
  } catch (...) {
    if (fp) {
      bwClose(fp);
    }
    bwCleanup();
    throw;
  }
}

std::string integrationTest::get_bedfile_content() const {
  return "chr1 499999 1199999 0\n"
         "chr1 1199999 1299999 0\n"
         "chr1 1299999 1499999 0\n"
         "chr1 1499999 2499999 0\n"
         "chr3 999999 1999999 0\n";
}

std::string integrationTest::get_bim_content() const {
  return "1 rs1 0 500000 A T\n"
         "1 rs2 0 1500000 C G\n"
         "3 rs3 0 1000000 A C\n";
}

std::string integrationTest::get_map_content() const {
  return "1 rs1 0 500000\n"
         "1 rs2 0 1500000\n"
         "3 rs3 0 1000000\n";
}

std::string integrationTest::get_bolt_content() const {
  return "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
         "1 1000000 0.1 0\n"
         "1 2000000 0.2 0.1\n"
         "1 3000000 0.25 0.3\n"
         "2 1000000 0.1 0\n"
         "2 3000000 0.25 0.3\n"
         "2 5000000 0.55 0.8\n";
}

std::string integrationTest::get_bedgraph_content() const {
  return "chr1 999999 1999999 0.1\n"
         "chr1 1999999 2999999 0.2\n"
         "chr1 2999999 3999999 0.25\n"
         "chr2 999999 2999999 0.1\n"
         "chr2 2999999 4999999 0.25\n"
         "chr2 4999999 5999999 0.55\n";
}

TEST_F(integrationTest, bedfileInputBoltOutputNoIncrement) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_bedfile_content());
  std::string input_gmap =
      create_plaintext_file(_in_gmap_tmpfile, get_bolt_content());
  std::string expected_output =
      "chr1\t499999\t0\t0\n"
      "chr1\t999999\t0.1\t0\n"
      "chr1\t1199999\t0.1\t0.02\n"
      "chr1\t1299999\t0.1\t0.03\n"
      "chr1\t1499999\t0.1\t0.05\n"
      "chr1\t1999999\t0.2\t0.1\n"
      "chr1\t2499999\t0\t0.2\n"
      "chr3\t999999\t0\t0\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "bed", _in_gmap_tmpfile, "bolt",
                 _out_tmpfile, "bed", false, 0.0, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

TEST_F(integrationTest, bedfileInputBedgraphOutputBoundaryIncrement) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_bedfile_content());
  std::string input_gmap =
      create_plaintext_file(_in_gmap_tmpfile, get_bolt_content());
  double boundary_increment = 1.0;
  std::string expected_output =
      "chr1\t499999\t0\t0\n"
      "chr1\t999999\t0.1\t0\n"
      "chr1\t1199999\t0.1\t1.02\n"
      "chr1\t1299999\t0.1\t2.03\n"
      "chr1\t1499999\t0.1\t3.05\n"
      "chr1\t1999999\t0.2\t3.1\n"
      "chr1\t2499999\t0\t4.2\n"
      "chr3\t999999\t0\t0\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "bed", _in_gmap_tmpfile, "bolt",
                 _out_tmpfile, "bed", false, boundary_increment, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

TEST_F(integrationTest, bimfileInputBimfileOutput) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_bim_content());
  std::string input_gmap =
      create_plaintext_file(_in_gmap_tmpfile, get_bolt_content());
  std::string expected_output =
      "1\trs1\t0\t500000\tA\tT\n"
      "1\trs2\t0.05\t1500000\tC\tG\n"
      "3\trs3\t0\t1000000\tA\tC\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "bim", _in_gmap_tmpfile, "bolt",
                 _out_tmpfile, "bim", false, 0.0, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

TEST_F(integrationTest, mapfileInputMapfileOutput) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_map_content());
  std::string input_gmap =
      create_plaintext_file(_in_gmap_tmpfile, get_bolt_content());
  std::string expected_output =
      "1\trs1\t0\t500000\n"
      "1\trs2\t0.05\t1500000\n"
      "3\trs3\t0\t1000000\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "map", _in_gmap_tmpfile, "bolt",
                 _out_tmpfile, "map", false, 0.0, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

TEST_F(integrationTest, gzippedMapfileInputMapfileOutput) {
  const std::string in_query_tmpfile =
      boost::filesystem::unique_path().native() + ".gz";
  try {
    std::string input_query =
        create_compressed_file(in_query_tmpfile, get_map_content());
    std::string input_gmap =
        create_plaintext_file(_in_gmap_tmpfile, get_bolt_content());
    std::string expected_output =
        "1\trs1\t0\t500000\n"
        "1\trs2\t0.05\t1500000\n"
        "3\trs3\t0\t1000000\n";
    igp::interpolator ip;
    ip.interpolate(in_query_tmpfile, "map", _in_gmap_tmpfile, "bolt",
                   _out_tmpfile, "map", false, 0.0, false);
    EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
    std::string observed_output = load_plaintext_file(_out_tmpfile);
    EXPECT_EQ(expected_output, observed_output);
    boost::filesystem::remove(in_query_tmpfile);
  } catch (...) {
    if (boost::filesystem::exists(in_query_tmpfile)) {
      boost::filesystem::remove(in_query_tmpfile);
    }
    throw;
  }
}

TEST_F(integrationTest, gzippedBoltGeneticMap) {
  const std::string in_gmap_tmpfile =
      boost::filesystem::unique_path().native() + ".gz";
  try {
    std::string input_query =
        create_plaintext_file(_in_query_tmpfile, get_map_content());
    std::string input_gmap =
        create_compressed_file(in_gmap_tmpfile, get_bolt_content());
    std::string expected_output =
        "1\trs1\t0\t500000\n"
        "1\trs2\t0.05\t1500000\n"
        "3\trs3\t0\t1000000\n";
    igp::interpolator ip;
    ip.interpolate(_in_query_tmpfile, "map", in_gmap_tmpfile, "bolt",
                   _out_tmpfile, "map", false, 0.0, false);
    EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
    std::string observed_output = load_plaintext_file(_out_tmpfile);
    EXPECT_EQ(expected_output, observed_output);
    boost::filesystem::remove(in_gmap_tmpfile);
  } catch (...) {
    if (boost::filesystem::exists(in_gmap_tmpfile)) {
      boost::filesystem::remove(in_gmap_tmpfile);
    }
  }
}

TEST_F(integrationTest, bedgraphGeneticMap) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_map_content());
  std::string input_gmap =
      create_plaintext_file(_in_gmap_tmpfile, get_bedgraph_content());
  std::string expected_output =
      "1\trs1\t0\t500000\n"
      "1\trs2\t0.05\t1500000\n"
      "3\trs3\t0\t1000000\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "map", _in_gmap_tmpfile, "bedgraph",
                 _out_tmpfile, "map", false, 0.0, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

TEST_F(integrationTest, bigwigGeneticMap) {
  std::string input_query =
      create_plaintext_file(_in_query_tmpfile, get_map_content());
  write_bigwig_content(_in_gmap_tmpfile);
  std::string expected_output =
      "1\trs1\t0\t500000\n"
      "1\trs2\t0.05\t1500000\n"
      "3\trs3\t0\t1000000\n";
  igp::interpolator ip;
  ip.interpolate(_in_query_tmpfile, "map", _in_gmap_tmpfile, "bigwig",
                 _out_tmpfile, "map", false, 0.0, false);
  EXPECT_TRUE(boost::filesystem::exists(_out_tmpfile));
  std::string observed_output = load_plaintext_file(_out_tmpfile);
  EXPECT_EQ(expected_output, observed_output);
}

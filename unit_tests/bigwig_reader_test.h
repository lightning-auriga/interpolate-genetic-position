/*!
 \file bigwig_reader_test.h
 \brief declaration of test suite for bigwig reader class
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#ifndef UNIT_TESTS_BIGWIG_READER_TEST_H_
#define UNIT_TESTS_BIGWIG_READER_TEST_H_

#include <cfloat>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "interpolate-genetic-position/bigwig_reader.h"

class bigwigReaderTest : public testing::Test {
 protected:
  bigwigReaderTest();
  ~bigwigReaderTest() throw();
  mpf_class _mpf_error_tolerance;
};
#endif  // UNIT_TESTS_BIGWIG_READER_TEST_H_

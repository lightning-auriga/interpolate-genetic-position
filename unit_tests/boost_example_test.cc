/*!
 \file boost_example_test.cc
 \brief example boost unit test outside of main source file.
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga

 Docs at <https://www.boost.org/doc/libs/1_83_0/libs/test/doc/html/index.html>
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(example_test_2) { BOOST_TEST(2 == 2); }

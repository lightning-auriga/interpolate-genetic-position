/*!
 \file boost_main.cc
 \brief example boost unit test module declaration and single test automatic
 registration. \copyright Released under the MIT License. Copyright 2023
 Lightning Auriga

 Docs at <https://www.boost.org/doc/libs/1_83_0/libs/test/doc/html/index.html>
 */

#define BOOST_TEST_MODULE interpolate - genetic - position unit test module
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(example_test_1) { BOOST_TEST(1 == 1); }

bin_PROGRAMS = interpolate-genetic-position.out test_suite.out

AM_CXXFLAGS = $(BOOST_CPPFLAGS) -ggdb -Wall -std=c++17

COMBINED_SOURCES = interpolate-genetic-position/cargs.cc interpolate-genetic-position/cargs.h interpolate-genetic-position/genetic_map.cc interpolate-genetic-position/genetic_map.h interpolate-genetic-position/interpolator.cc interpolate-genetic-position/interpolator.h interpolate-genetic-position/query_file.cc interpolate-genetic-position/query_file.h interpolate-genetic-position/utilities.cc interpolate-genetic-position/utilities.h
COMBINED_LDADD = $(BOOST_LDFLAGS) -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgmpxx -lgmp -lz

interpolate_genetic_position_out_SOURCES = interpolate-genetic-position/main.cc $(COMBINED_SOURCES)
interpolate_genetic_position_out_LDADD = $(COMBINED_LDADD)
UNIT_TEST_SOURCES = unit_tests/google_example_test.cc
test_suite_out_SOURCES = $(COMBINED_SOURCES) $(UNIT_TEST_SOURCES)
test_suite_out_LDADD = $(COMBINED_LDADD) -lgtest_main -lgtest -lgmock -lpthread

dist_doc_DATA = README
ACLOCAL_AMFLAGS = -I m4

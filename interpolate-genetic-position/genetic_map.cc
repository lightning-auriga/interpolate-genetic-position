/*!
 \file genetic_map.cc
 \brief implementation for reference genetic map file
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/genetic_map.h"

namespace igp = interpolate_genetic_position;

igp::genetic_map::genetic_map() {
  throw std::runtime_error("genetic_map: basic constructor disabled");
}
igp::genetic_map::genetic_map(const genetic_map &obj) {
  throw std::runtime_error(
      "genetic_map: copy constructor operation is invalid for this class");
}
igp::genetic_map::genetic_map(base_input_genetic_map_file *ptr)
    : _interface(ptr) {}
igp::genetic_map::~genetic_map() throw() { _interface->close(); }
void igp::genetic_map::open(const std::string &filename, format_type ft) {
  _interface->open(filename, ft);
}
void igp::genetic_map::query(const std::string &chr_query,
                             const mpz_class &pos1_query,
                             const mpz_class &pos2_query, bool verbose,
                             std::vector<query_result> *results) {
  results->clear();
  query_result result;
  // this was all very straightforward until I considered how to support
  // queries based on bedfiles. regions that cross boundaries in the
  // genetic map file need to return multiple results, each with their
  // own position and rate estimates. what a mess.
  mpz_class current_pos1 = pos1_query;
  while (true) {
    query(chr_query, current_pos1, pos2_query, verbose, &result);
    results->push_back(result);
    // if the initial query was a point estimate, or if the result reaches
    // the end of the query range, it's done.
    if (cmp(pos2_query, 0) == -1 || cmp(result.get_endpos(), pos2_query) == 0) {
      break;
    }
    // update the query to reflect the segment that was not interpolated
    // in this loop
    current_pos1 = result.get_endpos();
  }
}
void igp::genetic_map::query(const std::string &chr_query,
                             const mpz_class &pos1_query,
                             const mpz_class &pos2_query, bool verbose,
                             query_result *result) {
  // This comparison must cautiously detect whether a variant falls before,
  // after, or within the interval indicated by the currently loaded set of
  // genetic positions. Certain permutations of this comparison will
  // indicate that sorting assumptions on input data have been violated.
  // Additionally, the test for chromosome equality must be performed
  // carefully to deal with the many permutations of chromosome codes
  // that are encountered in datasets.
  direction query_vs_lower_bound = EQUAL;
  direction query_vs_upper_bound = EQUAL;
  mpf_class mb_adjustment = 1000000.0;
  if (verbose) {
    std::cout << "query: chr is " << chr_query << ", pos is " << pos1_query
              << std::endl;
  }
  result->set_chr(chr_query);
  result->set_startpos(pos1_query);
  std::string chr_lower_bound = _interface->get_chr_lower_bound();
  std::string chr_upper_bound = _interface->get_chr_upper_bound();
  if (chromosome_compare(chr_lower_bound, chr_upper_bound) == GREATER_THAN) {
    throw std::runtime_error(
        "genetic_map::query: your genetic map "
        "is unsorted. For now, the only solution to this issue is "
        "to sort your input data; in the future, there may be an "
        "option to have the program handle this for you, at the cost "
        "of RAM.");
  }
  while (!_interface->eof()) {
    chr_lower_bound = _interface->get_chr_lower_bound();
    chr_upper_bound = _interface->get_chr_upper_bound();
    query_vs_lower_bound = chromosome_compare(chr_query, chr_lower_bound);
    query_vs_upper_bound = chromosome_compare(chr_query, chr_upper_bound);
    if (query_vs_lower_bound == EQUAL && query_vs_upper_bound == EQUAL) {
      if (verbose) {
        std::cout << "\tsame chromosome, position comparing" << std::endl;
      }
      // proceed to position comparison
      // this early check enforces strict position ordering, and prevents
      // several later error conditions from ever occurring.
      mpz_class startpos_lower_bound = _interface->get_startpos_lower_bound();
      mpz_class startpos_upper_bound = _interface->get_startpos_upper_bound();
      if (cmp(startpos_lower_bound, startpos_upper_bound) != -1) {
        throw std::runtime_error(
            "genetic_map::query: your genetic map "
            "is unsorted. For now, the only solution to this issue is "
            "to sort your input data; in the future, there may be an "
            "option to have the program handle this for you, at the cost "
            "of RAM.");
      }
      int query_startpos_vs_lower_bound = cmp(pos1_query, startpos_lower_bound);
      int query_startpos_vs_upper_bound = cmp(pos1_query, startpos_upper_bound);
      if (query_startpos_vs_lower_bound == 0 &&
          query_startpos_vs_upper_bound == -1) {
        if (verbose) {
          std::cout << "\t\tmatches lower boundary exactly: "
                    << _interface->get_gpos_lower_bound() << std::endl;
        }
        // exact lower boundary rate
        result->set_gpos(_interface->get_gpos_lower_bound());
        result->set_endpos(cmp(pos2_query, startpos_upper_bound) < 0
                               ? pos2_query
                               : startpos_upper_bound);
        result->set_rate(_interface->get_rate_lower_bound());
        return;
      } else if (query_startpos_vs_lower_bound == 1 &&
                 query_startpos_vs_upper_bound == 0) {
        // exact upper boundary rate
        // since this is *waves hands* the same thing as the exact lower
        // boundary rate, get() once here and have the other logic deal with it
        if (verbose) {
          std::cout << "\t\tmatches upper boundary exactly, so incrementing to "
                       "use that logic"
                    << std::endl;
        }
        _interface->get();
      } else if (query_startpos_vs_lower_bound == -1 &&
                 query_startpos_vs_upper_bound == -1) {
        // beginning of chromosome
        if (verbose) {
          std::cout << "\t\tbeginning of chromosome, before rate estimates "
                       "start; setting to 0"
                    << std::endl;
        }
        result->set_gpos(0.0);
        result->set_endpos(cmp(pos2_query, startpos_lower_bound) < 0
                               ? pos2_query
                               : startpos_lower_bound);
        result->set_rate(0.0);
        return;
      } else if (query_startpos_vs_lower_bound == 1 &&
                 query_startpos_vs_upper_bound == -1) {
        // interpolate
        result->set_gpos(_interface->get_gpos_lower_bound() +
                         (pos1_query - startpos_lower_bound) / mb_adjustment *
                             _interface->get_rate_lower_bound());
        if (verbose) {
          std::cout << "interpolated: gpos_lower_bound = "
                    << _interface->get_gpos_lower_bound()
                    << ", pos_lower_bound = " << startpos_lower_bound
                    << ", pos_upper_bound = " << startpos_upper_bound
                    << ", rate_lower_bound = "
                    << _interface->get_rate_lower_bound()
                    << ", interpolation = " << result->get_gpos() << std::endl;
        }
        result->set_endpos(cmp(pos2_query, startpos_upper_bound) < 0
                               ? pos2_query
                               : startpos_upper_bound);
        result->set_rate(_interface->get_rate_lower_bound());
        return;
      } else {
        // query is greater than both positions; increment
        if (verbose) {
          std::cout << "past current loaded window, incrementing range"
                    << std::endl;
        }
        _interface->get();
      }
    } else if (query_vs_lower_bound == EQUAL &&
               query_vs_upper_bound == LESS_THAN) {
      // gpos extension beyond end of range
      mpf_class gpos_interpolated;
      // certain genetic maps terminate their regions with a non-zero rate
      // window. for those situations, the extension needs to adjust for the end
      // position's interpolation from the reported rate of the start position
      // of the range.
      mpz_class startpos_lower_bound = _interface->get_startpos_lower_bound();
      mpz_class endpos_lower_bound = _interface->get_endpos_lower_bound();
      if (cmp(endpos_lower_bound, 0) == -1) {
        // there is no end position; the position is fixed
        gpos_interpolated = _interface->get_gpos_lower_bound();
        result->set_endpos(pos2_query);
        result->set_rate(_interface->get_rate_lower_bound());
      } else {
        // there is an end position; partial interpolation is required,
        // though if the reported rate is 0, this will be the same
        // as the above condition
        gpos_interpolated =
            _interface->get_gpos_lower_bound() +
            ((cmp(pos1_query, endpos_lower_bound) == -1 ? pos1_query
                                                        : endpos_lower_bound) -
             startpos_lower_bound) /
                mb_adjustment * _interface->get_rate_lower_bound();
        result->set_endpos(cmp(pos1_query, endpos_lower_bound) == -1
                               ? endpos_lower_bound
                               : pos2_query);
        result->set_rate(cmp(pos1_query, endpos_lower_bound) == -1
                             ? _interface->get_rate_lower_bound()
                             : 0.0);
      }
      if (verbose) {
        std::cout << "\tchromosome beyond end of range, setting to "
                  << gpos_interpolated << std::endl;
      }
      result->set_gpos(gpos_interpolated);
      return;
    } else if (query_vs_upper_bound == LESS_THAN) {
      // no estimate for relevant chromosome; set to 0?
      if (verbose) {
        std::cout << "\tno estimate for entire chromosome, setting to 0"
                  << std::endl;
      }
      result->set_gpos(0.0);
      result->set_endpos(pos2_query);
      result->set_rate(0.0);
      return;
    } else if (query_vs_lower_bound == GREATER_THAN) {
      if (verbose) {
        std::cout << "\tgenetic map range needs incrementing" << std::endl;
      }
      _interface->get();
      if (query_vs_upper_bound == GREATER_THAN) {
        if (verbose) {
          std::cout << "\t\ttwice!" << std::endl;
        }
        _interface->get();
      }
    } else {
      // sort error
      throw std::runtime_error(
          "genetic_map::query: an impossible condition was encountered."
          " This likely means that your input queries or genetic map "
          "are unsorted. For now, the only solution to this issue is "
          "to sort your input data; in the future, there may be an "
          "option to have the program handle this for you, at the cost "
          "of RAM.");
    }
  }
  query_vs_upper_bound = chromosome_compare(chr_query, chr_upper_bound);
  if (_interface->eof()) {
    if (query_vs_upper_bound == EQUAL) {
      // For variants at the end of the chromosome's rate information,
      // what is the correct course of action? This current value sets
      // out-of-range recombination to 0; this may be swapped out for
      // linear interpolation using the last known good rate of the
      // chromosome, which might be more biological but violates the
      // traditional convention of analyses extending beyond a model's
      // estimation range.
      mpf_class gpos_interpolated;
      mpz_class startpos_upper_bound = _interface->get_startpos_upper_bound();
      mpz_class endpos_upper_bound = _interface->get_endpos_upper_bound();
      // certain genetic maps terminate their regions with a non-zero rate
      // window. for those situations, the extension needs to adjust for the end
      // position's interpolation from the reported rate of the start position
      // of the range.
      if (cmp(endpos_upper_bound, 0) == -1) {
        // there is no end position; the position is fixed
        gpos_interpolated = _interface->get_gpos_upper_bound();
        result->set_endpos(pos2_query);
        result->set_rate(_interface->get_rate_upper_bound());
      } else {
        // there is an end position; partial interpolation is required,
        // though if the reported rate is 0, this will be the same
        // as the above condition
        gpos_interpolated =
            _interface->get_gpos_upper_bound() +
            ((cmp(pos1_query, endpos_upper_bound) == -1 ? pos1_query
                                                        : endpos_upper_bound) -
             startpos_upper_bound) /
                mb_adjustment * _interface->get_rate_upper_bound();
        result->set_endpos(cmp(pos1_query, endpos_upper_bound) == -1 &&
                                   cmp(pos2_query, 0) >= 0 &&
                                   cmp(pos2_query, endpos_upper_bound) != -1
                               ? endpos_upper_bound
                               : pos2_query);
        result->set_rate(cmp(pos1_query, endpos_upper_bound) == -1
                             ? _interface->get_rate_upper_bound()
                             : 0.0);
      }
      if (verbose) {
        std::cout
            << "\tfinal value, beyond boundary of loaded data; setting to "
            << gpos_interpolated << std::endl;
      }
      result->set_gpos(gpos_interpolated);
      return;
    } else {
      // We may eventually want to let the user specify that this should
      // fall back to the naive estimator.
      if (verbose) {
        std::cout << "\tfinal loaded data don't match chromosome, "
                  << "setting to 0" << std::endl;
      }
      result->set_gpos(0.0);
      result->set_endpos(pos2_query);
      result->set_rate(0.0);
      return;
    }
  } else {
    throw std::runtime_error(
        "genetic_map::query: an impossible condition was encountered. This "
        "likely means that your input queries or genetic map are unsorted. "
        "For now, the only solution to this issue is to sort your input "
        "data; in the future, there may be an option to have the program "
        "handle this for you, at the cost of RAM.");
  }
}

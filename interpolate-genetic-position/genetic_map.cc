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
  _gzinput = NULL;
  _buffer = 0;
  _buffer_size = 1000;
  _ft = UNKNOWN;
  _chr_old = "";
  _pos_old = 0;
  _gpos_old = 0.0;
  _rate_old = 0.0;
  _chr_new = "";
  _pos_new = 0;
  _gpos_new = 0.0;
  _rate_new = 0.0;

  _buffer = new char[_buffer_size];
}
igp::genetic_map::genetic_map(const genetic_map &obj) {
  throw std::runtime_error(
      "genetic_map: copy constructor operation is invalid for this class");
}
igp::genetic_map::~genetic_map() throw() { close(); }
void igp::genetic_map::open(const std::string &filename, format_type ft) {
  _ft = ft;
  // function must chomp the expected header line of bolt format files
  if (filename.rfind(".gz") == filename.size() - 3) {
    _gzinput = gzopen(filename.c_str(), "rb");
    if (!_gzinput) {
      throw std::runtime_error("genetic_map: cannot read file \"" + filename +
                               "\"");
    }
    if (_ft == BOLT) {
      gzgets(_gzinput, _buffer, _buffer_size);
    }
  } else {
    std::string line = "";
    _input.open(filename.c_str());
    if (!_input.is_open()) {
      throw std::runtime_error("genetic_map: cannot read file \"" + filename +
                               "\"");
    }
    if (_ft == BOLT) {
      getline(_input, line);
    }
  }
  // Load the first two values, such that a valid range is available at the
  // beginning of iteration.
  get();
  get();
}
void igp::genetic_map::open(const char *filename, format_type ft) {
  open(std::string(filename), ft);
}
bool igp::genetic_map::get() {
  _chr_old = _chr_new;
  _pos_old = _pos_new;
  _gpos_old = _gpos_new;
  _rate_old = _rate_new;
  std::string line = "", pos1 = "", pos2 = "", gpos = "", rate = "";
  if (_input.is_open()) {
    if (_input.peek() == EOF) {
      return false;
    }
    getline(_input, line);
  } else {
    if (gzgets(_gzinput, _buffer, _buffer_size) == Z_NULL) {
      return false;
    }
    line = std::string(_buffer);
  }
  std::istringstream strm1(line);
  if (_ft == BOLT) {
    // BOLT format includes a first column chromosome indicator, without a "chr"
    // prefix independent of genome build. It includes 1-22 and X encoded as 23.
    if (!(strm1 >> _chr_new >> pos1 >> rate >> gpos)) {
      throw std::runtime_error(
          "genetic_map::get: cannot parse BOLT ratefile line \"" + line + "\"");
    }
    _pos_new = pos1;
    _gpos_new = gpos;
    _rate_new = rate;
  } else if (_ft == UCSC) {
    // UCSC format includes rate but not genetic position itself. Assuming the
    // first position in a rate file is 0 genetic position, as is conventional,
    // we need to track the accumulated genetic position each time we get a new
    // entry from the file.
    if (!(strm1 >> _chr_new >> pos1 >> pos2 >> rate)) {
      throw std::runtime_error(
          "genetic_map::get: cannot parse UCSC bedfile line \"" + line + "\"");
    }
    _pos_new = pos1;
    _rate_new = rate;
    _gpos_new = _gpos_old + _rate_new * (_pos_new - _pos_old);
  } else if (_ft == UNKNOWN) {
    throw std::runtime_error("genetic_map::get: format is unset");
  } else {
    throw std::runtime_error(
        "genetic_map::get: unrecognized file format (check format_type enum)");
  }
  return true;
}
int igp::genetic_map::chromosome_to_integer(const std::string &chr) const {
  int rep = 0;
  std::string stripped_chr = chr;
  if (stripped_chr.find("chr") == 0) {
    stripped_chr = stripped_chr.substr(3);
  }
  if (!stripped_chr.compare("X")) {
    rep = 23;
  } else if (!stripped_chr.compare("Y")) {
    rep = 24;
  } else if (!stripped_chr.compare("MT") || !stripped_chr.compare("M")) {
    rep = 26;
  } else {
    std::istringstream strm1(stripped_chr);
    if (!(strm1 >> rep)) {
      throw std::runtime_error(
          "chromosome_to_integer: unhandled chromosome code \"" + chr + "\"");
    }
  }
  return rep;
}
igp::direction igp::genetic_map::chromosome_compare(
    const std::string &chr1, const std::string &chr2) const {
  unsigned int1 = chromosome_to_integer(chr1);
  unsigned int2 = chromosome_to_integer(chr2);
  if (int1 == int2) return EQUAL;
  return int1 < int2 ? LESS_THAN : GREATER_THAN;
}
void igp::genetic_map::query(const std::string &chr_query,
                             const mpz_class &pos_query, bool verbose,
                             mpf_class *gpos_interpolated) {
  // This comparison must cautiously detect whether a variant falls before,
  // after, or within the interval indicated by the currently loaded set of
  // genetic positions. Certain permutations of this comparison will indicate
  // that sorting assumptions on input data have been violated. Additionally,
  // the test for chromosome equality must be performed carefully to deal with
  // the many permutations of chromosome codes that are encountered in datasets.
  direction query_vs_old = EQUAL;
  direction query_vs_new = EQUAL;
  if (verbose) {
    std::cout << "query: chr is " << chr_query << ", pos is " << pos_query
              << std::endl;
  }
  while (!eof()) {
    query_vs_old = chromosome_compare(chr_query, _chr_old);
    query_vs_new = chromosome_compare(chr_query, _chr_new);
    if (query_vs_old == EQUAL && query_vs_new == EQUAL) {
      if (verbose) {
        std::cout << "\tsame chromosome, position comparing" << std::endl;
      }
      // proceed to position comparison
      int query_pos_vs_old = cmp(pos_query, _pos_old);
      int query_pos_vs_new = cmp(pos_query, _pos_new);
      if (query_pos_vs_old == 0 && query_pos_vs_new == -1) {
        if (verbose) {
          std::cout << "\t\tmatches lower boundary exactly: " << _gpos_old
                    << std::endl;
        }
        // exact lower boundary rate
        *gpos_interpolated = _gpos_old;
        return;
      } else if (query_pos_vs_old == 1 && query_pos_vs_new == 0) {
        // exact upper boundary rate
        if (verbose) {
          std::cout << "\t\tmatches upper boundary exactly: " << _gpos_new
                    << std::endl;
        }
        *gpos_interpolated = _gpos_new;
        return;
      } else if (query_pos_vs_old == -1 && query_pos_vs_new == -1) {
        // beginning of chromosome
        if (verbose) {
          std::cout << "\t\tbeginning of chromosome, before rate estimates "
                       "start; setting to 0"
                    << std::endl;
        }
        *gpos_interpolated = 0.0;
        return;
      } else if (query_pos_vs_old == 1 && query_pos_vs_new == -1) {
        // interpolate
        *gpos_interpolated = _gpos_old + (_pos_new - _pos_old) * _rate_old;
        if (verbose) {
          std::cout << "interpolated: gpos_old = " << _gpos_old
                    << ", pos_old = " << _pos_old << ", pos_new = " << _pos_new
                    << ", rate_old = " << _rate_old
                    << ", interpolation = " << *gpos_interpolated << std::endl;
        }
        return;
      } else if (query_pos_vs_old == 1 && query_pos_vs_new == 1) {
        // increment
        if (verbose) {
          std::cout << "past current loaded window, incrementing range"
                    << std::endl;
        }
        get();
      } else {
        // impossible condition
        throw std::runtime_error(
            "genetic_map::query: an impossible condition was encountered. This "
            "likely means that your input queries or genetic map are unsorted. "
            "For now, the only solution to this issue is to sort your input "
            "data; in the future, there may be an option to have the program "
            "handle this for you, at the cost of RAM.");
      }
    } else if (query_vs_old == EQUAL && query_vs_new == LESS_THAN) {
      // gpos extension beyond end of range
      if (verbose) {
        std::cout << "\tchromosome beyond end of range, setting to "
                  << _gpos_old << std::endl;
      }
      *gpos_interpolated = _gpos_old;
      return;
    } else if (query_vs_new == LESS_THAN) {
      // no estimate for relevant chromosome; set to 0?
      if (verbose) {
        std::cout << "\tno estimate for entire chromosome, setting to 0"
                  << std::endl;
      }
      *gpos_interpolated = 0.0;
      return;
    } else if (query_vs_old == GREATER_THAN) {
      if (verbose) {
        std::cout << "\tgenetic map range needs incrementing" << std::endl;
      }
      get();
      if (query_vs_new == GREATER_THAN) {
        if (verbose) {
          std::cout << "\t\ttwice!" << std::endl;
        }
        get();
      }
    } else {
      // sort error
      throw std::runtime_error(
          "genetic_map::query: an impossible condition was encountered. This "
          "likely means that your input queries or genetic map are unsorted. "
          "For now, the only solution to this issue is to sort your input "
          "data; in the future, there may be an option to have the program "
          "handle this for you, at the cost of RAM.");
    }
  }
  if (eof()) {
    if (query_vs_new == EQUAL) {
      // For variants at the end of the chromosome's rate information, what is
      // the correct course of action? This current value sets out-of-range
      // recombination to 0; this may be swapped out for linear interpolation
      // using the last known good rate of the chromosome, which might be more
      // biological but violates the traditional convention of analyses
      // extending beyond a model's estimation range.
      if (verbose) {
        std::cout
            << "\tfinal value, beyond boundary of loaded data; setting to "
            << _gpos_new << std::endl;
      }
      *gpos_interpolated = _gpos_new;
      return;
    } else {
      // We may eventually want to let the user specify that this should fall
      // back to the naive estimator.
      if (verbose) {
        std::cout << "\tfinal loaded data don't match chromosome, setting to 0"
                  << std::endl;
      }
      *gpos_interpolated = 0;
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
void igp::genetic_map::close() {
  if (_input.is_open()) {
    _input.close();
  }
  if (_gzinput) {
    gzclose(_gzinput);
    _gzinput = 0;
  }
  if (_buffer) {
    delete[] _buffer;
    _buffer = 0;
  }
}
bool igp::genetic_map::eof() {
  if (_input.is_open()) {
    return _input.peek() == EOF;
  }
  if (_gzinput) {
    return gzeof(_gzinput);
  }
  return false;
}

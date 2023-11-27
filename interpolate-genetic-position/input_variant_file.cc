/*!
 \file input_variant_file.cc
 \brief implementation for input variant file base class
 \author Lightning Auriga
 \copyright Released under the MIT License.
 Copyright 2023 Lightning Auriga
 */

#include "interpolate-genetic-position/input_variant_file.h"

namespace igp = interpolate_genetic_position;

igp::base_input_variant_file::base_input_variant_file() {}
igp::base_input_variant_file::~base_input_variant_file() throw() {}

igp::input_variant_file::input_variant_file()
    : igp::base_input_variant_file::base_input_variant_file(),
      _gzinput(0),
      _sr(0),
      _fallback(0),
      _buffer(0),
      _buffer_size(1000),
      _chr_index(0),
      _pos1_index(0),
      _pos2_index(-1),
      _gpos_index(0),
      _chr(""),
      _varid(""),
      _pos1(0),
      _pos2(-1),
      _base0(false),
      _vcf_eof(false) {
  _buffer = new char[_buffer_size];
}

igp::input_variant_file::~input_variant_file() throw() {
  close();
  if (_buffer) {
    delete[] _buffer;
    _buffer = 0;
  }
}

void igp::input_variant_file::open(const std::string &filename) {
  // this needs to be updated to catch vcfs
  if (filename.rfind(".vcf.gz") == filename.size() - 7 ||
      filename.rfind(".vcf") == filename.size() - 4 ||
      filename.rfind(".bcf") == filename.size() - 4) {
    _sr = bcf_sr_init();
    hts_set_log_level(HTS_LOG_OFF);
    if (!bcf_sr_add_reader(_sr, filename.c_str())) {
      throw std::runtime_error("input_variant_file: " +
                               std::string(bcf_sr_strerror(_sr->errnum)));
    }
    hts_set_log_level(HTS_LOG_WARNING);
  } else if (filename.rfind(".gz") == filename.size() - 3) {
    _gzinput = gzopen(filename.c_str(), "rb");
    if (!_gzinput) {
      throw std::runtime_error("input_variant_file: cannot open gzfile \"" +
                               filename + "\"");
    }
  } else if (!filename.empty()) {
    _input.open(filename.c_str());
    if (!_input.is_open()) {
      throw std::runtime_error("input_variant_file: cannot open file \"" +
                               filename + "\"");
    }
  }  // otherwise, filename is empty, and the object will probe std::cin
}

void igp::input_variant_file::close() {
  if (_input.is_open()) {
    _input.close();
    _input.clear();
  }
  if (_gzinput) {
    gzclose(_gzinput);
    _gzinput = 0;
  }
  if (_sr) {
    bcf_sr_destroy(_sr);
    _sr = 0;
    _vcf_eof = false;
  }
}

void igp::input_variant_file::set_fallback_stream(std::istream *ptr) {
  _fallback = ptr;
}

std::istream *igp::input_variant_file::get_fallback_stream() const {
  if (!_fallback) {
    throw std::runtime_error(
        "ivf::get_fallback_stream: called on NULL pointer");
  }
  return _fallback;
}

void igp::input_variant_file::set_format_parameters(
    unsigned chr_index, unsigned pos1_index, int pos2_index,
    unsigned gpos_index, bool base0, unsigned n_tokens) {
  _chr_index = chr_index;
  _pos1_index = pos1_index;
  _pos2_index = pos2_index;
  _gpos_index = gpos_index;
  _base0 = base0;
  _line_contents.resize(n_tokens, "");
}

bool igp::input_variant_file::get_variant() {
  std::string line = "", catcher = "";
  // as best as possible, check streams for closed status
  if (eof()) {
    return false;
  }
  // vcf input only: extract fields with htslib and return
  if (_sr) {
    if (!bcf_sr_next_line(_sr)) {
      _vcf_eof = true;
      return false;
    }
    // use htslib internal accessors, and skip downstream logic
    // for other filetypes
    _chr = std::string(
        bcf_seqname_safe(bcf_sr_get_header(_sr, 0), bcf_sr_get_line(_sr, 0)));
    bcf_unpack(bcf_sr_get_line(_sr, 0), BCF_UN_STR);
    _varid = std::string(bcf_sr_get_line(_sr, 0)->d.id);
    _pos1 = bcf_sr_get_line(_sr, 0)->pos - 1;
    return true;
  }
  // other stream types are text line parsers of various kinds
  if (_input.is_open()) {
    getline(_input, line);
  } else if (_gzinput) {
    if (gzgets(_gzinput, _buffer, _buffer_size) == Z_NULL) {
      return false;
    }
    line = std::string(_buffer);
  } else {
    getline(*get_fallback_stream(), line);
  }

  std::istringstream strm1(line);
  for (unsigned i = 0; i < _line_contents.size(); ++i) {
    if (!(strm1 >> _line_contents.at(i))) {
      throw std::runtime_error("get_variant: insufficient tokens: \"" + line +
                               "\"");
    }
  }
  _chr = _line_contents.at(_chr_index);
  _pos1 = _line_contents.at(_pos1_index);
  if (_base0) {
    _pos1 = _pos1 + 1;
  }
  if (_pos2_index >= 0) {
    _pos2 = _line_contents.at(static_cast<unsigned>(_pos2_index));
    if (_base0) {
      _pos2 = _pos2 + 1;
    }
  }

  return true;
}

const std::string &igp::input_variant_file::get_chr() const { return _chr; }

const mpz_class &igp::input_variant_file::get_pos1() const { return _pos1; }

const mpz_class &igp::input_variant_file::get_pos2() const { return _pos2; }

const std::string &igp::input_variant_file::get_varid() const { return _varid; }

const std::vector<std::string> &igp::input_variant_file::get_line_contents()
    const {
  return _line_contents;
}

bool igp::input_variant_file::eof() {
  if (_input.is_open()) {
    return _input.peek() == EOF;
  }
  if (_gzinput) {
    return gzeof(_gzinput);
  }
  if (_sr) {
    return _vcf_eof;
  }
  return get_fallback_stream()->peek() == EOF;
}

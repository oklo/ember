#pragma once
#include <array>
#include <charconv>
#include <cmath>
#include <istream>
#include <stdexcept>
#include <string>

namespace ember::detail {
// Numeric table bodies use an ASCII, whitespace-separated decimal format.
// Parse a plane in one buffer, avoiding locale work for every coefficient.
class NumericTableData {
  std::string storage_;
  const char* next_;
  const char* end_;
  void space() {
    while(next_!=end_ && (*next_==' ' || *next_=='\n' || *next_=='\r'
          || *next_=='\t' || *next_=='\v' || *next_=='\f'))++next_;
  }
public:
  explicit NumericTableData(std::istream& in) {
    std::array<char,65536> block;
    while(in.read(block.data(),static_cast<std::streamsize>(block.size())) || in.gcount())
      storage_.append(block.data(),static_cast<std::size_t>(in.gcount()));
    if(in.bad() || !in.eof())throw std::runtime_error("numeric table: input read failed");
    next_=storage_.data();end_=next_+storage_.size();
  }
  NumericTableData(const NumericTableData&)=delete;
  NumericTableData& operator=(const NumericTableData&)=delete;
  template<class Number> Number read() {
    space();if(next_==end_)throw std::runtime_error("numeric table: missing value");
    if(*next_=='+') {
      ++next_;
      if(next_==end_ || *next_=='-' || *next_=='+')
        throw std::runtime_error("numeric table: invalid leading sign");
    }
    Number value{};const auto parsed=std::from_chars(next_,end_,value);
    if(parsed.ec!=std::errc{} || parsed.ptr==next_)
      throw std::runtime_error("numeric table: invalid or out-of-range number");
    next_=parsed.ptr;const auto stop=next_;space();
    if(next_==stop && next_!=end_)throw std::runtime_error("numeric table: invalid number separator");
    if(!std::isfinite(value))throw std::runtime_error("numeric table: nonfinite number");
    return value;
  }
  void finish() {
    space();if(next_!=end_)throw std::runtime_error("numeric table: trailing data");
  }
};
} // namespace ember::detail

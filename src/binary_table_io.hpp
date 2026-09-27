#pragma once
#include <bit>
#include <cmath>
#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace ember::detail::binary {
static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559);
inline void write_u64(std::ostream& out,std::uint64_t value) {
  if constexpr(std::endian::native!=std::endian::little)value=std::byteswap(value);
  out.write(reinterpret_cast<const char*>(&value),8);
  if(!out)throw std::runtime_error("binary EOS: write failed");
}
inline std::uint64_t read_u64(std::istream& in) {
  std::uint64_t value{};in.read(reinterpret_cast<char*>(&value),8);
  if(!in)throw std::runtime_error("binary EOS: truncated integer");
  if constexpr(std::endian::native!=std::endian::little)value=std::byteswap(value);
  return value;
}
inline void write_text(std::ostream& out,const std::string& text) {
  write_u64(out,text.size());out.write(text.data(),static_cast<std::streamsize>(text.size()));
  if(!out)throw std::runtime_error("binary EOS: text write failed");
}
inline std::string read_text(std::istream& in) {
  const auto size=read_u64(in);
  if(size==0 || size>65536)throw std::runtime_error("binary EOS: invalid text length");
  std::string text(static_cast<std::size_t>(size),'\0');
  in.read(text.data(),static_cast<std::streamsize>(size));
  if(!in)throw std::runtime_error("binary EOS: truncated text");
  return text;
}
inline void write_values(std::ostream& out,std::span<const double> values) {
  if constexpr(std::endian::native==std::endian::little)
    out.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(values.size_bytes()));
  else for(double value:values)write_u64(out,std::bit_cast<std::uint64_t>(value));
  if(!out)throw std::runtime_error("binary EOS: value write failed");
}
inline std::vector<double> read_values(std::istream& in,std::size_t size) {
  std::vector<double> values(size);
  in.read(reinterpret_cast<char*>(values.data()),static_cast<std::streamsize>(size*sizeof(double)));
  if(!in)throw std::runtime_error("binary EOS: truncated values");
  for(double& value:values) {
    if constexpr(std::endian::native!=std::endian::little)
      value=std::bit_cast<double>(std::byteswap(std::bit_cast<std::uint64_t>(value)));
    if(!std::isfinite(value))throw std::runtime_error("binary EOS: nonfinite value");
  }
  return values;
}
} // namespace ember::detail::binary

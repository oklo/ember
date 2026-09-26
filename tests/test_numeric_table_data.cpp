#include "../src/numeric_table_data.hpp"
#include <bit>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <chrono>

int main(int argc,char** argv) {
  using ember::detail::NumericTableData;
  std::size_t checked=0;
  auto compare=[&](const std::string& body) {
    std::istringstream slow(body),raw(body);NumericTableData fast(raw);double expected;
    while(slow>>expected) {
      const auto actual=fast.read<double>();
      if(std::bit_cast<std::uint64_t>(expected)!=std::bit_cast<std::uint64_t>(actual))
        throw std::runtime_error("numeric parser changed coefficient bits");
      ++checked;
    }
    fast.finish();
  };
  std::mt19937_64 generator(7241);std::ostringstream values;values<<std::setprecision(17);
  for(int i=0;i<10000;++i) {
    const auto bits=generator();const double x=std::bit_cast<double>(bits);
    // The formatted-stream reference rejects some underflowing subnormals.
    // Compare the finite normal coefficients used by these source tables.
    if(std::isnormal(x) || x==0)values<<x<<'\n';
  }
  compare(values.str());compare("+1.25\t-0 0\r\n1E-200 .5 5. 1e+300\v3\f4");
  for(const std::string bad:{"", "nan", "inf", "1e400", "1e-400", "1.2bad", "1e", "+", "+-1", "++1", "--1", "0x12"}) {
    bool rejected=false;try{std::istringstream raw(bad);NumericTableData d(raw);d.read<double>();}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("invalid input accepted: "+bad);
    ++checked;
  }
  {std::istringstream raw("0 1 2");NumericTableData d(raw);if(d.read<int>()!=0 || d.read<int>()!=1) return 1;
   bool rejected=false;try{d.finish();}catch(const std::exception&){rejected=true;}if(!rejected)return 1;}
  if(argc==2) {
    std::ifstream in(argv[1]);if(!in)throw std::runtime_error("cannot open test table");
    std::string line;while(std::getline(in,line) && line!="data"){}
    const std::string body(std::istreambuf_iterator<char>(in),{});
    const auto start=std::chrono::steady_clock::now();compare(body);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::printf("real table body %zu bytes, bitwise comparison %.4g s\n",body.size(),seconds);
  }
  std::printf("%zu numeric values/checks passed with identical coefficient bits\n",checked);
}

#include "ember/eos_variable_metal.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ember {
void VariableMetalHelmholtzEos::pack_binary(const std::filesystem::path& source,
                                           const std::filesystem::path& destination) {
  namespace fs=std::filesystem;
  const fs::path pending=destination.string()+".pending";
  if(fs::exists(destination) || fs::exists(pending))
    throw std::invalid_argument("binary EOS: output already exists");
  std::ifstream in(source);std::string magic;int version{};in>>magic>>version;
  if(!in || magic!="EMBER_VARIABLE_METAL_HELMHOLTZ" || (version!=1 && version!=2))
    throw std::runtime_error("binary EOS: expected a text variable-metal family");
  std::array<std::vector<double>,3> axes;
  const std::array<std::string,3> labels{"metals","hydrogen_share","helium3_share"};
  const std::array<std::size_t,3> minimum{version==2?4u:3u,4,3},maximum{version==2?32u:4u,100,4};
  std::size_t count=1;
  for(std::size_t k=0;k<axes.size();++k) {
    std::string label;std::size_t n{};in>>label>>n;
    if(!in || label!=labels[k] || n<minimum[k] || n>maximum[k])
      throw std::runtime_error("binary EOS: invalid family axis");
    axes[k].resize(n);count*=n;
    for(std::size_t i=0;i<n;++i) {
      auto& x=axes[k][i];in>>x;
      if(!in || !std::isfinite(x) || x<0 || x>1 || (i && x<=axes[k][i-1]))
        throw std::runtime_error("binary EOS: invalid family coordinate");
    }
  }
  if(axes[0].back()>=1)throw std::runtime_error("binary EOS: material needs H or helium");
  try {
    std::ofstream out(pending,std::ios::binary);out.exceptions(std::ios::badbit|std::ios::failbit);
    out<<std::setprecision(std::numeric_limits<double>::max_digits10)
       <<"EMBER_VARIABLE_METAL_HELMHOLTZ_BINARY "<<version<<'\n';
    for(std::size_t k=0;k<axes.size();++k) {
      out<<labels[k]<<' '<<axes[k].size();for(double x:axes[k])out<<' '<<x;out<<'\n';
    }
    out<<"planes_binary_v1\n";
    // Read one unmodified plane at a time. The runtime family subtracts ionic
    // mixing and constructs its slopes after loading, identically in both formats.
    for(std::size_t i=0;i<count;++i) {
      std::string file;in>>std::quoted(file);
      if(!in || file.empty())throw std::runtime_error("binary EOS: missing source plane");
      HelmholtzTableEos plane(source.parent_path()/file,HelmholtzTableEos::Mixture::allow_documented_proxy);
      plane.write_binary(out);
    }
    std::string extra;if(in>>extra)throw std::runtime_error("binary EOS: trailing source manifest data");
    out.flush();out.close();fs::rename(pending,destination);
  }catch(...) {
    std::error_code ignored;fs::remove(pending,ignored);throw;
  }
}
} // namespace ember

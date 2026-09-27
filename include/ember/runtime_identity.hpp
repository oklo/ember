#pragma once
#include "ember/evolution_checkpoint.hpp"

namespace ember::driver {

// File locations are provenance, not physics. Bind each input's role and
// family coordinates to its bytes so a relocated dataset remains restartable.
class RuntimeIdentity {
  std::map<fs::path,std::string> files_;
public:
  Identities values{{"identity_schema","ember-lifetime-inputs-v1"}};

  void file(const std::string& role,const fs::path& path) {
    const auto canonical=fs::canonical(path);
    auto found=files_.find(canonical);
    if(found==files_.end())found=files_.emplace(canonical,file_identity(canonical)).first;
    values[role]=found->second;
  }

  void number(const std::string& role,double value) {
    if(!std::isfinite(value))throw std::invalid_argument("nonfinite runtime identity: "+role);
    std::ostringstream text;text.imbue(std::locale::classic());
    text<<std::setprecision(std::numeric_limits<double>::max_digits10)<<value;
    values[role]=text.str();
  }

  void family(const std::string& role,const fs::path& path,bool eos) {
    std::ifstream in(path);std::string magic,label,name;int version{};std::size_t count=1;
    in>>magic>>version;
    if(!in)throw std::runtime_error("missing runtime family: "+role);
    values[role+".format"]=magic+" "+std::to_string(version);
    if(eos && magic=="EMBER_VARIABLE_METAL_HELMHOLTZ_BINARY" && (version==1 || version==2)) {
      file(role+".archive",path);return;
    }
    if(eos) {
      if(magic!="EMBER_VARIABLE_METAL_HELMHOLTZ" || (version!=1 && version!=2))
        throw std::runtime_error("invalid EOS identity format");
      for(const auto* axis:{"metals","hydrogen_share","helium3_share"}) {
        std::size_t n{};in>>label>>n;
        if(!in || label!=axis || n<2 || n>100)throw std::runtime_error("invalid EOS identity axis");
        count*=n;number(role+"."+axis+".size",static_cast<double>(n));
        for(std::size_t i=0;i<n;++i)
          number(role+"."+axis+"."+std::to_string(i),read_representable_double(in));
      }
    }else {
      in>>count>>label>>name;
      if(!in || magic!="EMBER_OPACITY_MIXTURE" || version!=1 || count<2 || count>100
          || (label!="logR" && label!="logRho"))
        throw std::runtime_error("invalid opacity identity family");
      values[role+".density_axis"]=label;
      values[role+".source"]=name;
      number(role+".size",static_cast<double>(count));
    }
    for(std::size_t i=0;i<count;++i) {
      const auto plane=role+".plane."+std::to_string(i);
      if(!eos)number(plane+".metallicity",read_representable_double(in));
      in>>std::quoted(name);
      if(!in || name.empty())throw std::runtime_error("missing family identity source");
      file(plane,path.parent_path()/name);
    }
    if(in>>name)throw std::runtime_error("trailing runtime family data");
  }
};
} // namespace ember::driver

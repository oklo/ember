#pragma once
#include "ember/model.hpp"
#include "ember/energy_grid.hpp"
#include "ember/nuclear_cn.hpp"
#include <algorithm>
#include <charconv>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ember::driver {
namespace fs=std::filesystem;
inline double read_representable_double(std::istream& in) {
  std::string word;if(!(in>>word))throw std::runtime_error("missing floating checkpoint token");
  double value{};const auto parsed=std::from_chars(word.data(),word.data()+word.size(),value);
  if(parsed.ec!=std::errc{} || parsed.ptr!=word.data()+word.size() || !std::isfinite(value))
    throw std::runtime_error("invalid or unrepresentable checkpoint number");
  return value;
}
using Identities=std::map<std::string,std::string>;
using Selections=std::array<std::string,5>;

// FNV-1a detects accidental input changes. It is not a cryptographic digest;
// the snapshot runner separately records SHA-256 provenance. Hash the actual
// executable and table bytes, never their modification times.
inline std::string file_identity(const fs::path& path) {
  std::ifstream in(path,std::ios::binary);
  if(!in)throw std::runtime_error("cannot fingerprint checkpoint input: "+path.string());
  std::uint64_t value=14695981039346656037ULL;
  std::array<char,65536> buffer;
  while(in) {
    in.read(buffer.data(),buffer.size());
    for(std::streamsize i=0;i<in.gcount();++i) {
      value^=static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]);
      value*=1099511628211ULL;
    }
  }
  if(!in.eof())throw std::runtime_error("cannot read checkpoint input: "+path.string());
  std::ostringstream text;text<<std::hex<<std::setw(16)<<std::setfill('0')<<value;
  return text.str();
}

inline Identities input_identities(const fs::path& executable,const fs::path& data,
                                  const std::string& atmosphere,const std::string& eos,
                                  const fs::path& opacity_directory) {
  std::set<fs::path> paths;
  auto directory=[&](const fs::path& path) {
    for(const auto& entry:fs::directory_iterator(path))
      if(entry.is_regular_file() && entry.path().extension()==".dat")paths.insert(fs::canonical(entry.path()));
  };
  for(const auto* name:{"eos","opacity","conduction","atmosphere"})directory(data/name);
  // Explicit selection and every referenced plane are part of restart
  // identity, including planes outside the selected manifest's directory.
  for(const auto* name:{"aesopus21_gs98_mixture.dat","tops_gs98_mixture_low.dat","tops_gs98_mixture_high.dat"}) {
    const auto manifest=fs::canonical(opacity_directory/name);paths.insert(manifest);
    std::ifstream in(manifest);std::string magic,axis,label;int version{};std::size_t count{};
    in>>magic>>version>>count>>axis>>label;
    if(!in || magic!="EMBER_OPACITY_MIXTURE" || version!=1 || count<2 || count>100)
      throw std::runtime_error("invalid selected opacity manifest");
    for(std::size_t i=0;i<count;++i) {
      double z{};std::string file;in>>z>>std::quoted(file);
      if(!in || file.empty())throw std::runtime_error("invalid selected opacity source path");
      paths.insert(fs::canonical(manifest.parent_path()/file));
    }
  }
  if(atmosphere.starts_with("nongrey:"))paths.insert(fs::canonical(atmosphere.substr(8)));
  if(eos.starts_with("metal:")) {
    const auto family=fs::canonical(eos.substr(6));paths.insert(family);directory(family.parent_path());
  }
  Identities identities;
  // The binary can be copied to a new snapshot directory without changing
  // its identity. Data paths are kept explicit to avoid ambiguous families.
  identities["executable"]=file_identity(executable);
  identities["opacity_selection"]=fs::canonical(opacity_directory).string();
  for(const auto& path:paths)identities[path.string()]=file_identity(path);
  return identities;
}

struct Checkpoint {
  Model model;
  double next_dt{};
  std::size_t accepted{},rejected{};
  // Accepted total H1, He3 and metal face rates (g/s), used to recover the
  // previous material heat when lagging finite convection coefficients.
  std::vector<std::array<double,3>> metal_heat_rates{};
};

inline void check_state(const Checkpoint& state,std::size_t points,const Composition& composition) {
  const auto& model=state.model;
  (void)face_luminosities(model);
  if(model.size()!=points || model.m.size()!=points || model.comp.size()!=points
      || !std::isfinite(model.M) || model.M<=0 || !std::isfinite(model.age) || model.age<0
      || !std::isfinite(state.next_dt) || state.next_dt<=0)
    throw std::runtime_error("invalid checkpoint dimensions or clock");
  double previous_mass=0,previous_radius=0;
  for(std::size_t i=0;i<points;++i) {
    const auto& c=model.comp[i];const auto& y=model.y[i];
    if(!std::isfinite(model.m[i]) || model.m[i]<=previous_mass || model.m[i]>model.M
        || !std::isfinite(y.lnr) || !std::isfinite(y.lnrho) || !std::isfinite(y.lnT) || !std::isfinite(y.L)
        || !std::isfinite(model.r(i)) || model.r(i)<=previous_radius
        || !std::isfinite(model.rho(i)) || model.rho(i)<=0 || !std::isfinite(model.T(i)) || model.T(i)<=0
        || c.basis!=composition.basis || c.metal_inventory!=composition.metal_inventory
        || std::abs(c.sum()-1)>2e-12 || c.cn_mass_convention!=composition.cn_mass_convention
        || (c.cn_mass_convention==CNMassConvention::fixed_metal_proxy && std::abs(c.Z()-composition.Z())>2e-12))
      throw std::runtime_error("invalid checkpoint structure or composition");
    for(double x:c.X)if(!std::isfinite(x) || x<0 || x>1)throw std::runtime_error("invalid checkpoint abundance");
    if(c.cn_molality.has_value()!=composition.cn_molality.has_value())
      throw std::runtime_error("checkpoint CN inventory selection differs");
    if(c.cn_molality)(void)cn_physical_ledger(c,*c.cn_molality);
    previous_mass=model.m[i];previous_radius=model.r(i);
  }
  if(model.m.back()!=model.M)throw std::runtime_error("checkpoint surface mass differs");
  for(double value:model.Lsurf_hist)if(!std::isfinite(value))throw std::runtime_error("invalid checkpoint diagnostic");
  if(!state.metal_heat_rates.empty()) {
    if(!face_luminosities(model) || state.metal_heat_rates.size()+1!=points
        || composition.cn_mass_convention!=CNMassConvention::explicit_metal_mass || !composition.cn_molality)
      throw std::runtime_error("checkpoint material heat requires matching physical-metal faces");
    for(const auto& c:model.comp)if(c[Species::H2]!=0)
      throw std::runtime_error("checkpoint material heat with D requires an isotope-aware rate inventory");
    for(const auto& face:state.metal_heat_rates)for(double rate:face)
      if(!std::isfinite(rate))throw std::runtime_error("nonfinite checkpoint species heat rate");
  }
}

inline void write_checkpoint(const fs::path& path,const Checkpoint& state,
                             const Selections& selections,double tolerance,const Identities& identities) {
  if(state.model.comp.empty())throw std::runtime_error("empty checkpoint composition");
  check_state(state,state.model.size(),state.model.comp.front());
  const bool has_cn=state.model.comp.front().cn_molality.has_value();
  const bool physical_metals=state.model.comp.front().cn_mass_convention==CNMassConvention::explicit_metal_mass;
  const bool has_deuterium=std::any_of(state.model.comp.begin(),state.model.comp.end(),
      [](const Composition& c){return c[Species::H2]!=0;});
  const bool faces=face_luminosities(state.model);
  const bool material_heat=!state.metal_heat_rates.empty();
  if(physical_metals && !has_cn)throw std::runtime_error("physical-metal checkpoint requires explicit CN inventory");
  const fs::path temporary=path.string()+".pending";
  if(fs::exists(temporary))throw std::runtime_error("checkpoint temporary file already exists");
  std::ofstream out(temporary);out.imbue(std::locale::classic());
  out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<std::setprecision(std::numeric_limits<double>::max_digits10)<<"EMBER_EVOLUTION_CHECKPOINT "
     <<(material_heat?6:(faces?5:(has_deuterium?4:(physical_metals?3:(has_cn?2:1)))))<<'\n';
  if(faces)out<<has_cn<<' '<<physical_metals<<" volume_faces\n";
  else if(has_deuterium)out<<has_cn<<' '<<physical_metals<<'\n';
  for(const auto& value:selections)out<<std::quoted(value)<<'\n';
  out<<tolerance<<'\n'<<identities.size()<<'\n';
  for(const auto& [name,value]:identities)out<<std::quoted(name)<<' '<<std::quoted(value)<<'\n';
  const auto& model=state.model;
  out<<model.size()<<' '<<model.M<<' '<<model.age<<' '<<state.next_dt<<' '<<state.accepted<<' '<<state.rejected<<'\n';
  for(std::size_t i=0;i<model.size();++i) {
    const auto& y=model.y[i];const auto& c=model.comp[i];
    out<<model.m[i]<<' '<<y.lnr<<' '<<y.lnrho<<' '<<y.lnT<<' '<<y.L<<' '
       <<static_cast<int>(c.basis)<<' '<<static_cast<int>(c.metal_inventory);
    for(std::size_t j=0;j<(faces || has_deuterium?NSPEC:TABLE_NSPEC);++j)out<<' '<<c.X[j];
    if(has_cn)for(double ycn:*c.cn_molality)out<<' '<<ycn;
    out<<'\n';
  }
  out<<model.Lsurf_hist.size()<<'\n';
  for(double value:model.Lsurf_hist)out<<value<<'\n';
  if(material_heat) {
    out<<"METAL_HEAT_RATES "<<state.metal_heat_rates.size()<<'\n';
    for(const auto& face:state.metal_heat_rates)out<<face[0]<<' '<<face[1]<<' '<<face[2]<<'\n';
  }
  out<<"END\n";out.flush();out.close();
  fs::rename(temporary,path);
}

inline Checkpoint read_checkpoint(const fs::path& path,std::size_t expected_points,double expected_mass,
                                  const Composition& composition,const Selections& selections,
                                  double tolerance,const Identities& identities,
                                  LuminosityGrid expected_grid=LuminosityGrid::mass_nodes,
                                  bool expected_metal_heat_rates=false) {
  std::ifstream in(path);in.imbue(std::locale::classic());
  if(!in)throw std::runtime_error("cannot open checkpoint");
  std::string marker;int version{};in>>marker>>version;
  if(marker!="EMBER_EVOLUTION_CHECKPOINT" || version<1 || version>6)throw std::runtime_error("unsupported checkpoint format");
  if((version==6)!=expected_metal_heat_rates)
    throw std::runtime_error("checkpoint material-heat inventory differs");
  bool has_cn=version==2 || version==3,physical_metals=version==3;
  LuminosityGrid grid=LuminosityGrid::mass_nodes;
  if(version>=4) {
    int cn{},metals{};in>>cn>>metals;
    if(!in || (cn!=0 && cn!=1) || (metals!=0 && metals!=1) || (metals && !cn))
      throw std::runtime_error("invalid extended checkpoint inventory flags");
    has_cn=cn!=0;physical_metals=metals!=0;
    if(version>=5) {
      std::string placement;in>>placement;
      if(placement=="volume_faces")grid=LuminosityGrid::volume_faces;
      else if(placement!="mass_nodes")throw std::runtime_error("unknown checkpoint luminosity grid");
    }
  }
  if(grid!=expected_grid)throw std::runtime_error("checkpoint luminosity grid differs");
  Selections saved;for(auto& value:saved)in>>std::quoted(value);
  double saved_tolerance{};in>>saved_tolerance;
  if(saved!=selections || saved_tolerance!=tolerance)throw std::runtime_error("checkpoint physics or tolerances differ");
  std::size_t count{};in>>count;
  if(!in || count==0 || count>10000)throw std::runtime_error("invalid checkpoint input inventory");
  Identities recorded;
  for(std::size_t i=0;i<count;++i) {
    std::string name,value;in>>std::quoted(name)>>std::quoted(value);
    if(!recorded.emplace(name,value).second)throw std::runtime_error("duplicate checkpoint input identity");
  }
  if(!recorded.contains("executable"))throw std::runtime_error("checkpoint executable identity missing");
  for(const auto& [name,value]:recorded) {
    const auto found=identities.find(name);
    if(found==identities.end() || found->second!=value)
      throw std::runtime_error("checkpoint executable or input tables differ");
  }
  Checkpoint state;auto& model=state.model;std::size_t points{};
  model.luminosity_grid=grid;
  // These are lifetime counters, not per-invocation execution limits. Parse
  // unsigned values strictly: formatted extraction would accept a minus sign.
  std::string accepted,rejected;
  in>>points>>model.M>>model.age>>state.next_dt>>accepted>>rejected;
  auto counter=[](const std::string& value) {
    if(value.empty() || value.find_first_not_of("0123456789")!=std::string::npos)
      throw std::runtime_error("invalid checkpoint counter");
    std::size_t result=0;
    for(char digit:value) {
      const auto n=static_cast<std::size_t>(digit-'0');
      if(result>(std::numeric_limits<std::size_t>::max()-n)/10)
        throw std::runtime_error("checkpoint counter overflow");
      result=10*result+n;
    }
    return result;
  };
  state.accepted=counter(accepted);state.rejected=counter(rejected);
  if(points!=expected_points || model.M!=expected_mass)throw std::runtime_error("checkpoint mesh or mass differs");
  model.m.resize(points);model.y.resize(points);model.comp.resize(points);
  for(std::size_t i=0;i<points;++i) {
    auto& y=model.y[i];auto& c=model.comp[i];int basis{},inventory{};
    in>>model.m[i]>>y.lnr>>y.lnrho>>y.lnT>>y.L>>basis>>inventory;
    if(basis!=static_cast<int>(composition.basis) || inventory!=static_cast<int>(composition.metal_inventory))
      throw std::runtime_error("checkpoint abundance basis or elemental inventory differs");
    c.basis=composition.basis;c.metal_inventory=composition.metal_inventory;
    for(std::size_t j=0;j<(version>=4?NSPEC:TABLE_NSPEC);++j)c.X[j]=read_representable_double(in);
    if(has_cn) {
      CNAbundances ycn{};for(auto& value:ycn)value=read_representable_double(in);c.cn_molality=ycn;
    }
    if(physical_metals)c.cn_mass_convention=CNMassConvention::explicit_metal_mass;
  }
  in>>count;if(!in || count>20000)throw std::runtime_error("invalid checkpoint diagnostics length");
  model.Lsurf_hist.resize(count);for(auto& value:model.Lsurf_hist)in>>value;
  if(version==6) {
    in>>marker>>count;
    if(!in || marker!="METAL_HEAT_RATES" || count+1!=points)
      throw std::runtime_error("invalid checkpoint material-heat inventory");
    state.metal_heat_rates.resize(count);
    for(auto& face:state.metal_heat_rates)for(auto& rate:face)rate=read_representable_double(in);
  }
  in>>marker;if(!in || marker!="END")throw std::runtime_error("truncated checkpoint");
  in>>std::ws;if(!in.eof())throw std::runtime_error("unexpected trailing checkpoint content");
  check_state(state,points,composition);return state;
}
} // namespace ember::driver

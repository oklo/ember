#pragma once
#include "ember/model.hpp"
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
                                  const std::string& atmosphere,const std::string& eos) {
  std::set<fs::path> paths;
  auto directory=[&](const fs::path& path) {
    for(const auto& entry:fs::directory_iterator(path))
      if(entry.is_regular_file() && entry.path().extension()==".dat")paths.insert(fs::canonical(entry.path()));
  };
  for(const auto* name:{"eos","opacity","conduction","atmosphere"})directory(data/name);
  if(atmosphere.starts_with("nongrey:"))paths.insert(fs::canonical(atmosphere.substr(8)));
  if(eos.starts_with("metal:")) {
    const auto family=fs::canonical(eos.substr(6));paths.insert(family);directory(family.parent_path());
  }
  Identities identities;
  // The binary can be copied to a new snapshot directory without changing
  // its identity. Data paths are kept explicit to avoid ambiguous families.
  identities["executable"]=file_identity(executable);
  for(const auto& path:paths)identities[path.string()]=file_identity(path);
  return identities;
}

struct Checkpoint {
  Model model;
  double next_dt{};
  std::size_t accepted{},rejected{};
};

inline void check_state(const Checkpoint& state,std::size_t points,const Composition& composition) {
  const auto& model=state.model;
  if(model.size()!=points || model.m.size()!=points || model.comp.size()!=points
      || !std::isfinite(model.M) || model.M<=0 || !std::isfinite(model.age) || model.age<0
      || !std::isfinite(state.next_dt) || state.next_dt<=0 || state.accepted>10000 || state.rejected>101)
    throw std::runtime_error("invalid checkpoint dimensions, clock or counters");
  double previous_mass=0,previous_radius=0;
  for(std::size_t i=0;i<points;++i) {
    const auto& c=model.comp[i];const auto& y=model.y[i];
    if(!std::isfinite(model.m[i]) || model.m[i]<=previous_mass || model.m[i]>model.M
        || !std::isfinite(y.lnr) || !std::isfinite(y.lnrho) || !std::isfinite(y.lnT) || !std::isfinite(y.L)
        || !std::isfinite(model.r(i)) || model.r(i)<=previous_radius
        || !std::isfinite(model.rho(i)) || model.rho(i)<=0 || !std::isfinite(model.T(i)) || model.T(i)<=0
        || c.basis!=composition.basis || c.metal_inventory!=composition.metal_inventory
        || std::abs(c.sum()-1)>2e-12 || std::abs(c.Z()-composition.Z())>2e-12)
      throw std::runtime_error("invalid checkpoint structure or composition");
    for(double x:c.X)if(!std::isfinite(x) || x<0 || x>1)throw std::runtime_error("invalid checkpoint abundance");
    previous_mass=model.m[i];previous_radius=model.r(i);
  }
  if(model.m.back()!=model.M)throw std::runtime_error("checkpoint surface mass differs");
  for(double value:model.Lsurf_hist)if(!std::isfinite(value))throw std::runtime_error("invalid checkpoint diagnostic");
}

inline void write_checkpoint(const fs::path& path,const Checkpoint& state,
                             const Selections& selections,double tolerance,const Identities& identities) {
  const fs::path temporary=path.string()+".pending";
  if(fs::exists(temporary))throw std::runtime_error("checkpoint temporary file already exists");
  std::ofstream out(temporary);out.imbue(std::locale::classic());
  out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<std::setprecision(std::numeric_limits<double>::max_digits10)<<"EMBER_EVOLUTION_CHECKPOINT 1\n";
  for(const auto& value:selections)out<<std::quoted(value)<<'\n';
  out<<tolerance<<'\n'<<identities.size()<<'\n';
  for(const auto& [name,value]:identities)out<<std::quoted(name)<<' '<<std::quoted(value)<<'\n';
  const auto& model=state.model;
  out<<model.size()<<' '<<model.M<<' '<<model.age<<' '<<state.next_dt<<' '<<state.accepted<<' '<<state.rejected<<'\n';
  for(std::size_t i=0;i<model.size();++i) {
    const auto& y=model.y[i];const auto& c=model.comp[i];
    out<<model.m[i]<<' '<<y.lnr<<' '<<y.lnrho<<' '<<y.lnT<<' '<<y.L<<' '
       <<static_cast<int>(c.basis)<<' '<<static_cast<int>(c.metal_inventory);
    for(double x:c.X)out<<' '<<x;
    out<<'\n';
  }
  out<<model.Lsurf_hist.size()<<'\n';
  for(double value:model.Lsurf_hist)out<<value<<'\n';
  out<<"END\n";out.flush();out.close();
  fs::rename(temporary,path);
}

inline Checkpoint read_checkpoint(const fs::path& path,std::size_t expected_points,double expected_mass,
                                  const Composition& composition,const Selections& selections,
                                  double tolerance,const Identities& identities) {
  std::ifstream in(path);in.imbue(std::locale::classic());
  if(!in)throw std::runtime_error("cannot open checkpoint");
  std::string marker;int version{};in>>marker>>version;
  if(marker!="EMBER_EVOLUTION_CHECKPOINT" || version!=1)throw std::runtime_error("unsupported checkpoint format");
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
  in>>points>>model.M>>model.age>>state.next_dt>>state.accepted>>state.rejected;
  if(points!=expected_points || model.M!=expected_mass)throw std::runtime_error("checkpoint mesh or mass differs");
  model.m.resize(points);model.y.resize(points);model.comp.resize(points);
  for(std::size_t i=0;i<points;++i) {
    auto& y=model.y[i];auto& c=model.comp[i];int basis{},inventory{};
    in>>model.m[i]>>y.lnr>>y.lnrho>>y.lnT>>y.L>>basis>>inventory;
    if(basis!=static_cast<int>(composition.basis) || inventory!=static_cast<int>(composition.metal_inventory))
      throw std::runtime_error("checkpoint abundance basis or elemental inventory differs");
    c.basis=composition.basis;c.metal_inventory=composition.metal_inventory;
    for(auto& x:c.X)in>>x;
  }
  in>>count;if(!in || count>20000)throw std::runtime_error("invalid checkpoint diagnostics length");
  model.Lsurf_hist.resize(count);for(auto& value:model.Lsurf_hist)in>>value;
  in>>marker;if(!in || marker!="END")throw std::runtime_error("truncated checkpoint");
  in>>std::ws;if(!in.eof())throw std::runtime_error("unexpected trailing checkpoint content");
  check_state(state,points,composition);return state;
}
} // namespace ember::driver

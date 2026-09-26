#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ember::pms {
// A narrowly defined restart operation: append high-T/high-g nodes, or fill cells
// that were previously missing, without changing any old axis, supported data
// value or physical declaration.
// The atmosphere's normal parser independently validates each complete table.
inline std::size_t check_atmosphere_extension_files(
    const std::filesystem::path& old_path, const std::filesystem::path& new_path) {
  struct Table {
    std::vector<std::string> metadata;
    std::array<std::vector<double>,4> axes;
    std::vector<std::vector<double>> rows;
  };
  const auto read=[](const std::filesystem::path& path) {
    Table out;std::ifstream in(path);std::string line;
    if(!std::getline(in,line) || line!="EMBER_COMPOSITION_ATMOSPHERE 2")
      throw std::runtime_error("PMS extension requires atmosphere version 2");
    out.metadata.push_back(line);bool data=false;
    while(std::getline(in,line)) {
      if(line=="data"){data=true;break;}
      std::istringstream row(line);std::string key;row>>key;
      const std::array<std::string,4> names{"hydrogen","helium3","log_teff","log_g"};
      auto at=std::find(names.begin(),names.end(),key);
      if(at==names.end()){out.metadata.push_back(line);continue;}
      auto& axis=out.axes[at-names.begin()];std::size_t n=0;row>>n;
      if(!axis.empty() || n<2)throw std::runtime_error("invalid extension axis");
      axis.resize(n);for(double& v:axis)if(!(row>>v) || !std::isfinite(v))
        throw std::runtime_error("invalid extension axis value");
      std::string extra;if(row>>extra)throw std::runtime_error("extra extension axis token");
    }
    std::size_t expected=1;for(const auto& a:out.axes)expected*=a.size();
    if(!data || !expected)throw std::runtime_error("missing extension axes/data");
    while(std::getline(in,line)) {
      if(line.empty())continue;
      std::istringstream row(line);std::vector<double> values;double v;
      while(row>>v){if(!std::isfinite(v))throw std::runtime_error("nonfinite extension data");values.push_back(v);}
      if(!row.eof() || !((values.size()==1 && values[0]==0) || (values.size()==3 && values[0]==1)))
        throw std::runtime_error("invalid extension data row");
      out.rows.push_back(values);
    }
    if(out.rows.size()!=expected)throw std::runtime_error("extension row count mismatch");
    return out;
  };
  const auto old=read(old_path),next=read(new_path);
  if(old.metadata!=next.metadata || old.axes[0]!=next.axes[0] || old.axes[1]!=next.axes[1])
    throw std::runtime_error("atmosphere extension changes physics/composition");
  for(std::size_t k:{2u,3u})
    if(next.axes[k].size()<old.axes[k].size() ||
       !std::equal(old.axes[k].begin(),old.axes[k].end(),next.axes[k].begin()))
      throw std::runtime_error("atmosphere extension changes an old axis");
  std::size_t added=0;
  for(std::size_t h=0;h<next.axes[0].size();++h)
    for(std::size_t he=0;he<next.axes[1].size();++he)
      for(std::size_t t=0;t<next.axes[2].size();++t)
        for(std::size_t g=0;g<next.axes[3].size();++g) {
          const auto i=((h*next.axes[1].size()+he)*next.axes[2].size()+t)*next.axes[3].size()+g;
          if(t>=old.axes[2].size() || g>=old.axes[3].size()){added+=next.rows[i][0]==1;continue;}
          const auto j=((h*old.axes[1].size()+he)*old.axes[2].size()+t)*old.axes[3].size()+g;
          // A cell inside the old axes may go from MISSING to supported. That is an
          // addition, not a change: a missing cell carries no value, so no evolution
          // can have consumed one - a query there aborts the run instead. Supported
          // cells stay frozen, and un-supporting one still fails the comparison below.
          // This admits lowering a supported temperature floor, which the axes cannot
          // express on their own, since an old axis must remain a prefix of the new.
          if(old.rows[j].size()==1 && next.rows[i].size()==3){++added;continue;}
          if(next.rows[i]!=old.rows[j])throw std::runtime_error("atmosphere extension changes an old value/mask");
        }
  if(!added)throw std::runtime_error("atmosphere extension adds no supported states");
  return added;
}
} // namespace ember::pms

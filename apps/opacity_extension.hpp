#pragma once
#include "evolution_checkpoint.hpp"
#include <vector>

namespace ember::driver {
// A narrowly defined input extension: retain every original table entry and
// axis, adding only lower-H planes to the high-temperature TOPS family.
// This is not a general permission to change physics on a saved trajectory.
inline std::size_t check_opacity_extension(const fs::path& old_directory,const fs::path& new_directory) {
  struct Family {std::string axis,label;std::vector<double> z;std::vector<fs::path> files;};
  auto family=[](const fs::path& path) {
    Family out;std::ifstream in(path);std::string magic;int version{};std::size_t n{};
    in>>magic>>version>>n>>out.axis>>out.label;
    if(!in || magic!="EMBER_OPACITY_MIXTURE" || version!=1 || n<2 || n>100)
      throw std::runtime_error("opacity extension: invalid family manifest");
    for(std::size_t i=0;i<n;++i) {
      double z{};std::string file;in>>z>>std::quoted(file);
      if(!in || !std::isfinite(z) || z<=0 || z>=1 || (i && z<=out.z.back()) || file.empty())
        throw std::runtime_error("opacity extension: invalid source coordinate");
      out.z.push_back(z);out.files.push_back(fs::canonical(path.parent_path()/file));
    }
    if(in>>magic)throw std::runtime_error("opacity extension: trailing manifest data");
    return out;
  };
  struct Plane {double x{};std::vector<std::string> lines;};
  struct Table {std::size_t nt{},nr{};std::string density,temperature;std::vector<Plane> planes;};
  auto table=[](const fs::path& path) {
    Table out;std::ifstream in(path);std::string line;std::size_t nx{};
    std::getline(in,line);std::istringstream header(line);header>>nx>>out.nt>>out.nr;
    if(!header || nx<2 || nx>100 || out.nt<2 || out.nt>1000 || out.nr<2 || out.nr>1000)
      throw std::runtime_error("opacity extension: invalid source dimensions");
    std::getline(in,out.density);std::getline(in,out.temperature);
    for(std::size_t i=0;i<nx;++i) {
      Plane plane;std::getline(in,line);std::istringstream coordinates(line);double z{};
      coordinates>>plane.x>>z;
      if(!coordinates || !std::isfinite(plane.x) || plane.x<0 || plane.x>1 || !std::isfinite(z)
          || (i && plane.x<=out.planes.back().x))
        throw std::runtime_error("opacity extension: invalid hydrogen coordinates");
      plane.lines.push_back(line);
      for(std::size_t j=0;j<out.nt;++j) {std::getline(in,line);plane.lines.push_back(line);}
      if(!in)throw std::runtime_error("opacity extension: truncated source table");
      out.planes.push_back(std::move(plane));
    }
    in>>std::ws;if(!in.eof())throw std::runtime_error("opacity extension: trailing source data");
    return out;
  };
  std::size_t additions=0;
  for(const auto* name:{"aesopus21_gs98_mixture.dat","tops_gs98_mixture_low.dat","tops_gs98_mixture_high.dat"}) {
    const auto old=family(old_directory/name), next=family(new_directory/name);
    if(old.axis!=next.axis || old.label!=next.label || old.z!=next.z)
      throw std::runtime_error("opacity extension: source families differ");
    for(std::size_t i=0;i<old.files.size();++i) {
      if(std::string(name)!="tops_gs98_mixture_high.dat") {
        if(file_identity(old.files[i])!=file_identity(next.files[i]))
          throw std::runtime_error("opacity extension: cooler opacity changed");
        continue;
      }
      const auto a=table(old.files[i]),b=table(next.files[i]);
      if(a.nt!=b.nt || a.nr!=b.nr || a.density!=b.density || a.temperature!=b.temperature
          || b.planes.size()<a.planes.size())
        throw std::runtime_error("opacity extension: original opacity axes changed");
      const auto added=b.planes.size()-a.planes.size();additions+=added;
      for(std::size_t j=0;j<a.planes.size();++j)
        if(a.planes[j].lines!=b.planes[j+added].lines)
          throw std::runtime_error("opacity extension: original opacity entries changed");
    }
  }
  if(additions==0)throw std::runtime_error("opacity extension: no lower-hydrogen planes added");
  return additions;
}
} // namespace ember::driver

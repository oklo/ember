#pragma once
#include "ember/atmosphere_metal_piecewise.hpp"
#include "ember/atmosphere_metal_response.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>

namespace ember {

// Consecutive measured metal intervals share their physical reference
// atmosphere. Adding a lower interval leaves every higher interval intact.
// A single response file retains the original one-interval convention.
class MetalAtmosphereChain final : public Atmosphere {
public:
  MetalAtmosphereChain(const Eos& eos, const Atmosphere& reference,
                      const std::filesystem::path& path, double first_reference)
      : selected_(&reference) {
    std::ifstream in(path);
    std::string magic;
    in >> magic;
    std::vector<std::pair<double,std::filesystem::path>> entries;
    if (magic == "EMBER_METAL_ATMOSPHERE_CHAIN") {
      unsigned version{}, count{};
      in >> version >> count;
      if (!in || version != 1 || count < 1 || count > 16)
        throw std::runtime_error("MetalAtmosphereChain: invalid manifest");
      for (unsigned i=0; i<count; ++i) {
        double z{}; std::string file;
        in >> z >> std::quoted(file);
        if (!in || file.empty() || !std::isfinite(z) || z<=0 || z>=1)
          throw std::runtime_error("MetalAtmosphereChain: invalid interval");
        entries.emplace_back(z,path.parent_path()/file);
      }
      if (in >> magic)
        throw std::runtime_error("MetalAtmosphereChain: trailing manifest data");
    } else if (magic == "EMBER_METAL_ATMOSPHERE_RESPONSE") {
      entries.emplace_back(first_reference,path);
    } else throw std::runtime_error("MetalAtmosphereChain: unknown source format");
    double next=first_reference;
    for (const auto& [join,file] : entries) {
      std::ifstream source(file);
      std::string line; double reference_z=-1, source_z=-1;
      unsigned references=0, sources=0;
      while (std::getline(source,line)) {
        std::istringstream row(line); std::string key; row >> key;
        if (key=="reference_Z") {row >> reference_z; ++references;}
        if (key=="source_Z") {row >> source_z; ++sources;}
      }
      if (references!=1 || sources!=1 || !std::isfinite(reference_z+source_z) ||
          source_z<0 || !(source_z<reference_z) || std::abs(join-next)>1e-14 ||
          std::abs(join-reference_z)>1e-14)
        throw std::runtime_error("MetalAtmosphereChain: intervals do not meet");
      auto response=std::make_unique<MetalResponseAtmosphere>(eos,*selected_,file,
          MetalResponseAtmosphere::Approximation::separable_gs98_response,
          MetalResponseAtmosphere::Options{0});
      auto joined=std::make_unique<PiecewiseMetalAtmosphere>(*selected_,*response,join);
      selected_=joined.get();
      responses_.push_back(std::move(response)); joins_.push_back(std::move(joined));
      next=source_z;
    }
  }
  AtmosphereState eval(double T,double g,const Composition& c) const override {
    return selected_->eval(T,g,c);
  }
  const char* name() const override {return "consecutive measured metal-atmosphere intervals";}
private:
  const Atmosphere* selected_;
  std::vector<std::unique_ptr<MetalResponseAtmosphere>> responses_;
  std::vector<std::unique_ptr<PiecewiseMetalAtmosphere>> joins_;
};
} // namespace ember

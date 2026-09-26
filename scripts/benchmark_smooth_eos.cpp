#include "ember/eos_smooth_mixture.hpp"
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <vector>

int main(int argc,char** argv) {
  if(argc!=3)return 2;
  try {
    struct State {double T,rho;ember::Composition c;};std::vector<State> states;
    std::ifstream in(argv[2]);std::string line;std::getline(in,line);
    auto split=[](const std::string& s) {
      std::istringstream stream(s);std::vector<std::string> v;std::string item;
      while(std::getline(stream,item,','))v.push_back(item);return v;
    };
    const auto header=split(line);
    auto column=[&](const char* key) {
      for(std::size_t i=0;i<header.size();++i)if(header[i]==key)return i;
      throw std::runtime_error("missing profile column");
    };
    const auto it=column("temperature_K"),ir=column("density_g_cm3"),ix=column("X"),iy=column("Y3");
    while(std::getline(in,line)) {
      const auto v=split(line);const double x=std::stod(v[ix]),y=std::stod(v[iy]);
      auto c=ember::solar_scaled(x,.02);c.X[1]=y;c.X[2]-=y;
      c.basis=ember::AbundanceBasis::baryon_mass;c.metal_inventory=ember::MetalInventory::gs98;
      states.push_back({std::stod(v[it]),std::stod(v[ir]),c});
    }
    if(states.empty())throw std::runtime_error("empty profile");
    const auto mixture=ember::HelmholtzTableEos::Mixture::allow_documented_proxy;
    const auto begin=std::chrono::steady_clock::now();const auto cpu=std::clock();
    ember::MetalHelmholtzEos linear(argv[1],mixture);
    const double linear_cpu=double(std::clock()-cpu)/CLOCKS_PER_SEC;
    const auto middle=std::clock();ember::SmoothMetalHelmholtzEos smooth(argv[1],mixture);
    const double smooth_cpu=double(std::clock()-middle)/CLOCKS_PER_SEC;
    constexpr int repeats=200;
    std::cout<<std::setprecision(17)<<"{\"profile_zones\":"<<states.size()
      <<",\"calls_per_block\":"<<states.size()*repeats<<",\"linear_load_cpu_seconds\":"<<linear_cpu
      <<",\"smooth_load_cpu_seconds\":"<<smooth_cpu<<",\"blocks\":[";
    for(int block=0;block<4;++block) {
      const bool use_smooth=block==1 || block==2;
      const ember::Eos& eos=use_smooth?static_cast<const ember::Eos&>(smooth):linear;
      const auto start=std::chrono::steady_clock::now();const auto c0=std::clock();double checksum=0;
      for(int repeat=0;repeat<repeats;++repeat)for(const auto& s:states) {
        const auto e=eos.eval_with_derivatives(s.T,s.rho,s.c);
        checksum+=e.state.P/(s.rho*s.T)+e.state.cp+e.dcp_dlnT;
      }
      std::cout<<(block?",":"")<<"{\"mode\":\""<<(use_smooth?"smooth":"linear")
        <<"\",\"cpu_seconds\":"<<double(std::clock()-c0)/CLOCKS_PER_SEC
        <<",\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
        <<",\"checksum\":"<<checksum<<"}";
    }
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"],\"peak_rss_platform_units\":"<<usage.ru_maxrss
      <<",\"total_elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<"}\n";
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

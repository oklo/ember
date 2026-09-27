#include "ember/opacity_radiative.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
using namespace ember;
int main(int argc,char**argv) {
 if(argc!=2)return 2;
 try {
  const std::filesystem::path data=argv[1];std::ifstream in(data/"production/configuration.txt");
  std::map<std::string,std::string> cfg;std::string k,v;while(in>>k>>std::quoted(v))cfg[k]=v;
  const RadiativeOpacity::Tables files{data/"production"/cfg.at("opacity_low"),
      data/"production"/cfg.at("opacity_warm"),data/"production"/cfg.at("opacity_bridge"),
      data/"production"/cfg.at("opacity_hot")};
  RadiativeOpacity base(files),opacity(files,{data/"opacity/lifetime/hydrogen_response.dat",0,.16,1});
  const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
  const auto composition=[](double h,double y,double z){
    auto c=solar_scaled(h,z);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=y;c.X[2]-=y;return c;
  };
  // An ordinary source-node composition, two nearly pure-H layers and a
  // metal-rich core from the preserved late model exercise different branches.
  for(const auto& values:std::vector<std::array<double,5>>{
      {1.013e4,1e-4,.7,0,.02},{1.013e4,1e-4,1,0,0},
      {2.017e5,.01,.999,.0001,1e-9},{7e5,8.,.998,.001,1e-30},
      {6.443e6,3.982e4,3.528e-8,8.537e-15,.1312}}) {
    const auto [T,rho,h,y,z]=values;const auto c=composition(h,y,z);const auto s=opacity.eval(T,rho,c);
    require(std::isfinite(s.kappa) && s.kappa>0,"supported layer must have positive radiative opacity");
    constexpr double e=1e-6;
    const auto logk=[&](double t,double r){return std::log(opacity.eval(t,r,c).kappa);};
    require(std::abs((logk(T*std::exp(e),rho)-logk(T*std::exp(-e),rho))/(2*e)-s.dlnk_dlnT)<2e-5,
        "assembled temperature derivative disagrees with finite differences");
    require(std::abs((logk(T,rho*std::exp(e))-logk(T,rho*std::exp(-e)))/(2*e)-s.dlnk_dlnRho)<2e-5,
        "assembled density derivative disagrees with finite differences");
  }
  auto source_node=solar_scaled(.7,.02);source_node.metal_inventory=MetalInventory::gs98;
  source_node.basis=AbundanceBasis::atomic_mass;
  require(std::abs(opacity.eval(1.013e4,1e-4,source_node).kappa/
      base.eval(1.013e4,1e-4,source_node).kappa-1)<1e-13,
      "the actual atomic source-node opacity must be retained");
  const double T=7e5;const auto hrich=composition(.998,.001,1e-30);
  for(double logr:{1.25,1.45}) {
    const double rho=std::pow(10.,logr)*std::pow(T/1e6,3),e=1e-9;
    const auto a=opacity.eval(T,rho*std::exp(-e),hrich),b=opacity.eval(T,rho*std::exp(e),hrich);
    require(std::abs(std::log(a.kappa/b.kappa))<1e-7,"density join is not continuous");
    require(std::abs(a.dlnk_dlnRho-b.dlnk_dlnRho)<1e-5,"density join derivative is not continuous");
  }
  bool rejected=false;try{opacity.eval(T,100.,hrich);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"dense hydrogen approximation must reject outside its declared density range");
  rejected=false;try{opacity.eval(6e6,4e4,composition(.1,0,.17));}catch(const std::domain_error&){rejected=true;}
  require(rejected,"upper metallicity bound must be enforced");
  rejected=false;try{opacity.eval(10.,1e-4,hrich);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"composition extension must not extrapolate temperature");
  std::cout<<"PASS joined radiative opacity: physical layers, derivatives, continuous density joins and domain rejections\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}

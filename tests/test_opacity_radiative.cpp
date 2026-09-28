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
  RadiativeOpacity expanded(files,{data/"opacity/lifetime/hydrogen_response.dat",0,.16,1,2.2});
  const double newT=4.15e5,newrho=6.17,e=1e-6;
  const auto more=expanded.eval(newT,newrho,hrich);
  require(more.kappa>0,"expanded approximation must query retained source support");
  const auto lk=[&](double t,double r){return std::log(expanded.eval(t,r,hrich).kappa);};
  require(std::abs((lk(newT*std::exp(e),newrho)-lk(newT*std::exp(-e),newrho))/(2*e)-more.dlnk_dlnT)<2e-5,
      "expanded temperature derivative mismatch");
  require(std::abs((lk(newT,newrho*std::exp(e))-lk(newT,newrho*std::exp(-e)))/(2*e)-more.dlnk_dlnRho)<2e-5,
      "expanded density derivative mismatch");
  const auto newrange=expanded.density_range(newT,hrich);
  require(newrange && newrange->min<newrho && newrange->max>newrho,"expanded density range missing");
  RadiativeOpacity warmer(files,{data/"opacity/lifetime/hydrogen_response.dat",0,.16,1,2.2,6.3});
  const double hotT=1.4e6,hotrho=60.;
  const auto warm=warmer.eval(hotT,hotrho,hrich);
  const auto hotlog=[&](double t,double r){return std::log(warmer.eval(t,r,hrich).kappa);};
  require(std::abs((hotlog(hotT*std::exp(e),hotrho)-hotlog(hotT*std::exp(-e),hotrho))/(2*e)-warm.dlnk_dlnT)<2e-5,
      "warmer source-slope temperature derivative mismatch");
  require(std::abs((hotlog(hotT,hotrho*std::exp(e))-hotlog(hotT,hotrho*std::exp(-e)))/(2*e)-warm.dlnk_dlnRho)<2e-5,
      "warmer source-slope density derivative mismatch");
  const auto hotrange=warmer.density_range(hotT,hrich);
  require(hotrange && hotrange->min<hotrho && hotrange->max>hotrho,
      "warmer source-slope density range missing");
  bool temperature_rejected=false;
  try{expanded.eval(hotT,hotrho,hrich);}catch(const std::domain_error&){temperature_rejected=true;}
  require(temperature_rejected,"default temperature limit must remain unchanged");
  RadiativeOpacity hotter(files,{data/"opacity/lifetime/hydrogen_response.dat",0,.16,1,2.5,6.6});
  for(const auto& [t,r]:std::array<std::pair<double,double>,3>{{{2.1e6,200.},{3e6,600.},{3.95e6,1500.}}}) {
    const auto state=hotter.eval(t,r,hrich);
    const auto logk=[&](double tt,double rr){return std::log(hotter.eval(tt,rr,hrich).kappa);};
    require(std::isfinite(state.kappa) && state.kappa>0,"hot envelope approximation needs valid source values");
    require(std::abs((logk(t*std::exp(e),r)-logk(t*std::exp(-e),r))/(2*e)-state.dlnk_dlnT)<2e-5,
        "hot envelope temperature derivative mismatch");
    require(std::abs((logk(t,r*std::exp(e))-logk(t,r*std::exp(-e)))/(2*e)-state.dlnk_dlnRho)<2e-5,
        "hot envelope density derivative mismatch");
    const auto range=hotter.density_range(t,hrich);
    require(range && range->min<r && range->max>r,"hot envelope source range missing");
  }
  temperature_rejected=false;
  try{warmer.eval(2.1e6,200.,hrich);}catch(const std::domain_error&){temperature_rejected=true;}
  require(temperature_rejected,"previously selected upper temperature still rejects unsupported requests");
  RadiativeOpacity denser(files,{data/"opacity/lifetime/hydrogen_response.dat",0,.16,1,2.5,6.3});
  for(double logr:{2.21,2.35,2.49}) {
    const double t=4.15e5,r=std::pow(10.,logr)*std::pow(t/1e6,3);
    const auto state=denser.eval(t,r,hrich);
    const auto logk=[&](double tt,double rr){return std::log(denser.eval(tt,rr,hrich).kappa);};
    require(std::isfinite(state.kappa) && state.kappa>0,"denser approximation must remain in source support");
    require(std::abs((logk(t*std::exp(e),r)-logk(t*std::exp(-e),r))/(2*e)-state.dlnk_dlnT)<2e-5,
        "denser temperature derivative mismatch");
    require(std::abs((logk(t,r*std::exp(e))-logk(t,r*std::exp(-e)))/(2*e)-state.dlnk_dlnRho)<2e-5,
        "denser density derivative mismatch");
    const auto range=denser.density_range(t,hrich);
    require(range && range->min<r && range->max>r,"denser source density range missing");
    bool previous_rejected=false;
    try{warmer.eval(t,r,hrich);}catch(const std::domain_error&){previous_rejected=true;}
    require(previous_rejected,"existing selected density limit must remain unchanged");
  }
  bool denser_rejected=false;
  try{denser.eval(4.15e5, std::pow(10.,2.51)*std::pow(.415,3),hrich);}
  catch(const std::domain_error&){denser_rejected=true;}
  require(denser_rejected,"expanded approximation must retain its declared density limit");
  RadiativeOpacity cooling(files,{data/"opacity/lifetime/hydrogen_response.dat",0,.16,1,3.5,6.6});
  for(double logr:{2.51,2.9,3.49}) {
    const double t=4.05e5,r=std::pow(10.,logr)*std::pow(t/1e6,3);
    const auto state=cooling.eval(t,r,hrich);
    const auto logk=[&](double tt,double rr){return std::log(cooling.eval(tt,rr,hrich).kappa);};
    require(std::abs((logk(t*std::exp(e),r)-logk(t*std::exp(-e),r))/(2*e)-state.dlnk_dlnT)<2e-5,
        "cooling-envelope temperature derivative mismatch");
    require(std::abs((logk(t,r*std::exp(e))-logk(t,r*std::exp(-e)))/(2*e)-state.dlnk_dlnRho)<2e-5,
        "cooling-envelope density derivative mismatch");
    const auto range=cooling.density_range(t,hrich);
    require(range && range->min<r && range->max>r,"cooling-envelope source range missing");
    bool prior_rejected=false;
    try{denser.eval(t,r,hrich);}catch(const std::domain_error&){prior_rejected=true;}
    require(prior_rejected,"previous selected density bound must still reject denser states");
  }
  bool source_rejected=false;
  try{cooling.eval(3.89e6,1.9e4,hrich);}catch(const std::domain_error&){source_rejected=true;}
  require(source_rejected,"larger approximation domain must not bypass actual source support");
  bool old_rejected=false;try{opacity.eval(newT,newrho,hrich);}catch(const std::domain_error&){old_rejected=true;}
  require(old_rejected,"default approximation domain must remain unchanged");
  bool rejected=false;try{opacity.eval(T,100.,hrich);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"dense hydrogen approximation must reject outside its declared density range");
  rejected=false;try{opacity.eval(6e6,4e4,composition(.1,0,.17));}catch(const std::domain_error&){rejected=true;}
  require(rejected,"upper metallicity bound must be enforced");
  rejected=false;try{opacity.eval(10.,1e-4,hrich);}catch(const std::domain_error&){rejected=true;}
  require(rejected,"composition extension must not extrapolate temperature");
  std::cout<<"PASS joined radiative opacity: physical layers, derivatives, continuous density joins and domain rejections\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}

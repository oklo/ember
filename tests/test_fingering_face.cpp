#include "ember/fingering_transport.hpp"
#include "ember/convection.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
using namespace ember;
struct UnusedRadiation final:Opacity {
 OpacityState eval(double,double,const Composition&)const override{throw std::runtime_error("stable face evaluated unsupported radiative transport");}
 const char* name()const override{return "unused radiation";}
};
struct ConstantRadiation final:Opacity {
 OpacityState eval(double,double,const Composition&)const override{return {1e4,0,0};}
 const char* name()const override{return "constant test opacity";}
};
Composition material(double X,double Y3,double Y4){auto c=solar_scaled(X,0);c.X[1]=Y3;c.X[2]=Y4;c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;return c;}
int main(int argc,char**argv)try{
 const char* family=argc>1?argv[1]:std::getenv("EMBER_COLD_HELIUM_FAMILY");
 const char* collisions=argc>2?argv[2]:std::getenv("EMBER_COLLISION_TABLE");
 if(!family||!collisions){std::cout<<"requires an EOS family and collision table\n";return 77;}
 ColdHeliumOptions o;o.mixture_phase=true;o.liquid_continuation_gamma=200;o.dense_transition=true;
 VariableMetalHelmholtzEos eos(family,HelmholtzTableEos::Mixture::allow_documented_proxy,
  VariableMetalHelmholtzEos::LowMetalInterpolation::quadratic,{},true,o,VariableMetalHelmholtzEos::IsotopeInterpolation::number_density);
 ScreenedCollisionTransport col(collisions);ScreenedMetalMicroscopicTransport base(eos,col,true,5e4);
 UnusedRadiation radiation;BrownFingeringTransport transport(eos,radiation,base,col,OscillatoryMixing::growth_squared);
 // Near-pure-H faces. Separate density inversions can reverse a stable
 // contrast; rounding dominant H can also leave unresolved He4 gradients.
 // The third case reproduces an evolution failure at a trace-He boundary.
 const std::array<std::array<double,10>,3> cases{{
 {11.194442944388559,4.0140572793315208,0.99999999999995637,4.3446745120768268e-14,1.8501974700038406e-16,11.176610325107879,3.9787128997547798,0.99999999999995948,4.041346672672761e-14,1.0967367209060386e-16},
 {11.096837321990273,3.837836243491247,0.99999999999996747,3.2247387874146398e-14,2.8214674737068859e-16,11.074627321571418,3.8030760731722948,0.99999999999996847,3.1194851428973045e-14,3.3548247038140087e-16},
 {11.088509243733608,3.8385093208614696,0.99999999999999623,3.7452899641124248e-15,2.9468319613107479e-17,11.067293516777376,3.8037129752764316,0.99999999999999623,3.6023675090078456e-15,1.7239077471768668e-16}
 }};
 unsigned failures=0,queries=0;set_composition_buoyancy_reuse(1e-4);
 for(const auto& x:cases)for(int it=-5;it<=5;++it)for(int ir=-3;ir<=3;++ir){
  Point lo{std::log(1e9),x[1]+ir*7e-6,x[0]+it*3e-6,0},hi{std::log(1.01e9),x[6]+ir*7e-6,x[5]+it*3e-6,0};
  set_composition_buoyancy_reuse(1e-4);
  auto a=material(x[2],x[3],x[4]),b=material(x[7],x[8],x[9]);++queries;
  try{const auto f=transport.face(0,1e32,1.001e32,lo,a,hi,b);
    if(f.mass_conductance!=0||f.heat_luminosity!=0)throw std::runtime_error("stable face acquired mixing");
  }catch(const std::exception&e){if(++failures<4)std::cerr<<e.what()<<'\n';}
 }
 // A resolved destabilizing gradient outside the material range must still
 // refuse; this change supplies no cool-atmosphere mixing prescription.
 for(double amplitude:{1e-9,1e-7}) {
  const auto&x=cases[0];Point lo{std::log(1e9),x[1],x[0],0},hi{std::log(1.01e9),x[6],x[5],0};
  auto a=material(1-2*amplitude,amplitude,amplitude),b=material(1-3*amplitude,2*amplitude,amplitude);
  bool refused=false;try{transport.face(0,1e32,1.001e32,lo,a,hi,b);}catch(const std::domain_error&){refused=true;}
  if(!refused){++failures;std::cerr<<"resolved instability was suppressed\n";}
 }
 // A dense, ionized H-rich face below 120 kK remains an active liquid
 // instability. The 100 kK floor must still reject a colder inverse gradient.
 ConstantRadiation constant;BrownFingeringTransport cold(eos,constant,base,col,OscillatoryMixing::growth_squared);
 for(double scale:{.9,.8}) {
  Point lo{21.165603919388833,6.833333611894318,11.726500180576595+std::log(scale),6.1127900473684186e26};
  Point hi{21.170643797255362,6.8052882453780255,11.724535611766834+std::log(scale),6.1278682066601316e26};
  auto a=material(.99994980028111979,4.9735371667703596e-5,4.643472125025384e-7);
  auto b=material(.99994929786654618,5.0235425225314747e-5,4.6670822850238615e-7);
  bool refused=false;
  try{const auto f=cold.face(409,1.9223496009829335e32,1.9245518375458182e32,lo,a,hi,b);
    if(!(f.diffusivity>0 && f.mass_conductance>0 && f.heat_luminosity<0)){
      ++failures;std::cerr<<"cold inverse gradient has no mixing response\n";
    }
  }catch(const std::domain_error&){refused=true;}
  if(refused!=(scale<.9)){++failures;std::cerr<<"cold mixing material range is wrong\n";}
 }
 set_composition_buoyancy_reuse(0);
 std::cout<<queries<<" stable or unresolved faces, "<<failures<<" failures\n";return failures?1:0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}

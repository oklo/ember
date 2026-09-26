#include "ember/material_flux.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace ember;
namespace {
void close(double a,double b,double tolerance=2e-12) {
  if(!std::isfinite(a+b) || std::abs(a-b)>tolerance*std::max({std::abs(a),std::abs(b),1.}))
    throw std::runtime_error("material flux comparison failed");
}
}
int main() {
  try {
    MaterialMatrix reduced{{{{2,.4,.2}},{{.4,1,-.15}},{{.2,-.15,3}}}};
    const std::array<double,2> enthalpy{7,-3},a{.7,-.3},b{1.1,-.5};
    const auto full=full_material_mobility(reduced,enthalpy,5,10,4);
    const auto reference=material_face_flux(full,a,b,3,6,2,7,11,10,4);
    close(reference.material_luminosity,reference.carried_luminosity+reference.conductive_luminosity);
    if(reference.entropy_production<0)throw std::runtime_error("negative face entropy");
    const auto reverse=material_face_flux(full,b,a,6,3,2,7,11,10,4);
    close(reference.material_luminosity,-reverse.material_luminosity);
    close(reference.entropy_production,reverse.entropy_production);
    for(std::size_t i=0;i<2;++i)close(reference.species_rate[i],-reverse.species_rate[i]);
    // Unit changes must preserve the complete physical flux, including all
    // species/heat cross terms. This checks both sides of the transformation.
    for(double e0:{.01,1.,100.})for(double s0:{.01,1.,1e8}) {
      const auto matrix=full_material_mobility(reduced,enthalpy,5,e0,s0);
      const auto result=material_face_flux(matrix,a,b,3,6,2,7,11,e0,s0);
      for(std::size_t i=0;i<2;++i)close(reference.species_rate[i],result.species_rate[i]);
      close(reference.material_luminosity,result.material_luminosity);
      close(reference.carried_luminosity,result.carried_luminosity);
      close(reference.conductive_luminosity,result.conductive_luminosity);
      close(reference.entropy_production,result.entropy_production);
    }
    const auto isothermal=material_face_flux(full,a,b,3,3,2,7,11,10,4);
    close(isothermal.conductive_luminosity,0);close(isothermal.material_luminosity,isothermal.carried_luminosity);
    // Set the chemical force to cancel species motion. The remaining energy
    // flux must be pure conduction, even with nonzero reciprocal cross terms.
    const auto split=split_material_heat(full,std::sqrt(18.),10,4);
    std::array<double,2> phi{};
    for(std::size_t i=0;i<2;++i)phi[i]=-split.transported_enthalpy[i]*(1.0/3-1.0/6);
    const auto heat=material_face_flux(full,{},phi,3,6,2,7,11,10,4);
    close(heat.species_rate[0],0);close(heat.species_rate[1],0);
    close(heat.material_luminosity,heat.conductive_luminosity);
    std::cout<<"material face flux controls passed\n";return 0;
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

#include "ember/ion_classical.hpp"
#include "ember/detail/ion_free_energy.hpp"
#include "ember/detail/ion_mixture.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace ember {
namespace {
using J=detail::Taylor3<1>;
constexpr double charge=4.80320471257e-10,maximum_gamma=200;
HelmholtzJet per_mass(double logT,double logne,double A,double Z) {
  const double lg=std::log(charge*charge/constants::kB)
      +(std::log(4*std::numbers::pi/3)+logne)/3+(5./3)*std::log(Z)-logT;
  const auto g=detail::exp(J::variable(lg,0));J f;
  if(g.value()<=maximum_gamma)f=detail::classical_ocp_free_energy(g);
  else {
    const auto anchor=detail::classical_ocp_free_energy(
        detail::exp(J::variable(std::log(maximum_gamma),0)));
    f=anchor.value()+anchor.derivative({1})*(g/maximum_gamma-1);
  }
  HelmholtzJet out{};
  for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)
    out[i][j]=constants::R_gas/A*(i%2?-1:1)*std::pow(1./3,j)*f.derivative({i+j});
  return out;
}
}
std::array<HelmholtzJet,10> ion_classical_liquid_jets(
    double T,double rho,const Composition& c,std::size_t channels) {
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0)
      ||(channels!=1&&channels!=4&&channels!=10))
    throw std::domain_error("classical ion EOS: invalid state or derivative request");
  if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98
      ||c[Species::H2]!=0)
    throw std::domain_error("classical ion EOS: requires baryonic GS98 material after D mapping");
  for(double x:c.X)if(!std::isfinite(x)||x<0)
    throw std::domain_error("classical ion EOS: invalid abundance");
  if(std::abs(c.sum()-1)>1e-10)
    throw std::domain_error("classical ion EOS: mass fractions do not sum to one");
  return detail::common_density_ion_jets(T,rho,c,channels,per_mass);
}
}

#include "ember/detail/ion_mixture.hpp"
#include "ember/ion_quantum.hpp"
#include "ember/constants.hpp"
#include "ember/detail/taylor3.hpp"
#include <cmath>
#include <numbers>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ember {
namespace {
using J=detail::Taylor3<1>;
constexpr double pi=std::numbers::pi;
constexpr double charge=4.80320471257e-10;
constexpr double hbar=constants::h/(2*pi);
constexpr double electron_metal=gs98_ion_moment(1);

// Equation 34 is ln[sinh(x/2)/(x/2)]. Within the explicitly assessed
// theta_H<=1.5 range its series through x^14 has a small truncation error
// per oscillator. Density coefficients retain their derivatives. The
// temperature dependence of each term is exactly T^(-2m), avoiding a
// costly multivariable expansion at every stellar mesh point.
HelmholtzJet per_mass(double logT,double logne,double A,double Z) {
  const auto n=J::variable(logne,0);
  const auto ai=detail::exp((std::log(3*Z/(4*pi))-n)/3);
  const auto rs=ai*(A*constants::amu*Z*Z*charge*charge/(hbar*hbar));
  const auto c1=.351*rs/(90+rs);
  const auto b1=c1*c1,b2=J(.294*.294),b3=1-b1-b2;
  const auto theta2=detail::exp(n-2*logT)*
      (hbar*hbar/(constants::kB*constants::kB)*4*pi*Z*charge*charge/(A*constants::amu));
  constexpr std::array<double,7> coefficient{1./24,-1./2880,1./181440,
      -1./9676800,1./479001600,-691./15692092416000,1./1046139494400};
  J power(1),p1(1),p2(1),p3(1);HelmholtzJet out{};
  for(unsigned m=1;m<=coefficient.size();++m) {
    power=power*theta2;p1=p1*b1;p2=p2*b2;p3=p3*b3;
    const auto term=(constants::R_gas/A*coefficient[m-1])*power*(m==1?J(1):p1+p2+p3);
    double thermal=1;
    for(unsigned i=0;i<4;++i) {
      for(unsigned j=0;i+j<=3;++j)out[i][j]+=thermal*term.derivative({j});
      thermal*=-2.*m;
    }
  }
  return out;
}
}
std::array<HelmholtzJet,10> ion_quantum_liquid_jets(
    double T,double rho,const Composition& c,std::size_t channels) {
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0) ||
      (channels!=1&&channels!=4&&channels!=10))
    throw std::domain_error("quantum ion EOS: invalid state or derivative request");
  if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98 || c[Species::H2]!=0)
    throw std::domain_error("quantum ion EOS: requires baryonic GS98 material after D mapping");
  for(double x:c.X)if(!std::isfinite(x)||x<0)
    throw std::domain_error("quantum ion EOS: invalid abundance");
  if(std::abs(c.sum()-1)>1e-10)throw std::domain_error("quantum ion EOS: mass fractions do not sum to one");
  const double Z=c.Z(),Ye=c.X[0]+(2./3)*c.X[1]+.5*c.X[2]+electron_metal*Z;
  const double ne=rho*Ye/constants::amu;
  const double thetaH=hbar/(constants::kB*T)*std::sqrt(4*pi*ne*charge*charge/constants::amu);
  // Common-electron-density linear mixing reproduces the leading quantum
  // term. The assessed helium-rich core range also admits theta_H<=1.5;
  // its higher-order mixture dependence is an approximation, not a phase EOS.
  // FreeEOS with variable ionization supports the dense H-rich extension:
  // over this range the electron deficit is <2.8e-5, including pressure
  // ionization. Elsewhere below 300 kK the largest ideal-ion heat correction
  // must remain < theta_H^2/18 <= 5.556e-4. These checks limit the domain;
  // they do not switch or weight the free energy.
  const bool dense_hydrogen=T>=2e5 && rho>=50 && rho<=150 &&
      c.X[0]>=.97 && c.Z()<=1e-8 && c.X[1]<=.5*(c.X[1]+c.X[2]);
  // In the cooler H-rich layers, direct equilibrium-ionization comparisons
  // find a maximum sampled difference of 0.135% of classical ion Cv in this
  // small correction. On the 3200 K profiles, the mass-weighted uncertainty
  // estimate is <2.3e-6 of whole-star classical ion Cv. This extends the
  // assessed domain, without tapering the potential or altering its derivatives.
  // At 2488 K, the additional 21 envelope shells carry 0.03875% of the
  // resolved mass. Equilibrium-ionization derivatives change the integrated
  // heat capacity by less than 1e-8. Controls through theta_H=.7 differ
  // by at most 0.4033% of classical ion Cv and 2.253e-5 of source pressure.
  // At T>=200 kK, the assessed pressure-ionized H-rich range extends to
  // rho=500 g/cm3. Direct equilibrium-ionization controls change the quantum
  // Cv by less than 3.1e-10 of classical ion Cv, at the differencing floor.
  // The colder envelope keeps its narrower density range.
  const bool cool_hydrogen=T>=5e4 && rho>=10 && (rho<=200 || (T>=2e5 && rho<=500)) && thetaH<=.7 &&
      c.X[0]>=.98 && c.Z()<=1e-8 && c.X[1]<=.5*(c.X[1]+c.X[2]);
  const bool assessed_core=T>=5e5 && rho>=1e4 && rho<=1e5 && c.X[0]<=.05 && c.Z()<=.16;
  if(thetaH>(assessed_core?1.5:1.) || (T<3e5 && thetaH>.1 && !dense_hydrogen && !cool_hydrogen)) {
    std::ostringstream message;
    message << std::setprecision(4) << "quantum ion EOS: outside assessed ionization/quantum range"
      << " (T=" << T << ", rho=" << rho << ", X=" << c.X[0]
      << ", Y3=" << c.X[1] << ", Z=" << Z << ", theta_H=" << thetaH << ')';
    throw std::domain_error(message.str());
  }
  double ion_number=c.X[0]+c.X[1]/3+c.X[2]/4;
  double charge_moment=c.X[0]+(c.X[1]/3+c.X[2]/4)*std::pow(2.,5./3);
  for(const auto& m:gs98_metals) {
    ion_number+=Z*m.fraction/m.mass_number;
    charge_moment+=Z*m.fraction/m.mass_number*std::pow(m.charge,5./3);
  }
  const double ae=std::cbrt(3/(4*pi*ne));
  if(charge*charge/(ae*constants::kB*T)*charge_moment/ion_number>100)
    throw std::domain_error("quantum ion EOS: strong-coupling phase treatment required");
  return detail::common_density_ion_jets(T,rho,c,channels,per_mass);
}
} // namespace ember

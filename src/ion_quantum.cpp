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

// Equation 34 is log[sinh(x/2)/(x/2)]. Use the exact expression when
// quantum effects grow, and its series at small x to avoid cancellation.
// Thermal responses are analytic; the density expansion differentiates
// the density-dependent coefficients as well as the plasma frequency.
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
  HelmholtzJet out{};
  if(theta2.value()>.25) {
    for(const auto b:{b1,b2,b3}) {
      const auto x2=theta2*b;
      if(x2.value()<=.25) {
        J power(1);
        for(unsigned m=1;m<=coefficient.size();++m) {
          power=power*x2;
          const auto term=(constants::R_gas/A*coefficient[m-1])*power;
          double thermal=1;
          for(unsigned i=0;i<4;++i) {
            for(unsigned j=0;i+j<=3;++j)out[i][j]+=thermal*term.derivative({j});
            thermal*=-2.*m;
          }
        }
      } else {
        const auto x=detail::sqrt(x2),e=detail::exp(-x),q=1-e;
        const auto f=detail::log(q)+x/2-detail::log(x);
        const auto u=x*e/q+x/2-1;
        const auto a=x2*e/(q*q);
        const auto d2=u+1-a;
        const std::array<J,4> thermal{f,-u,d2,-d2-2*a*u};
        for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)
          out[i][j]+=constants::R_gas/A*thermal[i].derivative({j});
      }
    }
    return out;
  }
  J power(1),p1(1),p2(1),p3(1);
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
HelmholtzJet detail::bc22_liquid_quantum_per_mass(double logT,double logne,double A,double Z) {
  return per_mass(logT,logne,A,Z);
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
  // Common-density linear mixing has the correct leading quantum limit.
  // Each ion uses its own mass, so He3 is included at any share of the
  // helium: its quantum ratio is sqrt(2/3) of the hydrogen ratio bounded below.
  // Larger corrections and partial ionization require the bounded core or
  // hydrogen domains below. These guards never taper the free energy.
  // Physical comparisons and remaining mixture limits are in docs/DENSE_EOS.md.
  const bool dense_hydrogen=T>=2e5 && rho>=50 && rho<=150 &&
      c.X[0]>=.97 && c.Z()<=1e-8;
  const bool cool_hydrogen=T>=5e4 && rho>=10 && (rho<=200 || (T>=2e5 && rho<=500)) && thetaH<=.7 &&
      c.X[0]>=.98 && c.Z()<=1e-8;
  // Dense hydrogen envelope at X >= 0.95: hydrogen fully ionized in direct
  // source checks; trace-helium recombination below ~130 kK is a declared term.
  const bool cold_hydrogen=T>=3.5e4 && rho>=10 && rho<=1000 && thetaH<=1.7 &&
      c.X[0]>=.95 && c.Z()<=1e-8;
  // More strongly quantum, nearly pure hydrogen. Direct equilibrium-ionization
  // controls assess the trace-He contribution separately; the classical EOS
  // retains its ionization physics. This only extends the admitted domain.
  const bool nearly_pure_hydrogen=T>=7e4 && rho>=10 && rho<=1000 && thetaH<=2.5 &&
      c.X[0]>=.999 && c.Z()<=1e-8;
  // Dense, fully ionized H/He layers of any hydrogen fraction (the He-rich
  // mantle and the H/He transition) from 160 to 300 kK.
  const bool dense_mixture=T>=1.6e5 && rho>=300 && rho<=4000 && thetaH<=1.85 &&
      c.Z()<=1e-10;
  // The helium-dominated cold overlap still requires its warm source and
  // anchors. Direct source checks retain full ionization in this interval.
  const bool dense_helium_overlap=T>=1.6e5 && rho>=300 && rho<=6000
      && c.X[0]<=.20 && c.Z()<=1e-8 && thetaH<=2.5;
  // The liquid mixture does not determine metal separation or freezing;
  // the mean-coupling guard below is not a mixture phase boundary.
  const bool assessed_core=T>=5e5 && rho>=1e3 && rho<=1e5 && c.X[0]<=.05 && c.Z()<=.16;
  // Cooler, nearly pure helium layers remain pressure ionized in direct
  // source checks. This does not extend the metal-rich core's phase range.
  const bool cool_helium=T>=4e5 && rho>=1e3 && rho<=1e5 && c.X[0]<=.01 && c.Z()<=1e-3;
  const double theta_limit=(assessed_core || cool_helium)?4.:dense_helium_overlap || nearly_pure_hydrogen?2.5:dense_mixture?1.85:cold_hydrogen?1.7:1.;
  if(thetaH>theta_limit ||
      (T<3e5 && thetaH>.1 && !dense_hydrogen && !cool_hydrogen && !cold_hydrogen
          && !dense_mixture && !dense_helium_overlap && !nearly_pure_hydrogen)) {
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

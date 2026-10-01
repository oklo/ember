#include "ember/dense_nuclear.hpp"
#include "ember/constants.hpp"
#include "ember/detail/differential.hpp"
#include <cmath>
#include <atomic>
#include <stdexcept>

namespace ember {
namespace { std::atomic<DenseNuclearModel> selected_model{DenseNuclearModel::none}; }
void set_dense_nuclear_model(DenseNuclearModel model) {
  if(model!=DenseNuclearModel::none&&model!=DenseNuclearModel::uniform_optimal&&model!=DenseNuclearModel::uniform_high)
    throw std::invalid_argument("unknown dense nuclear model");
  selected_model.store(model,std::memory_order_relaxed);
}
DenseNuclearModel dense_nuclear_model() { return selected_model.load(std::memory_order_relaxed); }

NuclearRateResponse dense_nuclear_rate(double T,double rho,const Composition& c,
    const NuclearPair& pair,DenseNuclearModel model) {
  using namespace constants;
  using D=detail::Differential<NSPEC+2>;
  using detail::exp;using detail::log;using detail::sqrt;using detail::cbrt;
  using detail::pow;
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0))
    throw std::domain_error("dense nuclear rate: positive finite T and density required");
  for(double v:c.X)if(!std::isfinite(v)||v<0)
    throw std::domain_error("dense nuclear rate: invalid composition");
  for(double v:{pair.charge1,pair.charge2,pair.mass1,pair.mass2,pair.s0})
    if(!std::isfinite(v)||v<=0)throw std::domain_error("dense nuclear rate: invalid reaction");
  if(model!=DenseNuclearModel::uniform_optimal&&model!=DenseNuclearModel::uniform_high)
    throw std::invalid_argument("dense nuclear rate: explicit uniform-mixture model required");
  const bool high=model==DenseNuclearModel::uniform_high;
  // Printed MCP entries of Table II. Equation 38 differs for the non-optimal
  // C_T entry; the parameter-table choice is retained explicitly here.
  const double Cexp=high?2.450:2.638,Cpyc=high?50.:3.90,Cpl=1.25;
  const double CT=high?.840:.724,al=high?1.05:1.,aw=high?.95:1.,Lambda=high?.35:.5;
  constexpr double e2=4.803204673e-10*4.803204673e-10,mev=1.602176634e-6;
  const double hb=h/(2*M_PI),Zi=pair.charge1,Zj=pair.charge2,Ai=pair.mass1,Aj=pair.mass2;
  const double mu=amu*Ai*Aj/(Ai+Aj),rB=hb*hb/(2*mu*Zi*Zj*e2);
  const double Ea=2*mu*Zi*Zi*Zj*Zj*e2*e2/(hb*hb);
  const double zsum=std::cbrt(Zi)+std::cbrt(Zj);
  const double Csc=.9*(std::pow(Zi+Zj,5./3)-std::pow(Zi,5./3)-std::pow(Zj,5./3))*zsum/(2*Zi*Zj);
  const D temp=exp(D::variable(std::log(T),0)),density=exp(D::variable(std::log(rho),1));
  D ions,electrons,mass;
  for(std::size_t j=0;j<NSPEC;++j) {
    const auto x=D::variable(c.X[j],j+2);mass+=x;
    if(c.metal_inventory==MetalInventory::gs98&&is_metal_species(j)) {
      ions+=x*c.metal_ion_moment(0);electrons+=x*c.metal_ion_moment(1);
    } else {
      ions+=x/c.abundance_weight(j);electrons+=x*nuclides[j].Z/c.abundance_weight(j);
    }
  }
  if(!(ions.value>0&&electrons.value>0))throw std::domain_error("dense nuclear rate: empty mixture");
  const D ne=density*NA*electrons,n=density*NA*ions;
  const D aij=.5*zsum*cbrt(3/(4*M_PI*ne));
  const D nij=3/(4*M_PI*aij*aij*aij),coupling=Zi*Zj*e2/(aij*kB*temp);
  const D Tp=hb/kB*sqrt(4*M_PI*Zi*Zj*e2*nij/(2*mu)),Tpw=aw*Tp;
  const D lam=rB*cbrt(nij/2),lt=al*lam;
  const D tau=3*std::pow(M_PI/2,2./3)*cbrt(Ea/(kB*temp)),zeta=3*coupling/tau;
  const D Tt=sqrt(temp*temp+CT*CT*Tp*Tp);
  const D taut=3*std::pow(M_PI/2,2./3)*cbrt(Ea/(kB*Tt)),Gt=Zi*Zj*e2/(aij*kB*Tt);
  const D phi=sqrt(coupling)/pow(pow(Csc/zeta,4)+coupling*coupling,.25);
  const D gamma=(temp*temp*(2./3)+Tpw*Tpw*(2./3)*(Cpl+.5))/(temp*temp+Tpw*Tpw);
  const D exponent=-taut+Csc*Gt*phi*exp(-Lambda*Tpw/temp)-Lambda*Tpw/temp;
  // Eq.36 divided by n_i n_j/(1+delta_ij), then multiplied by N_A.
  const double pref=8*std::cbrt(M_PI)/(std::sqrt(3)*std::cbrt(2.));
  const D thermal=exp(log(D(NA*pair.s0*mev*1e-24/hb*rB*pref))
                       +gamma*log(Ea/(kB*Tt))+exponent);
  // Eq.33 in its published numerical units. X_N=sum X and <A>=sum X/sum(X/A)
  // also define consistent unconstrained abundance derivatives off sum X=1.
  const D pyc=exp(log(1e46*Cpyc*8*density*mass*Ai*Aj*(mass/ions)*Zi*Zi*Zj*Zj
                           /((Ai+Aj)*(Ai+Aj))*pair.s0*NA/(n*n))
                   +(3-Cpl)*log(lt)-Cexp/sqrt(lt));
  const D combined=thermal+pyc;
  if(!std::isfinite(combined.value))throw std::domain_error("dense nuclear rate: nonfinite coefficient");
  NuclearRateResponse out;out.molar_rate=combined.value;
  if(combined.value>0) {
    out.dlnrate_dlnT=combined.d[0]/combined.value;
    out.dlnrate_dlnRho=combined.d[1]/combined.value;
    for(std::size_t j=0;j<NSPEC;++j)out.dlnrate_dX[j]=combined.d[j+2]/combined.value;
  }
  return out;
}
} // namespace ember

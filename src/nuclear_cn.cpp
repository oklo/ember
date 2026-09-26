#include "ember/nuclear_cn.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
constexpr std::array<CNReaction,3> reactions{CNReaction::c12_p,CNReaction::c13_p,CNReaction::n14_p};
// Moles produced per mole of the corresponding reaction, in Species order.
constexpr std::array<std::array<double,NSPEC>,3> stoich{{
  {-1,0,0,-1,1,0,0,0}, {-1,0,0,0,-1,1,0,0}, {-2,0,1,1,0,-1,0,0}
}};
void check_lookup(const Composition& c) {
  if(c.basis!=AbundanceBasis::baryon_mass)
    throw std::invalid_argument("CNNetwork: baryon-mass lookup required");
  for(double x:c.X)if(!std::isfinite(x) || x<0)
    throw std::domain_error("CNNetwork: nonnegative finite lookup fractions required");
}
void check_catalysts(const CNAbundances& y) {
  for(double v:y)if(!std::isfinite(v) || v<0)
    throw std::domain_error("CNNetwork: nonnegative finite catalyst molalities required");
}
}
CNAbundances initial_gs98_cn(const Composition& c) {
  check_lookup(c);
  if(c.metal_inventory!=MetalInventory::gs98)
    throw std::invalid_argument("CNNetwork: initial GS98 inventory required");
  return {c.Z()*gs98_metals[0].fraction/12,0,c.Z()*gs98_metals[1].fraction/14};
}
CNNetwork::CNNetwork(PPRates rates,PPScreening screening):rates_(rates),screening_(screening) {
  if(rates!=PPRates::solar_fusion_ii && rates!=PPRates::solar_fusion_iii)
    throw std::invalid_argument("CNNetwork: Solar Fusion II or III required");
  if(screening!=PPScreening::debye_fermi && screening!=PPScreening::salpeter_van_horn)
    throw std::invalid_argument("CNNetwork: finite-degeneracy screening required");
}
CNNetworkResponse CNNetwork::response(double T,double rho,const Composition& c,const CNAbundances& y) const {
  using namespace constants;
  check_lookup(c);check_catalysts(y);
  if(!std::isfinite(T) || !std::isfinite(rho) || T<=0 || rho<=0 || T>2e7)
    throw std::domain_error("CNNetwork: finite positive T<=20 MK and density required");
  CNNetworkResponse out;
  if(T<1e5)return out;
  constexpr double mev=1.602176634e-6;
  constexpr std::array<double,3> neutrino{.706*mev*NA,0,.996*mev*NA};
  double epsT=0,epsR=0;
  for(std::size_t k=0;k<3;++k) {
    const auto bare=cn_bare_rate(T,reactions[k],rates_);
    const auto screen=cn_screening(T,rho,c,reactions[k],screening_);
    const double coefficient=rho*bare.molar_rate*std::exp(screen.log_factor);
    const double frequency=coefficient*c.X[0],rate=frequency*y[k];
    out.frequency[k]=frequency;out.reaction_rate[k]=rate;
    double q=0;
    for(std::size_t i=0;i<NSPEC;++i)
      q-=stoich[k][i]*(nuclides[i].A-mass_numbers[i])*c_light*c_light;
    const double heat=q-neutrino[k];
    auto& n=out.physical;auto& s=n.state;
    s.eps+=rate*heat;s.eps_neutrino+=rate*neutrino[k];
    epsT+=rate*heat*(bare.dlnrate_dlnT+screen.dlog_dlnT);
    epsR+=rate*heat*(1+screen.dlog_dlnRho);
    out.deps_dY[k]=frequency*heat;
    for(std::size_t i=0;i<NSPEC;++i) {
      s.dXdt[i]+=stoich[k][i]*mass_numbers[i]*rate;
      out.d_dXdt_dY[i][k]=stoich[k][i]*mass_numbers[i]*frequency;
    }
    for(std::size_t j=0;j<NSPEC;++j) {
      const double dr=rate*screen.dlog_dX[j]+(j==0?coefficient*y[k]:0);
      n.deps_dX[j]+=heat*dr;
      for(std::size_t i=0;i<NSPEC;++i)
        n.d_dXdt_dX[i][j]+=stoich[k][i]*mass_numbers[i]*dr;
    }
  }
  if(out.physical.state.eps>0) {
    out.physical.state.dlneps_dlnT=epsT/out.physical.state.eps;
    out.physical.state.dlneps_dlnRho=epsR/out.physical.state.eps;
  }
  return out;
}
NuclearResponse PPCNNetwork::composition_response(double T,double rho,const Composition& c) const {
  return response(T,rho,c,true);
}
NuclearResponse PPCNNetwork::response(double T,double rho,const Composition& c,bool initial_deuterium) const {
  if(!c.cn_molality || c.metal_inventory!=MetalInventory::gs98)
    throw std::invalid_argument("PPCNNetwork: explicit GS98 CN inventory required");
  // Unconstrained composition partials also evaluate unnormalized X, so the
  // stricter physical-ledger validation belongs at evolution boundaries.
  auto out=initial_deuterium?pp_.composition_response(T,rho,c):pp_.pp().composition_response(T,rho,c);
  const auto n=cn_.response(T,rho,c,*c.cn_molality).physical;
  const double t=out.state.eps*out.state.dlneps_dlnT+n.state.eps*n.state.dlneps_dlnT;
  const double r=out.state.eps*out.state.dlneps_dlnRho+n.state.eps*n.state.dlneps_dlnRho;
  out.state.eps+=n.state.eps;out.state.eps_neutrino+=n.state.eps_neutrino;
  out.state.dlneps_dlnT=out.state.eps>0?t/out.state.eps:0;
  out.state.dlneps_dlnRho=out.state.eps>0?r/out.state.eps:0;
  if(c.cn_mass_convention==CNMassConvention::explicit_metal_mass) {
    // Capture adds real metal mass; the N14 cycle returns it to helium.
    // The material table uses the actual total Z with its declared pattern.
    const double Z=c.Z(),metal_rate=n.state.dXdt[3]+n.state.dXdt[4]+n.state.dXdt[5];
    if(Z==0 && std::any_of(c.cn_molality->begin(),c.cn_molality->end(),[](double y){return y!=0;}))
      throw std::domain_error("PPCNNetwork: catalysts exceed zero metal inventory");
    for(std::size_t i=0;i<3;++i) {
      out.state.dXdt[i]+=n.state.dXdt[i];
      for(std::size_t j=0;j<NSPEC;++j)out.d_dXdt_dX[i][j]+=n.d_dXdt_dX[i][j];
    }
    for(std::size_t j=0;j<NSPEC;++j)out.deps_dX[j]+=n.deps_dX[j];
    if(Z>0)for(std::size_t i=3;i<METAL_END;++i) {
      const double fraction=c.X[i]/Z;out.state.dXdt[i]+=fraction*metal_rate;
      for(std::size_t j=0;j<NSPEC;++j) {
        const double derivative=n.d_dXdt_dX[3][j]+n.d_dXdt_dX[4][j]+n.d_dXdt_dX[5][j];
        out.d_dXdt_dX[i][j]+=fraction*derivative;
        if(is_metal_species(j))out.d_dXdt_dX[i][j]+=((i==j?1.:0.)-fraction)/Z*metal_rate;
      }
    }
    return out;
  }
  out.state.dXdt[0]+=n.state.dXdt[0];
  out.state.dXdt[2]-=n.state.dXdt[0];
  for(std::size_t j=0;j<NSPEC;++j) {
    out.deps_dX[j]+=n.deps_dX[j];
    out.d_dXdt_dX[0][j]+=n.d_dXdt_dX[0][j];
    out.d_dXdt_dX[2][j]-=n.d_dXdt_dX[0][j];
  }
  return out;
}
NuclearState PPCNNetwork::eval(double T,double rho,const Composition& c) const {
  // No primordial D is regenerated by the reduced pp network. Do not query
  // its screening domain or quadrature for the zero-reservoir thermal source.
  return response(T,rho,c,c[Species::H2]!=0).state;
}
double PPCNNetwork::rest_energy_correction(const Composition& c) const {
  if(!c.cn_molality)throw std::invalid_argument("PPCNNetwork: missing physical CN inventory");
  return cn_physical_ledger(c,*c.cn_molality).rest_energy_difference;
}
CNPhysicalLedger cn_physical_ledger(const Composition& lookup,const CNAbundances& y) {
  const auto initial=initial_gs98_cn(lookup);check_catalysts(y);
  if(std::abs(lookup.sum()-1)>1e-10)
    throw std::domain_error("CN ledger: normalized lookup required");
  if(lookup.cn_mass_convention==CNMassConvention::explicit_metal_mass) {
    CNPhysicalLedger out;out.molality=y;out.hydrogen=lookup.X[0];out.helium3=lookup.X[1];
    out.helium4=lookup.X[2];out.metal_fraction=lookup.Z();
    const double cn_mass=12*y[0]+13*y[1]+14*y[2],inert=out.metal_fraction-cn_mass;
    if(inert<0)throw std::domain_error("CN ledger: catalysts exceed total metal mass");
    double fraction=0,ions=0,electrons=0,rest=0;
    for(std::size_t i=2;i<gs98_metals.size();++i) {
      const auto& m=gs98_metals[i];fraction+=m.fraction;
      ions+=m.fraction/m.mass_number;electrons+=m.fraction*m.charge/m.mass_number;
      rest+=m.fraction*(m.atomic_weight/m.mass_number-1);
    }
    ions*=inert/fraction;electrons*=inert/fraction;rest*=inert/fraction;
    for(std::size_t k=0;k<3;++k) {
      ions+=y[k];electrons+=nuclides[k+3].Z*y[k];
      rest+=(nuclides[k+3].A-mass_numbers[k+3])*y[k];
    }
    out.ion_molality_difference=ions-out.metal_fraction*gs98_ion_moment(0);
    out.electron_molality_difference=electrons-out.metal_fraction*gs98_ion_moment(1);
    // Evolution's generic audit counts the carried X slots, whose proxy
    // isotope labels are not the actual GS98 rest masses. Replace that part
    // by the physical CN rest masses and a conserved inert-metal reference.
    for(std::size_t i=3;i<METAL_END;++i)rest-=(nuclides[i].A/mass_numbers[i]-1)*lookup.X[i];
    out.rest_energy_difference=rest*constants::c_light*constants::c_light;
    return out;
  }
  // Catalyst number is conserved by nuclear reactions, but microscopic
  // transport may change its LOCAL value. Keep inert metals fixed and
  // close the physical mass with He4; conservation is checked by the solver.
  CNPhysicalLedger out;out.molality=y;
  double ion=0,electron=0,rest=0;
  for(std::size_t k=0;k<3;++k) {
    const double delta=y[k]-initial[k];const auto j=k+3;
    out.extra_metal_mass+=mass_numbers[j]*delta;
    ion+=delta;electron+=nuclides[j].Z*delta;
    rest+=(nuclides[j].A-mass_numbers[j])*delta;
  }
  out.hydrogen=lookup.X[0];out.helium3=lookup.X[1];
  out.helium4=lookup.X[2]-out.extra_metal_mass;
  out.metal_fraction=lookup.Z()+out.extra_metal_mass;
  if(out.helium4<0 || out.metal_fraction<0)
    throw std::domain_error("CN ledger: negative physical helium or metal fraction");
  out.ion_molality_difference=ion-out.extra_metal_mass/4;
  out.electron_molality_difference=electron-out.extra_metal_mass/2;
  out.rest_energy_difference=(rest-(nuclides[2].A/4-1)*out.extra_metal_mass)
    *constants::c_light*constants::c_light;
  return out;
}
Composition explicit_cn_material(const Composition& source) {
  if(!source.cn_molality)throw std::invalid_argument("explicit CN material: missing catalyst inventory");
  const auto physical=cn_physical_ledger(source,*source.cn_molality);
  if(source.cn_mass_convention==CNMassConvention::explicit_metal_mass)return source;
  auto result=source;const double Z=source.Z();
  if(Z==0 && physical.metal_fraction!=0)throw std::domain_error("explicit CN material: missing material metal pattern");
  result.X[2]=physical.helium4;
  if(Z>0)for(std::size_t i=3;i<METAL_END;++i)result.X[i]*=physical.metal_fraction/Z;
  result.cn_mass_convention=CNMassConvention::explicit_metal_mass;
  (void)cn_physical_ledger(result,*result.cn_molality);return result;
}
CNAbundances cn_backward_euler(const CNAbundances& old,const std::array<double,3>& rate,double dt) {
  check_catalysts(old);
  if(!std::isfinite(dt) || dt<0)throw std::domain_error("CN implicit update: invalid timestep");
  std::array<double,3> v{};
  for(std::size_t k=0;k<3;++k) {
    if(!std::isfinite(rate[k]) || rate[k]<0 || !std::isfinite(v[k]=dt*rate[k]))
      throw std::domain_error("CN implicit update: invalid or unrepresentable capture frequency");
  }
  // Positive adjugate, with the cubic term cancelled symbolically. Scaling
  // avoids cancellation and overflow even far beyond catalyst equilibrium.
  const double scale=std::max({1.,v[0],v[1],v[2]});
  const double one=1/scale,a=v[0]/scale,b=v[1]/scale,c=v[2]/scale;
  const double d=one*one+one*(a+b+c)+a*b+a*c+b*c;
  if(!(d>0))throw std::domain_error("CN implicit update: unrepresentable frequency ratio");
  return {((one+b)*(one+c)*old[0]+b*c*old[1]+c*(one+b)*old[2])/d,
    (a*(one+c)*old[0]+(one+a)*(one+c)*old[1]+a*c*old[2])/d,
    (a*b*old[0]+b*(one+a)*old[1]+(one+a)*(one+b)*old[2])/d};
}
} // namespace ember

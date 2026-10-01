#include "ember/fingering_transport.hpp"
#include "ember/fingering.hpp"
#include "ember/viscosity.hpp"
#include "ember/convection.hpp"
#include "ember/eos_component.hpp"
#include "ember/constants.hpp"
#include "thermal_transport.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
constexpr double pi=std::numbers::pi;
MetalSpeciesMatrix inverse(MetalSpeciesMatrix a) {
  MetalSpeciesMatrix b{};for(std::size_t i=0;i<3;++i)b[i][i]=1;
  for(std::size_t j=0;j<3;++j) {
    std::size_t pivot=j;for(std::size_t i=j+1;i<3;++i)if(std::abs(a[i][j])>std::abs(a[pivot][j]))pivot=i;
    if(a[pivot][j]==0)throw std::domain_error("fingering: singular diffusion response");
    std::swap(a[j],a[pivot]);std::swap(b[j],b[pivot]);
    const double divisor=a[j][j];
    for(std::size_t k=0;k<3;++k){a[j][k]/=divisor;b[j][k]/=divisor;}
    for(std::size_t i=0;i<3;++i)if(i!=j) {
      const double factor=a[i][j];
      for(std::size_t k=0;k<3;++k){a[i][k]-=factor*a[j][k];b[i][k]-=factor*b[j][k];}
    }
  }
  return b;
}
MetalCNVector mixing_rates(const BrownFingeringTransport::Face& f,
    const Composition& a,const Composition& b) {
  const auto left=metal_cn_abundances(a),right=metal_cn_abundances(b);
  const MetalSpeciesVector contrast{a.X[0]-b.X[0],a.X[1]-b.X[1],a.Z()-b.Z()};
  MetalSpeciesFaceResponse cross;
  for(std::size_t row=0;row<3;++row)for(std::size_t col=0;col<3;++col)
    if(row!=2 || col!=2)cross.rate[row]+=f.mixing_conductance[row][col]*contrast[col];
  // Direct turbulent diffusion mixes each CN isotope, even at zero total-Z
  // gradient. Cross terms advect the metal group with its existing upwind rule.
  auto result=common_metal_cn_flux(cross,a,b,false).rate;
  for(std::size_t j=2;j<METAL_CN_SIZE;++j)
    result[j]+=f.mixing_conductance[2][2]*(left[j]-right[j]);
  return result;
}
}

BrownFingeringTransport::BrownFingeringTransport(const VariableMetalHelmholtzEos& eos,
    const Opacity& radiation,const MetalMicroscopicTransport& base,
    const ScreenedCollisionTransport& collisions,OscillatoryMixing oscillatory)
    :eos_(eos),radiation_(radiation),base_(base),collisions_(collisions),oscillatory_(oscillatory) {
  if(radiation.includes_conduction())
    throw std::invalid_argument("fingering: opacity must contain radiation only");
}

BrownFingeringTransport::Face BrownFingeringTransport::face(std::size_t index,
    double mlo,double mhi,const Point& lo,const Composition& a,const Point& hi,
    const Composition& b) const {
  if(a.X==b.X)return {};
  const double Tlo=std::exp(lo.lnT),Thi=std::exp(hi.lnT);
  const double rlo=std::exp(lo.lnrho),rhi=std::exp(hi.lnrho);
  const double T=.5*(Tlo+Thi),rho=.5*(rlo+rhi),mass=.5*(mlo+mhi);
  const auto e0=eos_.eval(Tlo,rlo,a),e1=eos_.eval(Thi,rhi,b);
  const double P=.5*(e0.P+e1.P),cp=.5*(e0.cp+e1.cp),ad=.5*(e0.grad_ad+e1.grad_ad);
  const double dp=std::log(e1.P/e0.P),delta=.5*(e0.delta+e1.delta);
  const double B=composition_buoyancy(eos_,T,P,delta,dp,a,b,rho).B;
  Face result;
  if(dp==0 || !(cp>0 && mass>0 && mhi>mlo))
    throw std::domain_error("fingering: invalid face geometry or thermodynamics");
  const double grad=(hi.lnT-lo.lnT)/dp;
  if(grad>=ad || grad>=ad+B)return result; // ordinary convection
  result.density_ratio=B==0?std::numeric_limits<double>::infinity():(ad-grad)/(-B);

  const auto c=mean_composition(a,b);
  if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98
      || c[Species::H2]!=0 || c.Z()>1e-6 || !(c.X[0]>0 && c.X[1]>0 && c.X[2]>0)
      || rho<300 || T<1e5) {
    if(B<0) {
      // Tiny composition contrasts can be smaller than the error in two
      // separate density inversions. Resolve their sign from the EOS response
      // at fixed pressure; a negative result still requires supported mixing.
      const std::array<double,3> dx{b.X[0]-a.X[0],b.X[1]-a.X[1],b.Z()-a.Z()};
      if(c.basis==AbundanceBasis::baryon_mass && c.metal_inventory==MetalInventory::gs98
          && c[Species::H2]==0 && c.Z()<=1e-6 && c.X[0]>.999
          && std::max({std::abs(dx[0]),std::abs(dx[1]),std::abs(dx[2])})
             <std::sqrt(std::numeric_limits<double>::epsilon())) {
        const double rp=eos_.rho_from_PT(T,P,c,rho);
        const auto response=eos_.isobaric_composition_response(T,rp,c,{c.X[0]>0,c.X[1]>0,c.Z()>0});
        double contrast=0;
        for(unsigned j=0;j<3;++j)if(dx[j]!=0)
          contrast=std::fma(response.dlnRho[j],dx[j],contrast);
        const double resolved=contrast/(delta*dp);
        // He4 is the dependent fraction. Near pure H, rounding the
        // dominant H abundance leaves an uncertainty of order epsilon in
        // the residual He4 and therefore in its density contrast. Do not
        // request an unsupported mixing law for an unresolved sign. Require
        // the thermal restoring contrast to exceed that bound by 64 times.
        const double density_roundoff=std::numeric_limits<double>::epsilon()
          *(std::abs(a.X[0])+std::abs(b.X[0]))*std::abs(response.dlnRho[0]);
        const bool unresolved=std::abs(contrast)<=density_roundoff
          && -delta*dp*(ad-grad)>64*density_roundoff;
        if(resolved>=0 || unresolved) {
          result.density_ratio=resolved==0 || unresolved?std::numeric_limits<double>::infinity():(ad-grad)/(-resolved);
          return result;
        }
      }
      throw std::domain_error("fingering coefficients require cold dense baryonic H/He liquid with negligible metals");
    }
    return result; // no multicomponent prescription is supplied outside this material domain
  }
  const double rp=eos_.rho_from_PT(T,P,c,rho);
  const auto p=eos_.isobaric_composition_response(T,rp,c,{c.X[0]>0,c.X[1]>0,c.Z()>0});
  const std::array<double,3> dx{b.X[0]-a.X[0],b.X[1]-a.X[1],b.Z()-a.Z()};
  const double thermal_buoyancy=-delta*dp*(ad-grad);
  if(!(thermal_buoyancy>0))throw std::domain_error("fingering: pressure must decrease outwards");
  if(p.dlnRho[0]*dx[0]<=0 && p.dlnRho[1]*dx[1]<=0 && B>=0)return result;
  const double hhe_drive=std::abs(p.dlnRho[0]*dx[0])+std::abs(p.dlnRho[1]*dx[1]);
  if(c.Z()>0 && std::abs(p.dlnRho[2]*dx[2])>1e-3*hhe_drive)
    throw std::domain_error("fingering: metal buoyancy needs an additional active field");
  const double ne=rho*c.mu_elec_inv()/constants::amu;
  // Add electron-ion collision frequencies at the common electron density.
  // The OCP ion viscosity of the average ion is a small, explicitly classical
  // mixture approximation. It is assessed separately from electron momentum.
  double frequency=0,ee=0,numerator=0;
  for(std::size_t j=0;j<3;++j)if(c.X[j]>0) {
    const double Z=nuclides[j].Z,A=mass_numbers[j];
    const auto v=ocp_liquid_viscosity(T,ne*A*constants::amu/Z,Z,A).value;
    frequency+=c.X[j]*Z/A/c.mu_elec_inv()*v.electron_ion_frequency;
    ee=v.electron_electron_frequency;
    numerator=v.electron*(v.electron_ion_frequency+ee);
  }
  const double Zi=c.mu_elec_inv()/c.mu_ions_inv(),Ai=1/c.mu_ions_inv();
  const auto average=ocp_liquid_viscosity(T,rho,Zi,Ai).value;
  const double viscosity=(numerator/(frequency+ee)+average.ion_classical)/rho;

  const auto electron=ElectronGas{}.eval_with_derivatives(T,rho,c);
  const double stiffness=electron.dP_dlnRho/(c.mu_elec_inv()*constants::NA*rho);
  const double screening=ScreenedCollisionTransport::screening_length(T,rho,c.X[0],c.X[1],c.Z(),stiffness,true);
  const auto kinetic=collisions_.bulk_metal_eval(T,rho,c.X[0],c.X[1],c.Z(),screening);
  const std::array<bool,3> active{c.X[0]>0,c.X[1]>0,c.Z()>0};
  MetalSpeciesMatrix chemical{};
  for(std::size_t j=0;j<3;++j)if(active[j])for(std::size_t k=0;k<3;++k)if(active[k])
    for(std::size_t l=0;l<3;++l)if(active[l])
      chemical[j][l]+=kinetic.mobility[j][k]*p.potential_hessian[k][l]/rho;
  const double gap2=std::pow(chemical[0][0]-chemical[1][1],2)+4*chemical[0][1]*chemical[1][0];
  if(gap2<0)throw std::domain_error("fingering: unresolved real composition diffusion modes");
  const double gap=std::sqrt(gap2),trace=chemical[0][0]+chemical[1][1];
  const std::array<double,2> eigenvalues{.5*(trace-gap),.5*(trace+gap)};
  if(!(eigenvalues[0]>0))throw std::domain_error("fingering: composition diffusion is not positive");
  std::array<double,2> drive{};
  const double net=-B/(ad-grad);
  if(gap>1e-10*trace) {
    // Spectral projector onto the slower chemical mode. Its complement takes
    // the remaining finite-face Ledoux contrast, retaining the same total B
    // as the thermal equations. This is a local linearization of the EOS.
    double slow=0;
    for(std::size_t j=0;j<2;++j)for(std::size_t k=0;k<2;++k)
      slow+=p.dlnRho[j]*(chemical[j][k]-(j==k?eigenvalues[1]:0))*dx[k]/(-gap);
    drive[0]=slow/thermal_buoyancy;drive[1]=net-drive[0];
  } else drive={net,0};
  const auto heat=base_.heat(index,mlo,mhi,lo,a,hi,b,false);
  const double opacity=.5*(radiation_.eval(Tlo,rlo,a).kappa+radiation_.eval(Thi,rhi,b).kappa);
  const double conductivity=4*constants::a_rad*constants::c*T*T*T/(3*rho*opacity)+heat.conductivity;
  const double thermal=conductivity/(rho*cp);
  result.prandtl=viscosity/thermal;
  result.diffusivity_ratio=eigenvalues[0]/thermal;
  if(eigenvalues[1]>=thermal)
    throw std::domain_error("fingering: chemical diffusion is outside the heat-fast closure domain");
  const auto f=two_composition_fingering(result.prandtl,
      {eigenvalues[0]/thermal,eigenvalues[1]/thermal},drive,oscillatory_);
  if(f.regime!=FingeringRegime::fingering)return result;
  result.thermal_nusselt_excess=f.thermal_nusselt_excess;
  MetalSpeciesMatrix relaxation{};
  for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k)
    relaxation[j][k]=f.wavenumber_squared*chemical[j][k]/thermal+(j==k?f.growth_rate:0);
  auto response=inverse(relaxation);
  if(f.oscillation_frequency>0) {
    MetalSpeciesMatrix square{},real_response{};
    for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k) {
      for(std::size_t l=0;l<3;++l)square[j][k]+=relaxation[j][l]*relaxation[l][k];
      if(j==k)square[j][k]+=f.oscillation_frequency*f.oscillation_frequency;
    }
    const auto inv=inverse(square);
    for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k)for(std::size_t l=0;l<3;++l)
      real_response[j][k]+=inv[j][l]*relaxation[l][k];
    response=real_response;
  }
  const double scale=thermal*f.velocity_squared;
  result.diffusivity=scale*.5*(response[0][0]+response[1][1]);
  const double Tlog=detail::logarithmic_temperature<0>(lo.lnT,hi.lnT).value;
  const double conductance=4*pi*constants::G*mass*rho*conductivity*Tlog/P;
  result.heat_luminosity=conductance*f.thermal_nusselt_excess*(grad-ad);
  const double r=.5*(std::exp(lo.lnr)+std::exp(hi.lnr)),area_mass=4*pi*r*r*rho;
  const double geometry=area_mass*area_mass/(mhi-mlo);
  result.mass_conductance=geometry*result.diffusivity;
  for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k)
    result.mixing_conductance[j][k]=geometry*scale*response[j][k];
  if(!std::isfinite(result.heat_luminosity+result.mass_conductance) || result.heat_luminosity>0)
    throw std::domain_error("fingering: invalid heat or composition flux");
  return result;
}

void BrownFingeringTransport::add_heat(MicroscopicHeatResponse& out,std::size_t i,
    double ml,double mh,const Point& lo,const Composition& a,const Point& hi,
    const Composition& b,bool derivatives) const {
  const auto f=face(i,ml,mh,lo,a,hi,b);out.carried_luminosity+=f.heat_luminosity;
  if(!derivatives || f.heat_luminosity==0)return;
  constexpr double step=2e-5;
  for(int side=0;side<2;++side)for(int k=0;k<3;++k) {
    auto plus=side?hi:lo,minus=plus;
    auto change=[&](Point& p,double h) {if(k==0)p.lnr+=h;else if(k==1)p.lnrho+=h;else p.lnT+=h;};
    change(plus,step);change(minus,-step);
    const auto fp=side?face(i,ml,mh,lo,a,plus,b):face(i,ml,mh,plus,a,hi,b);
    const auto fm=side?face(i,ml,mh,lo,a,minus,b):face(i,ml,mh,minus,a,hi,b);
    (side?out.dcarried_hi:out.dcarried_lo)[k]+=(fp.heat_luminosity-fm.heat_luminosity)/(2*step);
  }
}
MetalMicroscopicFaceResponse BrownFingeringTransport::metal_eval(std::size_t i,double ml,double mh,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  auto out=base_.metal_eval(i,ml,mh,lo,a,hi,b,derivatives);add_heat(out,i,ml,mh,lo,a,hi,b,derivatives);return out;
}
MicroscopicHeatResponse BrownFingeringTransport::heat(std::size_t i,double ml,double mh,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  auto out=base_.heat(i,ml,mh,lo,a,hi,b,derivatives);add_heat(out,i,ml,mh,lo,a,hi,b,derivatives);return out;
}
MicroscopicHeatResponse BrownFingeringTransport::heat_with_total_metal_rate(std::size_t i,double ml,double mh,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,
    const MetalSpeciesVector& rates,bool derivatives) const {
  auto out=base_.heat_with_total_metal_rate(i,ml,mh,lo,a,hi,b,rates,derivatives);
  add_heat(out,i,ml,mh,lo,a,hi,b,derivatives);return out;
}
void BrownFingeringTransport::add_mixing_flux(MetalCNFaceResponse& out,std::size_t i,double ml,double mh,
    const Point& lo,const Composition& a,const Point& hi,const Composition& b,bool derivatives) const {
  const auto f=face(i,ml,mh,lo,a,hi,b);
  if(f.mass_conductance==0)return;
  const auto flux=mixing_rates(f,a,b);
  for(std::size_t row=0;row<METAL_CN_SIZE;++row)out.rate[row]+=flux[row];
  if(!derivatives)return;
  const auto left=metal_cn_abundances(a),right=metal_cn_abundances(b);
  for(int side=0;side<2;++side)for(std::size_t col=0;col<METAL_CN_D;++col) {
    const auto& composition=side?b:a;const auto values=side?right:left;
    const double step=1e-6*std::min(values[col],composition.X[2]);
    if(step==0)continue;
    auto plus=values,minus=values;plus[col]+=step;minus[col]-=step;
    const auto cp=metal_cn_composition(composition,plus),cm=metal_cn_composition(composition,minus);
    const auto fp=side?face(i,ml,mh,lo,a,hi,cp):face(i,ml,mh,lo,cp,hi,b);
    const auto fm=side?face(i,ml,mh,lo,a,hi,cm):face(i,ml,mh,lo,cm,hi,b);
    const auto rp=side?mixing_rates(fp,a,cp):mixing_rates(fp,cp,b);
    const auto rm=side?mixing_rates(fm,a,cm):mixing_rates(fm,cm,b);
    for(std::size_t row=0;row<METAL_CN_SIZE;++row)
      (side?out.dright:out.dleft)[row][col]+=(rp[row]-rm[row])/(2*step);
  }
}
} // namespace ember

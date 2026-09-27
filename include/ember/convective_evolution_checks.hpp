#pragma once
#include "ember/eos_variable_metal.hpp"
#include "ember/constants.hpp"
#include "ember/metal_cn_transport.hpp"
#include "ember/metal_microscopic_transport.hpp"
#include "ember/controller.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ember::driver {
using HeatDirections=std::array<double,4>; // H1, He3, total Z, initial D; He4 closes the mass
inline HeatDirections trace_deuterium_enthalpy(const VariableMetalHelmholtzEos& eos,
    double T,double rho,const Composition& composition) {
  const double f=1-composition[Species::H2]/2;auto c=composition;
  c.X[0]+=c[Species::H2]/2;c[Species::H2]=0;
  for(auto& x:c.X)x/=f;
  if(c.cn_molality)for(auto& x:*c.cn_molality)x/=f;
  const auto response=eos.eval_with_derivatives(T,f*rho,c);const auto& e=response.state;
  const auto potential=eos.composition_potential(T,f*rho,c,{true,true,true});
  HeatDirections h{};double ed=-.5*(e.E+response.dE_dlnRho),pd=-e.P*e.chiRho/(2*f);
  const std::array<double,3> direction{1+c.X[0],c.X[1],c.Z()};
  for(std::size_t k=0;k<3;++k) {
    const double ek=-T*potential.dgradient_dlnT[k],pk=f*rho*T*potential.dgradient_dlnRho[k];
    h[k]=ek+e.delta*pk/(f*rho);ed+=.5*direction[k]*ek;pd+=direction[k]*pk/(2*f);
  }
  h[3]=ed+e.delta*pd/rho;return h;
}
inline MetalCNVector physical_source(const PPCNNetwork& nuclear,double T,double rho,const Composition& c) {
  const auto light=nuclear.light().eval(T,rho,c);
  const auto cn=nuclear.cn().response(T,rho,c,*c.cn_molality).physical.state;
  return {light.dXdt[0]+cn.dXdt[0],light.dXdt[1]+cn.dXdt[1],cn.dXdt[3],cn.dXdt[4],cn.dXdt[5],0,light.dXdt[8]};
}
struct HomogeneousCheck {
  double maximum_heat_fraction{},maximum_gross_heat_fraction{},travel_years{},D_gradient_estimate{};
  double burn_gradient_estimate{},drift_gradient_proxy{},kinetic_heat_proxy{};
  double convective_mass_fraction{1},maximum_relative_mixing_gradient{};
  std::size_t radiative_boundaries{};
};
// Instantaneous mixing is selected only during initial-D contraction. A
// completely mixed star has no boundary through which microscopic settling
// can change its global inventory. This does not bound microscopic kinetic
// heat or establish a partially ionized diffusion law.
inline HomogeneousCheck check_initial_convection(const Model& m,const Physics& physics,
    const VariableMetalHelmholtzEos& eos,const PPCNNetwork& nuclear,bool omit_material_heat=true) {
  const auto regions=convective_mixing_regions(m,physics);
  if(regions.size()!=1 || regions.front()!=std::pair<std::size_t,std::size_t>{0,m.size()})
    throw std::domain_error("initial-D approximation requires whole-star convection");
  for(const auto& c:m.comp)if(c!=m.comp.front())
    throw std::domain_error("initial-D approximation requires homogeneous composition");
  if(m.comp.front()[Species::H2]>1e-4)throw std::domain_error("initial D exceeds material approximation");
  const auto w=nodal_mass_weights(m);const auto mixing=convective_mixing_faces(m,physics);
  std::vector<HeatDirections> source(m.size());std::array<long double,4> total{},flux{};
  double minimum_D_time=std::numeric_limits<double>::infinity();
  for(std::size_t i=0;i<m.size();++i) {
    const auto s=physical_source(nuclear,m.T(i),m.rho(i),m.comp[i]);
    source[i]={s[0],s[1],s[2]+s[3]+s[4]+s[5],s[6]};
    for(std::size_t k=0;k<4;++k)total[k]+=static_cast<long double>(w[i])*source[i][k];
    if(s[6]<0)minimum_D_time=std::min(minimum_D_time,-m.comp[i][Species::H2]/s[6]);
  }
  HomogeneousCheck out;
  HeatDirections burn_gradient{};
  for(std::size_t i=0;i+1<m.size();++i) {
    for(std::size_t k=0;k<4;++k)flux[k]+=static_cast<long double>(w[i])*(source[i][k]-total[k]/m.M);
    const auto h=trace_deuterium_enthalpy(eos,std::sqrt(m.T(i)*m.T(i+1)),.5*(m.rho(i)+m.rho(i+1)),m.comp[i]);
    double Q=0,gross=0;for(std::size_t k=0;k<4;++k){const double q=h[k]*static_cast<double>(flux[k]);Q+=q;gross+=std::abs(q);}
    const double L=m.y[i].L;
    if(!(L>0) || !(mixing[i].velocity>0))throw std::domain_error("initial-D convection lacks outward luminosity or finite mixing velocity");
    out.maximum_heat_fraction=std::max(out.maximum_heat_fraction,std::abs(Q)/L);
    out.maximum_gross_heat_fraction=std::max(out.maximum_gross_heat_fraction,gross/L);
    out.travel_years+=(m.r(i+1)-m.r(i))/mixing[i].velocity/31557600.;
    const double T=std::sqrt(m.T(i)*m.T(i+1)),rho=.5*(m.rho(i)+m.rho(i+1));
    const double r=.5*(m.r(i)+m.r(i+1)),mass=.5*(m.m[i]+m.m[i+1]),dr=m.r(i+1)-m.r(i);
    const auto e=physics.eos->eval(T,rho,m.comp[i]);const double Hp=e.P/(rho*constants::G*mass/(r*r));
    for(std::size_t k=0;k<4;++k)
      burn_gradient[k]+=std::abs(static_cast<double>(flux[k]))*dr/(4*M_PI*r*r*rho*mixing[i].diffusivity);
    // Explicit order-of-magnitude uncertainty bracket, not a cold collision
    // law. Use the larger neutral (sigma=1e-16 cm2) and Coulomb (ln Lambda=1)
    // mobility, even where only one applies. A=56 overestimates the settling
    // force of the bulk H/He mixture. Conductive/kinetic heat sensitivity is
    // tested separately on actual stellar structures.
    const double v=std::sqrt(constants::kB*T/constants::amu),number=rho/constants::amu,e2=2.307077552e-19;
    const double mobility=std::max(v/(number*1e-16),v*std::pow(constants::kB*T,2)/(number*e2*e2));
    out.drift_gradient_proxy+=56*mobility/mixing[i].diffusivity*dr/Hp;
    out.kinetic_heat_proxy=std::max(out.kinetic_heat_proxy,
        4*M_PI*r*r*rho*mobility/Hp*(20*constants::kB*T/constants::amu)/L);
  }
  out.burn_gradient_estimate=*std::max_element(burn_gradient.begin(),burn_gradient.end());
  out.D_gradient_estimate=m.comp.front()[Species::H2]*out.travel_years*31557600./minimum_D_time;
  if((omit_material_heat && out.maximum_gross_heat_fraction>1e-4) || out.D_gradient_estimate>1e-8)
    throw std::domain_error("initial-D homogeneous/heat approximation exceeds measured early-phase bounds");
  if(!omit_material_heat && (out.burn_gradient_estimate>1e-8
        || out.drift_gradient_proxy>1e-6 || out.kinetic_heat_proxy>.05))
    throw std::domain_error("whole-star convective approximation exceeds assessed transport range");
  return out;
}

// Assess the instantaneous-mixing approximation region by region. Species
// crossing a radiative boundary use the actual microscopic provider. Within
// a mixed region, flux/(rho D) estimates the gradient needed for convection
// to carry the reconstructed total rate. The omitted cool drift/heat retain
// the same order-of-magnitude estimates used for the wholly convective star.
inline HomogeneousCheck check_envelope_transport(const Model& m,const Physics& physics,
    const MetalMicroscopicTransport& transport,std::span<const MetalSpeciesVector> rates,
    double minimum_microscopic_T,double maximum_relative_gradient) {
  if(rates.size()+1!=m.size() || !(maximum_relative_gradient>0))
    throw std::domain_error("envelope transport assessment: missing rates or invalid mixing allowance");
  for(const auto& c:m.comp)if(c[Species::H2]!=0)
    throw std::domain_error("envelope transport assessment requires initial D exhaustion");
  const auto regions=convective_mixing_regions(m,physics);
  const auto mixing=convective_mixing_faces(m,physics);const auto weights=nodal_mass_weights(m);
  HomogeneousCheck out;out.convective_mass_fraction=0;
  for(const auto [begin,end]:regions) {
    if(end-begin>1) {
      std::array<double,3> gradient{};double travel=0;
      for(std::size_t i=begin;i<end;++i) {
        if(m.comp[i]!=m.comp[begin])throw std::domain_error("instantaneous convective region is not homogeneous");
        out.convective_mass_fraction+=weights[i]/m.M;
      }
      for(std::size_t i=begin;i+1<end;++i) {
        const double T=std::sqrt(m.T(i)*m.T(i+1)),rho=.5*(m.rho(i)+m.rho(i+1));
        const double r=.5*(m.r(i)+m.r(i+1)),mass=.5*(m.m[i]+m.m[i+1]),dr=m.r(i+1)-m.r(i);
        if(!(mixing[i].velocity>0 && mixing[i].diffusivity>0))
          throw std::domain_error("convective region lacks finite mixing");
        travel+=dr/mixing[i].velocity/31557600.;
        auto convection_rate=rates[i];
        if(std::min(m.T(i),m.T(i+1))>=minimum_microscopic_T) {
          const auto micro=transport.metal_eval(i,m.m[i],m.m[i+1],m.y[i],m.comp[i],m.y[i+1],m.comp[i+1],false);
          for(std::size_t k=0;k<3;++k)convection_rate[k]-=micro.species.rate[k];
        }else {
          const auto e=physics.eos->eval(T,rho,m.comp[i]);const double Hp=e.P/(rho*constants::G*mass/(r*r));
          const double v=std::sqrt(constants::kB*T/constants::amu),number=rho/constants::amu,e2=2.307077552e-19;
          const double mobility=std::max(v/(number*1e-16),v*std::pow(constants::kB*T,2)/(number*e2*e2));
          out.drift_gradient_proxy+=56*mobility/mixing[i].diffusivity*dr/Hp;
          out.kinetic_heat_proxy=std::max(out.kinetic_heat_proxy,
              4*M_PI*r*r*rho*mobility/Hp*(20*constants::kB*T/constants::amu)/std::max(std::abs(m.y[i].L),1.));
        }
        for(std::size_t k=0;k<3;++k)
          gradient[k]+=std::abs(convection_rate[k])*dr/(4*M_PI*r*r*rho*mixing[i].diffusivity);
      }
      out.travel_years=std::max(out.travel_years,travel);
      for(std::size_t k=0;k<3;++k) {
        out.burn_gradient_estimate=std::max(out.burn_gradient_estimate,gradient[k]);
        const double abundance=k==2?m.comp[begin].Z():m.comp[begin].X[k];
        out.maximum_relative_mixing_gradient=std::max(out.maximum_relative_mixing_gradient,
            gradient[k]/std::max(abundance,1e-12));
      }
    }
    if(end<m.size()) {
      const auto i=end-1;++out.radiative_boundaries;
      // No cool fallback for species: this verifies the actual law's domain.
      (void)transport.metal_eval(i,m.m[i],m.m[i+1],m.y[i],m.comp[i],m.y[i+1],m.comp[i+1],false);
    }
  }
  if(out.maximum_relative_mixing_gradient>maximum_relative_gradient)
    throw std::domain_error("instantaneous convection requires finite mixing at this abundance gradient");
  if(out.drift_gradient_proxy>1e-6 || out.kinetic_heat_proxy>.05)
    throw std::domain_error("cool convective envelope exceeds assessed drift/heat range");
  return out;
}

// Global baryonic species budgets are independent of internal redistribution.
// Evaluate physical CN isotope sources, not the material table's GS98 slots.
inline EvolutionAudit check_interval(const Model& old,const EvolutionStep& step,double dt,
    const PPCNNetwork& nuclear,double abundance_tolerance) {
  EvolutionAudit out;if(!step.converged || old.m!=step.model.m)return out;
  const auto& next=step.model;const auto w=nodal_mass_weights(next);
  std::array<long double,METAL_CN_SIZE> change{},source{},roundoff{};
  long double mass_power=0,mass_roundoff=0,heat=0,neutrinos=0;
  constexpr std::array<std::size_t,METAL_CN_SIZE> slot{0,1,3,4,5,7,8};
  const auto ulp=[](double x){return std::nextafter(x,std::numeric_limits<double>::infinity())-x;};
  const long double c2=static_cast<long double>(constants::c)*constants::c;
  for(std::size_t i=0;i<next.size();++i) {
    const auto a=metal_cn_abundances(old.comp[i]),b=metal_cn_abundances(next.comp[i]);
    const auto s=physical_source(nuclear,next.T(i),next.rho(i),next.comp[i]);
    long double rest_change=0,rest_roundoff=0;
    for(std::size_t k=0;k<METAL_CN_SIZE;++k) {
      const long double d=static_cast<long double>(b[k])-a[k];
      change[k]+=w[i]/next.M*d;source[k]+=static_cast<long double>(w[i])/next.M*dt*s[k];
      const double rounding=2*(ulp(a[k])+ulp(b[k]));
      roundoff[k]+=w[i]/next.M*rounding;
      // Inert mass is conserved exactly; its relative mixture is fixed.
      const long double binding=(nuclides[slot[k]].A/mass_numbers[slot[k]]-nuclides[2].A/4)*c2;
      rest_change+=binding*d;rest_roundoff+=std::abs(binding)*rounding;
    }
    mass_power-=w[i]*rest_change/dt;mass_roundoff+=w[i]*rest_roundoff/dt;
    const auto n=nuclear.eval(next.T(i),next.rho(i),next.comp[i]);heat+=w[i]*n.eps;neutrinos+=w[i]*n.eps_neutrino;
  }
  bool species=true;
  for(std::size_t k=0;k<METAL_CN_SIZE;++k) {
    const double error=static_cast<double>(std::abs(change[k]-source[k]));
    out.maximum_species_error=std::max(out.maximum_species_error,error);
    species &= error<=abundance_tolerance+roundoff[k];
  }
  const double L=next.y.back().L;const long double released=heat+neutrinos;
  out.mass_error_surface=static_cast<double>((mass_power-released)/L);
  out.source_error=static_cast<double>((heat-step.nuclear_luminosity)/L);
  out.first_law=std::abs(step.luminosity_balance);
  out.pass=species && std::abs(mass_power-released)<=2e-6L*std::abs(released)+mass_roundoff
    && std::abs(out.mass_error_surface)<2e-7 && std::abs(out.source_error)<2e-7 && out.first_law<2e-7;
  return out;
}
} // namespace ember::driver

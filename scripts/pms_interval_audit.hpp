#pragma once
#include "ember/deuterium.hpp"
#include "ember/evolution.hpp"
#include <cmath>
#include <limits>

// Fully mixed pp + initial-D intervals. This deliberately refuses partial
// mixing or metal transport; those need regional, rather than global, checks.
// Source powers and isotope changes are evaluated separately. Nuclear neutrino
// energy belongs in the mass defect, but not in the deposited heat.
struct PMSIntervalAudit {
  bool pass=false, species=false, deuterium=false, source_energy=false;
  double fuel_error=0, mass_error=0, roundoff=0, first_law=0;
  double mass_error_surface=0, pp_heat=0, D_heat=0, source_heat_error=0;
};
inline PMSIntervalAudit audit_pms_interval(const ember::Model& before,
    const ember::EvolutionStep& step,double dt,const ember::PPDeuterium&) {
  using namespace ember;
  PMSIntervalAudit out;
  if(!step.converged || before.m!=step.model.m || !std::isfinite(dt) || dt<=0)
    return out;
  const auto& after=step.model;const auto w=nodal_mass_weights(after);
  const PPChains pp(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn);
  std::array<long double,NSPEC> pp_increment{};
  long double fuel=0,Dpower=0,pp_heat=0,pp_neutrino=0,masspower=0,roundoff=0;
  bool species=step.convective_mass_fraction>.99999;
  auto ulp=[](double x){return std::nextafter(x,std::numeric_limits<double>::infinity())-x;};
  for(std::size_t i=0;i<after.size();++i) {
    const auto source=pp.eval(after.T(i),after.rho(i),after.comp[i]);
    for(std::size_t k=0;k<NSPEC;++k)
      pp_increment[k]+=(long double)dt*w[i]/after.M*source.dXdt[k];
    pp_heat+=(long double)w[i]*source.eps;
    pp_neutrino+=(long double)w[i]*source.eps_neutrino;
    auto other=after.comp[i];other[Species::H2]=0;
    Dpower+=(long double)w[i]*deuterium_capture(after.T(i),after.rho(i),other,
                                              after.comp[i][Species::H2]).source.eps;
  }
  for(std::size_t i=0;i<after.size();++i) {
    const auto& a=before.comp[i];const auto& b=after.comp[i];
    species &= a.X==before.comp.front().X && b.X==after.comp.front().X;
    species &= b.cn_molality==std::nullopt && b.basis==AbundanceBasis::baryon_mass
      && std::abs(b.sum()-1)<2e-12;
    for(std::size_t k=0;k<NSPEC;++k) {
      species &= std::isfinite(b.X[k]) && b.X[k]>=0;
      if(is_metal_species(k))species &= b.X[k]==a.X[k];
      const long double binding=(nuclides[k].A/mass_numbers[k]-1)*constants::c*constants::c;
      masspower-=(long double)w[i]*binding*(b.X[k]-a.X[k])/dt;
      roundoff+=(long double)w[i]*std::abs(binding)*.5*(ulp(a.X[k])+ulp(b.X[k]))/dt;
    }
    const double d=b[Species::H2]-a[Species::H2];fuel-=(long double)w[i]*d;
    species &= std::abs((long double)(b.X[0]-a.X[0])-.5L*d-pp_increment[0])
      <=.5*(ulp(a.X[0])+ulp(b.X[0]))+.25*(ulp(a[Species::H2])+ulp(b[Species::H2]));
    species &= std::abs((long double)(b.X[1]-a.X[1])+1.5L*d-pp_increment[1])
      <=.5*(ulp(a.X[1])+ulp(b.X[1]))+.75*(ulp(a[Species::H2])+ulp(b[Species::H2]));
    // He4 is the conserved complement. Include the stored rounding of each
    // changed constituent when checking its independent pp production.
    double helium_roundoff=0;
    for(auto k:{0u,1u,2u,8u})helium_roundoff+=.5*(ulp(a.X[k])+ulp(b.X[k]));
    species &= std::abs((long double)(b.X[2]-a.X[2])-pp_increment[2])<=helium_roundoff;
  }
  const double Q=(nuclides[0].A+nuclides[8].A-nuclides[1].A)*constants::c*constants::c/2;
  out.deuterium=fuel==0 && Dpower==0;
  if(fuel>0) {
    out.fuel_error=Dpower*dt/(fuel*Q)-1;
    out.deuterium=std::abs(out.fuel_error)<2e-6;
  }
  const double release=step.nuclear_luminosity+step.neutrino_luminosity;
  const long double error=masspower-release;
  out.pp_heat=pp_heat;out.D_heat=Dpower;
  out.source_heat_error=(pp_heat+Dpower-step.nuclear_luminosity)/std::max(release,std::numeric_limits<double>::min());
  out.source_energy=std::abs(out.source_heat_error)<2e-6
    && std::abs(pp_neutrino-step.neutrino_luminosity)<=2e-6*release;
  out.mass_error=error/std::max(release,std::numeric_limits<double>::min());
  out.roundoff=roundoff/std::max(release,std::numeric_limits<double>::min());
  out.mass_error_surface=error/after.y.back().L;out.first_law=step.luminosity_balance;
  out.species=species;
  out.pass=species && out.deuterium && out.source_energy
    && std::abs(error)<=2e-6*release+roundoff
    // Same quantity, same representation limit. `error` can only be formed to within `roundoff`, the
    // stored-abundance rounding already accumulated above, so both bounds on it carry that allowance;
    // the surface bound previously omitted it and so tightened without limit as the step shrank, since
    // the composition change per step falls with dt while its rounding does not. Measured at 220.72 Myr:
    // the allowance times dt is constant to five figures across 10000, 3750 and 2500 yr steps, which is
    // the signature of rounding rather than of a physical error. The 2e-7 coefficient is unchanged, and
    // near the main sequence, where release approaches the surface luminosity, this bound stays the
    // tighter of the two by an order of magnitude.
    && std::abs(error)<=2e-7*after.y.back().L+roundoff && std::abs(out.first_law)<2e-7;
  return out;
}

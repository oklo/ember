#pragma once
#include "ember/convective_evolution_checks.hpp"
#include "ember/convective_material_heat.hpp"
#include "ember/envelope_transport.hpp"
#include "ember/boundary.hpp"
#include <iomanip>
#include <sstream>

namespace ember::driver {
// Physical assessment for the lifetime driver's optional second-order state.
// Every constituent solve and the full/two-half time check remain mandatory.
struct LifetimeExtrapolation {
  const EvolutionControlHooks& hooks;
  const EvolutionOptions& options;
  const Atmosphere& atmosphere;
  const PPCNNetwork& nuclear;
  ConvectiveMaterialHeat& convective_heat;
  EnvelopeTransport& envelope_heat;
  std::size_t threads;
  double residual_limit,energy_tolerance,species_tolerance;
  std::ostream& log;
  std::string operator()(const EvolutionState& initial,const EvolutionStep& full,const EvolutionStep& h1,
      const EvolutionStep& h2,const Model& R,const std::vector<std::array<double,3>>& rates,double dt) const {
      const auto& start=initial.model;
      const auto has_D=[](const Model& m){return std::any_of(m.comp.begin(),m.comp.end(),[](const auto& c){return c[Species::H2]!=0;});};
      std::string reason;
      // The accounting below covers the metal-CN slots, not the initial-D isotope.
      if(has_D(start) || has_D(R))return "initial deuterium present: Richardson not assessed for that isotope";
      const auto with_rates=[&](std::span<const std::array<double,3>> selected_rates,auto action) {
        const auto old_convective=convective_heat.diagnostic_rates,old_envelope=envelope_heat.diagnostic_rates;
        struct Restore {
          ConvectiveMaterialHeat& c;EnvelopeTransport& e;
          std::span<const std::array<double,3>> old_c,old_e;
          ~Restore(){c.diagnostic_rates=old_c;e.diagnostic_rates=old_e;}
        } restore{convective_heat,envelope_heat,old_convective,old_envelope};
        convective_heat.diagnostic_rates=selected_rates;envelope_heat.diagnostic_rates=selected_rates;
        return action();
      };
      auto partitions=[&](const Model& m,std::span<const std::array<double,3>> r) {
        return with_rates(r,[&]{
          auto selected=options;hooks.configure_step(m,selected);const auto& p=hooks.physics(m);
          return std::pair{convective_mixing_regions(m,p,threads),instantaneous_mixing_regions(m,p,selected)};
        });
      };
      bool partitions_equal=false;
      try {
        const auto current=partitions(R,rates);
        partitions_equal=partitions(start,initial.metal_heat_rates)==current
          && partitions(full.model,full.total_metal_species_rates)==current
          && partitions(h1.model,h1.total_metal_species_rates)==current
          && partitions(h2.model,h2.total_metal_species_rates)==current;
        if(!partitions_equal)reason="convective membership changed during the interval";
      }catch(const std::exception& e){reason=std::string("partition evaluation failed: ")+e.what();}
      // 2. Support and algebraic residuals (the backward-Euler thermal rows are not gated):
      // hydrostatic and mass rows (relative), transport row scaled by zone mass width
      // (a ln T mismatch), surface ln T / ln P against the atmosphere, central radius row.
      struct Residuals {double hydro{},mass{},transport{},surface_T{},surface_P{},central_r{};};
      auto residuals=[&](const Model& m,double step,const Model& before) {
        const auto& phys=hooks.physics(m);Residuals out;
        for(std::size_t i=0;i+1<m.size();++i) {
          const auto f=zone_equations(m,i,phys,step,&before);
          const double r=.5*(m.r(i)+m.r(i+1)),rho=.5*(m.rho(i)+m.rho(i+1)),mb=.5*(m.m[i]+m.m[i+1]),dm=m.m[i+1]-m.m[i];
          const double P=.5*(phys.eos->eval(m.T(i),m.rho(i),m.comp[i]).P+phys.eos->eval(m.T(i+1),m.rho(i+1),m.comp[i+1]).P);
          out.mass=std::max(out.mass,std::abs(f[0])*4*M_PI*r*r*r*rho);
          out.hydro=std::max(out.hydro,std::abs(f[1])*4*M_PI*r*r*r*r*P/(constants::G*mb));
          out.transport=std::max(out.transport,std::abs(f[3])*dm);
        }
        const auto surface=surface_residual(m.y.back(),m.M,m.comp.back(),*phys.eos,atmosphere);
        out.surface_T=std::abs(surface.f[0]);out.surface_P=std::abs(surface.f[1]);
        out.central_r=std::abs(central_residual(m,phys,step,&before).f[0]);
        return out;
      };
      Residuals cand,half;
      try {
        const auto& phys=hooks.physics(R);
        for(std::size_t i=0;i<R.size();++i) {
          for(double v:{R.y[i].lnr,R.y[i].lnrho,R.y[i].lnT,R.y[i].L})
            if(!std::isfinite(v))throw std::domain_error("non-finite extrapolated state");
          (void)phys.eos->eval(R.T(i),R.rho(i),R.comp[i]);(void)phys.opacity->eval(R.T(i),R.rho(i),R.comp[i]);
          (void)nuclear.eval(R.T(i),R.rho(i),R.comp[i]);
        }
        cand=with_rates(rates,[&]{return residuals(R,dt,start);});
        half=with_rates(h2.total_metal_species_rates,[&]{return residuals(h2.model,dt/2,h1.model);});
        for(double v:{cand.hydro,cand.mass,cand.transport,cand.surface_T,cand.surface_P,cand.central_r})
          if(!std::isfinite(v))throw std::domain_error("non-finite residual");
        if(reason.empty() && std::max({cand.hydro,cand.mass,cand.transport,cand.surface_T,cand.surface_P,cand.central_r})>residual_limit)
          reason="an algebraic residual of the extrapolated state exceeds the limit";
        hooks.assess(R,rates);
      }catch(const std::exception& e){if(reason.empty())reason=std::string("extrapolated state assessment failed: ")+e.what();}
      // 3. Accounting.
      struct Budget {std::array<long double,METAL_CN_SIZE> inventory{},source{};long double energy{},nuclear_power{};double emitted{};};
      constexpr std::array<std::size_t,METAL_CN_SIZE> slot{0,1,3,4,5,7,8};
      const long double c2=static_cast<long double>(constants::c)*constants::c;
      auto budget=[&](const Model& m) {
        Budget b;const auto w=nodal_mass_weights(m);long double Lnu=0,Lth=0;
        for(std::size_t i=0;i<m.size();++i) {
          const auto a=metal_cn_abundances(m.comp[i]);const auto s=physical_source(nuclear,m.T(i),m.rho(i),m.comp[i]);
          const auto n=nuclear.eval(m.T(i),m.rho(i),m.comp[i]);
          const auto& ph=hooks.physics(m);
          for(std::size_t k=0;k<METAL_CN_SIZE;++k) {
            b.inventory[k]+=static_cast<long double>(w[i])/m.M*a[k];b.source[k]+=static_cast<long double>(w[i])/m.M*s[k];
            b.energy+=w[i]*(nuclides[slot[k]].A/mass_numbers[slot[k]]-nuclides[2].A/4)*c2*a[k];
          }
          b.energy+=w[i]*(static_cast<long double>(ph.eos->eval(m.T(i),m.rho(i),m.comp[i]).E)-constants::G*m.m[i]/m.r(i));
          b.nuclear_power+=w[i]*n.eps;Lnu+=w[i]*n.eps_neutrino;
          if(ph.neutrino_losses)Lth+=w[i]*evaluate_losses(ph.neutrino_losses,m.T(i),m.rho(i),m.comp[i]).eps;
        }
        b.emitted=static_cast<double>(m.y.back().L+Lnu+Lth);
        return b;
      };
      std::ostringstream line;line<<std::setprecision(6);
      try {
        const auto S=budget(start),F=budget(full.model),A=budget(h1.model),B=budget(h2.model),X=budget(R);
        auto species_error=[&](const Budget& begin,const Budget& end,double step) {
          long double e=0;for(std::size_t k=0;k<METAL_CN_SIZE;++k)
            e=std::max(e,std::abs(end.inventory[k]-begin.inventory[k]-step*end.source[k]));
          return static_cast<double>(e);
        };
        auto energy_error=[&](const Budget& begin,const Budget& end,double step) {
          return static_cast<double>(end.energy-begin.energy+step*end.emitted);};
        long double sR=0,sT=0,change=0;
        for(std::size_t k=0;k<METAL_CN_SIZE;++k) {
          const long double d=X.inventory[k]-S.inventory[k];change=std::max(change,std::abs(d));
          sR=std::max(sR,std::abs(d-(2*(dt/2*A.source[k]+dt/2*B.source[k])-dt*F.source[k])));
          sT=std::max(sT,std::abs(d-dt/2*(S.source[k]+X.source[k])));
        }
        const long double eR=X.energy-S.energy+(2*(dt/2*A.emitted+dt/2*B.emitted)-dt*F.emitted);
        const long double eT=X.energy-S.energy+dt/2*(S.emitted+X.emitted);
        const long double nonlinear=X.energy-(2*B.energy-F.energy);
        const double emitted=dt*static_cast<double>(S.emitted),nuclear_energy=dt*static_cast<double>(S.nuclear_power);
        if(!std::isfinite(static_cast<double>(sR)) || !std::isfinite(static_cast<double>(sT)) || !std::isfinite(static_cast<double>(eR))
           || !std::isfinite(static_cast<double>(eT)) || !std::isfinite(emitted) || !(emitted>0))
          throw std::domain_error("non-finite candidate accounting");
        if(reason.empty() && std::abs(static_cast<double>(eR))>energy_tolerance*std::max(emitted,std::abs(nuclear_energy)))
          reason="extrapolated global energy budget exceeds its allowance";
        if(reason.empty() && static_cast<double>(sT)>species_tolerance)
          reason="extrapolated species budget exceeds the time accuracy";
        line<<"{\"age_years\":"<<start.age/31557600.<<",\"dt_years\":"<<dt/31557600.<<",\"accepted\":"<<reason.empty()
          <<",\"reason\":"<<std::quoted(reason)<<",\"partitions_equal\":"<<partitions_equal
          <<",\"residuals\":{\"candidate\":["<<cand.hydro<<','<<cand.mass<<','<<cand.transport<<','<<cand.surface_T<<','<<cand.surface_P<<','<<cand.central_r
          <<"],\"h2\":["<<half.hydro<<','<<half.mass<<','<<half.transport<<','<<half.surface_T<<','<<half.surface_P<<','<<half.central_r
          <<"],\"order\":\"hydrostatic,mass,transport_dlnT,surface_lnT,surface_lnP,central_radius\"}"
          <<",\"inventory_change\":"<<static_cast<double>(change)
          <<",\"species_error\":{\"full\":"<<species_error(S,F,dt)<<",\"h1\":"<<species_error(S,A,dt/2)<<",\"h2\":"<<species_error(A,B,dt/2)
          <<",\"extrapolated_vs_richardson_quadrature\":"<<static_cast<double>(sR)<<",\"extrapolated_vs_trapezoid\":"<<static_cast<double>(sT)<<"}"
          <<",\"energy_error_over_emitted\":{\"full\":"<<energy_error(S,F,dt)/emitted<<",\"h1\":"<<energy_error(S,A,dt/2)/emitted
          <<",\"h2\":"<<energy_error(A,B,dt/2)/emitted<<",\"extrapolated_vs_richardson_quadrature\":"<<static_cast<double>(eR)/emitted
          <<",\"extrapolated_vs_trapezoid\":"<<static_cast<double>(eT)/emitted<<",\"nonlinearity\":"<<static_cast<double>(nonlinear)/emitted<<"}"
          <<",\"emitted_energy_erg\":"<<emitted<<",\"nuclear_energy_erg\":"<<nuclear_energy<<"}";
      }catch(const std::exception& e){
        if(reason.empty())reason=std::string("candidate accounting failed: ")+e.what();
        line.str("");line<<"{\"age_years\":"<<start.age/31557600.<<",\"accepted\":0,\"reason\":"<<std::quoted(reason)<<"}";
      }
      if(log)log<<line.str()<<'\n'<<std::flush;
      return reason;
  }
};
} // namespace ember::driver

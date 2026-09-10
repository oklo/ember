#include "ember/evolution.hpp"
#include "ember/constants.hpp"
#include "ember/convection.hpp"
#include "energy.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
void validate_composition(const Composition& c) {
  if(c.basis!=AbundanceBasis::baryon_mass || std::abs(c.sum()-1)>1e-10)
    throw std::invalid_argument("evolution requires normalized baryonic abundances");
  for(double x:c.X) if(!std::isfinite(x) || x<0)
    throw std::invalid_argument("evolution requires finite nonnegative abundances");
}
}
std::vector<double> nodal_mass_weights(const Model& m) {
  if(m.size()<2 || m.m.size()!=m.size() || !(m.m[0]>0) || m.m.back()!=m.M)
    throw std::invalid_argument("nodal_mass_weights: invalid mass mesh");
  std::vector<double> w(m.size());w[0]=m.m[0];
  for(std::size_t i=1;i<m.size();++i) {
    const double dm=m.m[i]-m.m[i-1];
    if(!(dm>0) || !std::isfinite(dm)) throw std::invalid_argument("nodal_mass_weights: unordered mesh");
    w[i-1]+=.5*dm;w[i]+=.5*dm;
  }
  return w;
}

MixingRegions schwarzschild_mixing_regions(const Model& m,const Physics& p) {
  auto selected=p;selected.criterion=ConvectiveCriterion::schwarzschild;
  return convective_mixing_regions(m,selected);
}

MixingRegions convective_mixing_regions(const Model& m,const Physics& p) {
  if(!p.eos || !p.opacity || m.comp.size()!=m.size())
    throw std::invalid_argument("mixing regions: invalid model or physics");
  nodal_mass_weights(m);
  struct Local {double P,kappa,ad,delta;};
  std::vector<Local> local;
  for(std::size_t i=0;i<m.size();++i) {
    const auto e=p.eos->eval(m.T(i),m.rho(i),m.comp[i]);
    local.push_back({e.P,p.opacity->eval(m.T(i),m.rho(i),m.comp[i]).kappa,e.grad_ad,e.delta});
  }
  MixingRegions regions;std::size_t first=0;
  for(std::size_t i=0;i+1<m.size();++i) {
    // Same arithmetic midpoint transport quantities as zone assembly.
    const double T=.5*(m.T(i)+m.T(i+1)),mass=.5*(m.m[i]+m.m[i+1]);
    const double P=.5*(local[i].P+local[i+1].P),k=.5*(local[i].kappa+local[i+1].kappa);
    const double L=.5*(m.y[i].L+m.y[i+1].L),ad=.5*(local[i].ad+local[i+1].ad);
    const double rad=3*k*L*P/(16*M_PI*constants::a_rad*constants::c*constants::G*mass*std::pow(T,4));
    double B=0;
    if(p.criterion==ConvectiveCriterion::ledoux)
      B=composition_buoyancy(*p.eos,T,P,.5*(local[i].delta+local[i+1].delta),
          std::log(local[i+1].P)-std::log(local[i].P),m.comp[i],m.comp[i+1],.5*(m.rho(i)+m.rho(i+1))).B;
    if(!(rad>ad+B)) {regions.emplace_back(first,i+1);first=i+1;}
  }
  regions.emplace_back(first,m.size());return regions;
}

std::vector<Composition> burn_and_mix(const Model& thermal,const Model& previous,
    const Nuclear& nuclear,const MixingRegions& regions,double dt,double tolerance) {
  if(!std::isfinite(dt) || dt<=0 || !std::isfinite(tolerance) || tolerance<=0
      || thermal.m!=previous.m || thermal.size()!=previous.size()
      || thermal.comp.size()!=thermal.size() || previous.comp.size()!=previous.size())
    throw std::invalid_argument("burn_and_mix: invalid step, tolerance or previous model");
  const auto weights=nodal_mass_weights(thermal);
  for(const auto& c:previous.comp) validate_composition(c);
  std::vector<Composition> output=previous.comp;
  std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())
      throw std::invalid_argument("burn_and_mix: regions must partition the mesh");
    next=end;
    Composition old{};old.basis=AbundanceBasis::baryon_mass;
    old.metal_inventory=previous.comp[begin].metal_inventory;double mass=0;
    for(std::size_t i=begin;i<end;++i) {
      if(previous.comp[i].metal_inventory!=old.metal_inventory)
        throw std::invalid_argument("burn_and_mix: region has inconsistent metal inventories");
      mass+=weights[i];for(std::size_t j=0;j<NSPEC;++j) old.X[j]+=weights[i]*previous.comp[i].X[j];
    }
    for(double& v:old.X) v/=mass;
    auto candidate=old;
    const double available=1-old.Z();
    candidate.X[2]=available-candidate.X[0]-candidate.X[1];
    auto equation=[&](const Composition& c,bool jacobian) {
      std::array<double,6> f{c.X[0]-old.X[0],c.X[1]-old.X[1],1,0,0,1};
      for(std::size_t i=begin;i<end;++i) {
        const double factor=dt*weights[i]/mass;
        NuclearResponse response;
        if(jacobian) response=nuclear.composition_response(thermal.T(i),thermal.rho(i),c);
        else response.state=nuclear.eval(thermal.T(i),thermal.rho(i),c);
        for(int row=0;row<2;++row) {
          f[row]-=factor*response.state.dXdt[row];
          if(jacobian) for(int col=0;col<2;++col)
            f[2+2*row+col]-=factor*(response.d_dXdt_dX[row][col]-response.d_dXdt_dX[row][2]);
        }
      }
      return f;
    };
    bool converged=false;
    for(int iteration=0;iteration<50;++iteration) {
      const auto f=equation(candidate,true);const double norm=std::max(std::abs(f[0]),std::abs(f[1]));
      if(!std::isfinite(norm)) throw std::runtime_error("burn_and_mix: nonfinite rates");
      if(norm<tolerance) {converged=true;break;}
      const double det=f[2]*f[5]-f[3]*f[4];
      if(!std::isfinite(det) || std::abs(det)<1e-30) throw std::runtime_error("burn_and_mix: singular abundance Jacobian");
      const double d0=(-f[0]*f[5]+f[1]*f[3])/det,d1=(-f[1]*f[2]+f[0]*f[4])/det;
      double damping=1;bool accepted=false;
      for(int trial=0;trial<40;++trial,damping*=.5) {
        auto c=candidate;c.X[0]+=damping*d0;c.X[1]+=damping*d1;c.X[2]=available-c.X[0]-c.X[1];
        if(c.X[0]<0 || c.X[1]<0 || c.X[2]<0) continue;
        const auto residual=equation(c,false);
        if(std::max(std::abs(residual[0]),std::abs(residual[1]))<norm) {candidate=c;accepted=true;break;}
      }
      if(!accepted) throw std::runtime_error("burn_and_mix: abundance line search failed");
    }
    if(!converged) throw std::runtime_error("burn_and_mix: abundance iteration limit");
    validate_composition(candidate);
    for(std::size_t i=begin;i<end;++i) output[i]=candidate;
  }
  if(next!=thermal.size()) throw std::invalid_argument("burn_and_mix: incomplete region partition");
  return output;
}

EvolutionStep evolve_step(const Model& previous,const Physics& p,const Atmosphere& atmosphere,
    double dt,const EvolutionOptions& options) {
  EvolutionStep result;result.model=previous;
  if(!p.eos || !p.nuclear || !p.opacity || !p.eos->has_internal_energy() || !(dt>0) || !std::isfinite(dt)
      || previous.comp.size()!=previous.size()
      || !std::isfinite(previous.age+dt) || previous.age<0 || !(previous.age+dt>previous.age)
      || !std::isfinite(options.abundance_tolerance) || options.abundance_tolerance<=0
      || !std::isfinite(options.max_abundance_change) || options.max_abundance_change<=0)
    throw std::invalid_argument("evolve_step: invalid physics, age or options");
  const auto weights=nodal_mass_weights(previous);
  for(const auto& c:previous.comp) validate_composition(c);
  Model current=previous;
  try {
    auto burning=[&](const MixingRegions& regions) {
      if(regions.size()==1 || (p.alpha_semiconvection==0 && p.alpha_thermohaline==0))
        return burn_and_mix(current,previous,*p.nuclear,regions,dt,options.abundance_tolerance*.1);
      return burn_and_transport(current,previous,*p.nuclear,regions,secular_mixing_diffusivities(current,p),
          dt,options.abundance_tolerance*.1);
    };
    for(std::size_t iteration=0;iteration<options.max_coupling_iterations;++iteration) {
      const auto regions=convective_mixing_regions(current,p);
      current.comp=burning(regions);
      double change=0;
      for(std::size_t i=0;i<current.size();++i) for(std::size_t j=0;j<NSPEC;++j)
        change=std::max(change,std::abs(current.comp[i].X[j]-previous.comp[i].X[j]));
      if(change>options.max_abundance_change) throw std::runtime_error("evolve_step: abundance change exceeds step limit");
      const auto structure=relax(current,p,atmosphere,options.relaxation,dt,&previous);
      result.coupling_iterations=iteration+1;result.residual=structure.residual;result.correction=structure.correction;
      if(!structure.converged) throw std::runtime_error("evolve_step: "+structure.message);
      current=structure.model;
      const auto next_regions=convective_mixing_regions(current,p);
      const auto next_comp=burning(next_regions);
      double residual=0;
      for(std::size_t i=0;i<current.size();++i) for(std::size_t j=0;j<NSPEC;++j)
        residual=std::max(residual,std::abs(current.comp[i].X[j]-next_comp[i].X[j]));
      result.abundance_residual=residual;
      if(regions!=next_regions || residual>options.abundance_tolerance) continue;
      double mass_release=0;
      for(std::size_t i=0;i<current.size();++i) {
        const auto e=p.eos->eval(current.T(i),current.rho(i),current.comp[i]);
        const auto old=p.eos->eval(previous.T(i),previous.rho(i),previous.comp[i]);
        const auto nuclear=p.nuclear->eval(current.T(i),current.rho(i),current.comp[i]);
        result.nuclear_luminosity+=weights[i]*nuclear.eps;
        result.neutrino_luminosity+=weights[i]*nuclear.eps_neutrino;
        result.thermal_neutrino_luminosity+=weights[i]*evaluate_losses(
            p.neutrino_losses,current.T(i),current.rho(i),current.comp[i]).eps;
        result.gravitational_luminosity+=weights[i]*detail::gravitational_heating(e.E,e.P,current.rho(i),old.E,previous.rho(i),dt);
        // Subtract the conserved baryon rest energy before differencing.
        // This reduces cancellation, while retaining nuclear binding energy.
        for(std::size_t j=0;j<NSPEC;++j)
          mass_release-=weights[i]*(nuclides[j].A/mass_numbers[j]-1)
            *(current.comp[i].X[j]-previous.comp[i].X[j])*constants::c*constants::c/dt;
      }
      result.luminosity_balance=(result.nuclear_luminosity+result.gravitational_luminosity
          -result.thermal_neutrino_luminosity)/current.y.back().L-1;
      const double release=result.nuclear_luminosity+result.neutrino_luminosity;
      result.nuclear_mass_balance=release>0 ? mass_release/release-1 : 0;
      for(auto [begin,end]:regions) if(end>begin+1) {
        ++result.mixed_regions;
        for(std::size_t i=begin;i<end;++i) result.convective_mass_fraction+=weights[i]/current.M;
      }
      current.age=previous.age+dt;
      result.model=std::move(current);result.converged=true;result.message="converged";return result;
    }
    result.message="evolve_step: coupling iteration limit";
  } catch(const std::exception& e) {result.message=e.what();}
  return result;
}
} // namespace ember

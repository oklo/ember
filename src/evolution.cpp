#include "ember/evolution.hpp"
#include "ember/constants.hpp"
#include "ember/convection.hpp"
#include "ember/energy_grid.hpp"
#include "energy.hpp"
#include "thermal_transport.hpp"
#include "ember/species_transport.hpp"
#include "ember/cn_burning.hpp"
#include "ember/deuterium_burning.hpp"
#include "ember/cn_transport.hpp"
#include "ember/metal_microscopic_transport.hpp"
#include "parallel_evaluate.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace ember {
namespace {
void validate_composition(const Composition& c) {
  if(c.basis!=AbundanceBasis::baryon_mass || std::abs(c.sum()-1)>1e-10)
    throw std::invalid_argument("evolution requires normalized baryonic abundances");
  for(double x:c.X) if(!std::isfinite(x) || x<0)
    throw std::invalid_argument("evolution requires finite nonnegative abundances");
  if(c.cn_molality)(void)cn_physical_ledger(c,*c.cn_molality);
}
double composition_difference(const Composition& a,const Composition& b) {
  if(a.cn_mass_convention!=b.cn_mass_convention)
    throw std::logic_error("evolution changed the metal-mass convention");
  if(a.cn_molality.has_value()!=b.cn_molality.has_value())
    throw std::logic_error("evolution changed the active CN inventory");
  double difference=0;
  for(std::size_t j=0;j<NSPEC;++j)difference=std::max(difference,std::abs(a.X[j]-b.X[j]));
  if(a.cn_molality)for(std::size_t j=0;j<3;++j)
    difference=std::max(difference,mass_numbers[j+3]*std::abs((*a.cn_molality)[j]-(*b.cn_molality)[j]));
  return difference;
}
// The total species rates depend on sources across a mixed region. Hold them
// fixed for a local thermal Newton solve and converge them in the outer loop.
class FrozenSpeciesHeat final:public MicroscopicTransport {
 public:
  explicit FrozenSpeciesHeat(const MicroscopicTransport& source):source_(source) {}
  std::span<const SpeciesVector> rates;
  std::span<const MetalSpeciesVector> metal_rates;
  MicroscopicFaceResponse eval(std::size_t face,double mlo,double mhi,const Point& lo,
      const Composition& a,const Point& hi,const Composition& b,bool derivatives)const override {
    return source_.eval(face,mlo,mhi,lo,a,hi,b,derivatives);
  }
  MicroscopicHeatResponse heat(std::size_t face,double mlo,double mhi,const Point& lo,
      const Composition& a,const Point& hi,const Composition& b,bool derivatives)const override {
    if(!metal_rates.empty()) {
      const auto* metal=dynamic_cast<const MetalMicroscopicTransport*>(&source_);
      if(!metal || face>=metal_rates.size())throw std::logic_error("evolve_step: missing metal heat provider or rate");
      return microscopic_heat_with_total_metal_rate(*metal,face,mlo,mhi,lo,a,hi,b,metal_rates[face],derivatives);
    }
    if(rates.empty())return source_.heat(face,mlo,mhi,lo,a,hi,b,derivatives);
    if(face>=rates.size())throw std::out_of_range("evolve_step: missing total species rate");
    return microscopic_heat_with_total_species_rate(source_,face,mlo,mhi,lo,a,hi,b,rates[face],derivatives);
  }
  bool requires_positive_species_guess()const override {return source_.requires_positive_species_guess();}
  const char* name()const override {return "material heat at prescribed total species rates";}
 private:
  const MicroscopicTransport& source_;
};
struct CompositionUpdate {
  std::vector<Composition> composition;
  std::vector<SpeciesVector> total_rates;
  std::vector<MetalSpeciesVector> total_metal_rates;
};
std::vector<BurningResponse> local_burning_response(const Model& thermal,
    const Model& previous,const Nuclear& nuclear,const MixingRegions& regions,
    double dt,double tolerance) {
  // Transport is deliberately absent from this approximate Jacobian. Each
  // singleton has an independent local burn; mixed regions receive no slope.
  if(std::none_of(regions.begin(),regions.end(),[](auto r){return r.second==r.first+1;}))return {};
  constexpr double h=1e-5;
  try {
    const auto base=burn_and_mix(thermal,previous,nuclear,regions,dt,tolerance);
    auto hot=thermal,dense=thermal;
    for(auto& y:hot.y)y.lnT+=h;
    for(auto& y:dense.y)y.lnrho+=h;
    const auto warm=burn_and_mix(hot,previous,nuclear,regions,dt,tolerance);
    const auto packed=burn_and_mix(dense,previous,nuclear,regions,dt,tolerance);
    std::vector<BurningResponse> response(thermal.size());
    for(const auto [begin,end]:regions)if(end==begin+1) {
      const double T=thermal.T(begin),rho=thermal.rho(begin);
      const double e0=nuclear.eval(T,rho,base[begin]).eps;
      auto& r=response[begin];
      r.dEps_dlnT=(nuclear.eval(T,rho,warm[begin]).eps-e0)/h;
      r.dEps_dlnRho=(nuclear.eval(T,rho,packed[begin]).eps-e0)/h;
      if(!std::isfinite(r.dEps_dlnT) || !std::isfinite(r.dEps_dlnRho))return {};
    }
    return response;
  } catch(const std::domain_error&) {
    // A derivative probe may leave the nuclear source domain. The original
    // solver remains available and still enforces the physical domain.
    return {};
  } catch(const std::runtime_error&) {
    // Failure of an auxiliary local burn need not reject the full coupled
    // solve. No composition or heating from these probes is retained.
    return {};
  }
}
}
std::vector<double> nodal_mass_weights(const Model& m) {
  if(m.size()<2 || m.m.size()!=m.size() || !(m.m[0]>0) || !m.valid_outer_mass())
    throw std::invalid_argument("nodal_mass_weights: invalid mass mesh");
  std::vector<double> w(m.size());w[0]=m.m[0];
  for(std::size_t i=1;i<m.size();++i) {
    const double dm=m.m[i]-m.m[i-1];
    if(!(dm>0) || !std::isfinite(dm)) throw std::invalid_argument("nodal_mass_weights: unordered mesh");
    w[i-1]+=.5*dm;w[i]+=.5*dm;
  }
  w.back()+=m.envelope_mass;
  return w;
}

MixingRegions schwarzschild_mixing_regions(const Model& m,const Physics& p) {
  auto selected=p;selected.criterion=ConvectiveCriterion::schwarzschild;
  return convective_mixing_regions(m,selected);
}

MixingRegions convective_mixing_regions(const Model& m,const Physics& p,std::size_t threads) {
  if(!p.eos || !p.opacity || m.comp.size()!=m.size())
    throw std::invalid_argument("mixing regions: invalid model or physics");
  detail::check_thermal_transport(p);
  nodal_mass_weights(m);
  struct Local {double P,kappa,ad,delta;};
  std::vector<Local> local(m.size());
  detail::independent_evaluations(m.size(),threads,[&](std::size_t i) {
    const auto e=p.eos->eval(m.T(i),m.rho(i),m.comp[i]);
    local[i]={e.P,p.opacity->eval(m.T(i),m.rho(i),m.comp[i]).kappa,e.grad_ad,e.delta};
  });
  std::vector<char> boundary(m.size()-1);
  detail::independent_evaluations(m.size()-1,threads,[&](std::size_t i) {
    // Same arithmetic midpoint transport quantities as zone assembly.
    const double T=.5*(m.T(i)+m.T(i+1)),mass=.5*(m.m[i]+m.m[i+1]);
    const double P=.5*(local[i].P+local[i+1].P),k=.5*(local[i].kappa+local[i+1].kappa);
    const double L=thermal_face_luminosity(m,i),ad=.5*(local[i].ad+local[i+1].ad);
    MicroscopicHeatResponse heat;
    if(p.microscopic) {
      heat=microscopic_heat(*p.microscopic,i,m.m[i],m.m[i+1],m.y[i],m.comp[i],m.y[i+1],m.comp[i+1],false);
    }
    const double rad=detail::thermal_transport<0>(T,.5*(m.rho(i)+m.rho(i+1)),P,mass,k,L,
        p.microscopic!=nullptr || face_luminosities(m),heat.carried_luminosity,heat.conductivity,
        m.y[i].lnT,m.y[i+1].lnT).gradient.value;
    double B=0;
    if(p.criterion==ConvectiveCriterion::ledoux)
      B=composition_buoyancy(*p.eos,T,P,.5*(local[i].delta+local[i+1].delta),
          std::log(local[i+1].P)-std::log(local[i].P),m.comp[i],m.comp[i+1],.5*(m.rho(i)+m.rho(i+1))).B;
    boundary[i]=!(rad>ad+B);
  });
  MixingRegions regions;std::size_t first=0;
  for(std::size_t i=0;i+1<m.size();++i)if(boundary[i]) {regions.emplace_back(first,i+1);first=i+1;}
  regions.emplace_back(first,m.size());return regions;
}

double abundance_change_after_mixing(const Model& previous,
    const std::vector<Composition>& next,const MixingRegions& regions) {
  if(previous.comp.size()!=previous.size() || next.size()!=previous.size())
    throw std::invalid_argument("abundance change: inconsistent mesh size");
  const auto weights=nodal_mass_weights(previous);
  double largest=0;std::size_t consumed=0;
  for(auto [begin,end]:regions) {
    if(begin!=consumed || end<=begin || end>previous.size())
      throw std::invalid_argument("abundance change: regions must partition the mesh");
    consumed=end;
    auto average=previous.comp[begin];
    if(end>begin+1) {
      long double mass=0;std::array<long double,NSPEC> inventory{};
      std::array<long double,3> cn{};
      for(std::size_t i=begin;i<end;++i) {
        const auto& c=previous.comp[i];
        if(c.basis!=average.basis || c.metal_inventory!=average.metal_inventory)
          throw std::invalid_argument("abundance change: incompatible composition conventions");
        (void)composition_difference(c,average);
        mass+=weights[i];
        for(std::size_t j=0;j<NSPEC;++j)inventory[j]+=static_cast<long double>(weights[i])*c.X[j];
        if(c.cn_molality)for(std::size_t j=0;j<3;++j)
          cn[j]+=static_cast<long double>(weights[i])*(*c.cn_molality)[j];
      }
      for(std::size_t j=0;j<NSPEC;++j)average.X[j]=static_cast<double>(inventory[j]/mass);
      if(average.cn_molality)for(std::size_t j=0;j<3;++j)
        (*average.cn_molality)[j]=static_cast<double>(cn[j]/mass);
    }
    for(std::size_t i=begin;i<end;++i)
      largest=std::max(largest,composition_difference(next[i],average));
  }
  if(consumed!=previous.size())throw std::invalid_argument("abundance change: incomplete partition");
  return largest;
}

std::vector<Composition> burn_and_mix(const Model& thermal,const Model& previous,
    const Nuclear& nuclear,const MixingRegions& regions,double dt,double tolerance) {
  if(!std::isfinite(dt) || dt<=0 || !std::isfinite(tolerance) || tolerance<=0
      || thermal.m!=previous.m || thermal.size()!=previous.size()
      || thermal.comp.size()!=thermal.size() || previous.comp.size()!=previous.size())
    throw std::invalid_argument("burn_and_mix: invalid step, tolerance or previous model");
  const auto weights=nodal_mass_weights(thermal);
  for(const auto& c:previous.comp) validate_composition(c);
  if(const auto* network=dynamic_cast<const PPDeuterium*>(&nuclear))
    return burn_deuterium_and_mix(thermal,previous,*network,regions,dt,tolerance);
  if(const auto* network=dynamic_cast<const PPCNNetwork*>(&nuclear);
      network && previous.comp.front().cn_mass_convention==CNMassConvention::explicit_metal_mass)
    return burn_metal_cn_and_transport(thermal,previous,*network,regions,{},dt,tolerance);
  for(const auto& c:previous.comp)if(c[Species::H2]!=0)
    throw std::invalid_argument("burn_and_mix: explicit D requires its coupled network");
  if(const auto* network=dynamic_cast<const PPCNNetwork*>(&nuclear)) {
    std::vector<CNAbundances> old_cn;old_cn.reserve(previous.size());
    for(const auto& c:previous.comp) {
      if(!c.cn_molality)throw std::invalid_argument("burn_and_mix: CN inventory missing");
      old_cn.push_back(*c.cn_molality);
    }
    return burn_cn_and_mix(thermal,previous,old_cn,network->pp(),network->cn(),regions,dt,tolerance).lookup;
  }
  for(const auto& c:previous.comp)if(c.cn_molality)
    throw std::invalid_argument("burn_and_mix: explicit CN requires its coupled network");
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

MixingRegions instantaneous_mixing_regions(const Model& model,const Physics& physics,
    const EvolutionOptions& options) {
  if(!std::isfinite(options.instantaneous_mixing_below_T) || options.instantaneous_mixing_below_T<0)
    throw std::invalid_argument("invalid instantaneous mixing temperature");
  if(options.convective_mixing!=ConvectiveMixing::instantaneous
      && options.convective_mixing!=ConvectiveMixing::finite_implicit
      && options.convective_mixing!=ConvectiveMixing::finite_lagged)
    throw std::invalid_argument("unknown convective mixing choice");
  if(options.convective_mixing==ConvectiveMixing::instantaneous)
    return convective_mixing_regions(model,physics,options.relaxation.zone_threads);
  MixingRegions regions;regions.reserve(model.size());
  if(options.instantaneous_mixing_below_T==0) {
    for(std::size_t i=0;i<model.size();++i)regions.emplace_back(i,i+1);
    return regions;
  }
  for(const auto [begin,end]:convective_mixing_regions(model,physics,options.relaxation.zone_threads)) {
    auto first=begin;
    for(std::size_t i=begin;i+1<end;++i)
      if(std::min(model.T(i),model.T(i+1))>=options.instantaneous_mixing_below_T) {
        regions.emplace_back(first,i+1);first=i+1;
      }
    regions.emplace_back(first,end);
  }
  return regions;
}

EvolutionStep evolve_step(const Model& previous,const Physics& p,const Atmosphere& atmosphere,
    double dt,const EvolutionOptions& options) {
  EvolutionStep result;result.model=previous;
  if(!p.eos || !p.nuclear || !p.opacity || !p.eos->has_internal_energy() || !(dt>0) || !std::isfinite(dt)
      || previous.comp.size()!=previous.size() || !p.burning_response.empty()
      || !std::isfinite(previous.age+dt) || previous.age<0 || !(previous.age+dt>previous.age)
      || !std::isfinite(options.abundance_tolerance) || options.abundance_tolerance<=0
      || !std::isfinite(options.homogeneous_abundance_tolerance) || options.homogeneous_abundance_tolerance<0
      || options.homogeneous_abundance_tolerance>options.abundance_tolerance
      || !std::isfinite(options.instantaneous_mixing_below_T) || options.instantaneous_mixing_below_T<0
      || (options.convective_mixing==ConvectiveMixing::instantaneous && options.instantaneous_mixing_below_T!=0)
      || !std::isfinite(options.coupling_stop_tolerance) || options.coupling_stop_tolerance<0
      || !std::isfinite(options.verification_residual_tolerance) || options.verification_residual_tolerance<0
      || !std::isfinite(options.verification_correction_tolerance) || options.verification_correction_tolerance<0
      || !std::isfinite(options.material_heat_tolerance) || options.material_heat_tolerance<=0
      || !std::isfinite(options.max_abundance_change) || options.max_abundance_change<=0)
    throw std::invalid_argument("evolve_step: invalid physics, age or options");
  const auto weights=nodal_mass_weights(previous);
  detail::check_thermal_transport(p);
  for(const auto& c:previous.comp) validate_composition(c);
  Model current=previous;
  if(const auto* guess=options.initial_structure_guess) {
    if(guess->M!=previous.M || guess->m!=previous.m || guess->size()!=previous.size()
        || guess->luminosity_grid!=previous.luminosity_grid)
      throw std::invalid_argument("evolve_step: initial structure uses a different mesh");
    for(const auto& point:guess->y)
      for(std::size_t k=0;k<NVAR;++k)
        if(!std::isfinite(point[static_cast<Var>(k)]))
          throw std::invalid_argument("evolve_step: nonfinite initial structure");
    current.y=guess->y;
  }
  try {
    for(const auto& c:previous.comp)if(c[Species::H2]!=0 && p.microscopic)
      throw std::invalid_argument("evolve_step: microscopic transport with D needs an isotope-aware provider");
    const auto* metal=dynamic_cast<const MetalMicroscopicTransport*>(p.microscopic);
    const auto* cn_network=dynamic_cast<const PPCNNetwork*>(p.nuclear);
    for(const auto& c:previous.comp) {
      const bool explicit_metals=c.cn_mass_convention==CNMassConvention::explicit_metal_mass;
      if((metal && !explicit_metals) || (explicit_metals && (!cn_network || (p.microscopic && !metal)
          || (!p.microscopic && !p.explicit_metal_mixing_only))))
        throw std::invalid_argument("evolve_step: physical metal abundances require their CN network and compatible transport provider");
    }
    if(metal && !cn_network)
      throw std::invalid_argument("evolve_step: physical metal transport requires the explicit CN network");
    Physics coupled=p;
    std::optional<FrozenSpeciesHeat> frozen;
    if(p.microscopic && p.microscopic->uses_total_species_heat()) {
      frozen.emplace(*p.microscopic);coupled.microscopic=&*frozen;
    }
    const bool finite=options.convective_mixing!=ConvectiveMixing::instantaneous;
    if(options.convective_mixing!=ConvectiveMixing::instantaneous
        && options.convective_mixing!=ConvectiveMixing::finite_implicit
        && options.convective_mixing!=ConvectiveMixing::finite_lagged)
      throw std::invalid_argument("evolve_step: unknown convective mixing choice");
    if(finite && (!metal || !cn_network || !frozen || !face_luminosities(previous)))
      throw std::invalid_argument("evolve_step: finite convection requires physical metal CN transport, total-species heat and volume-face luminosities");
    auto check_rates=[&](auto rates) {
      if(!rates.empty() && rates.size()+1!=previous.size())
        throw std::invalid_argument("evolve_step: previous species heat rates have wrong dimensions");
      for(const auto& face:rates)for(double rate:face)if(!std::isfinite(rate))
        throw std::invalid_argument("evolve_step: nonfinite previous species heat rate");
    };
    check_rates(options.previous_species_heat_rates);check_rates(options.previous_metal_heat_rates);
    if((!options.previous_species_heat_rates.empty() && (!frozen || metal))
        || (!options.previous_metal_heat_rates.empty() && (!frozen || !metal)))
      throw std::invalid_argument("evolve_step: previous heat rates do not match the material provider");
    if(frozen) {
      frozen->rates=options.previous_species_heat_rates;
      frozen->metal_rates=options.previous_metal_heat_rates;
    }
    std::vector<double> lagged_mixing;
    if(options.convective_mixing==ConvectiveMixing::finite_lagged) {
      const bool supplied=metal?!options.previous_metal_heat_rates.empty():!options.previous_species_heat_rates.empty();
      if(!supplied)for(const auto& c:previous.comp)
        if(composition_difference(c,previous.comp.front())!=0)
          throw std::invalid_argument("evolve_step: lagged convection needs the previous total species heat rates for a stratified model");
      // Freeze the whole mass conductance, including radius and density.
      // Freezing D alone changes the coefficient during the thermal iteration.
      lagged_mixing=finite_mixing_conductances(previous,coupled);
    }
    auto mixing_regions=[&]() {
      return instantaneous_mixing_regions(current,coupled,options);
    };
    auto abundance_tolerance=[&](const MixingRegions& regions) {
      return regions.size()==1 && options.homogeneous_abundance_tolerance>0
          ?options.homogeneous_abundance_tolerance:options.abundance_tolerance;
    };
    auto burning=[&](const MixingRegions& regions)->CompositionUpdate {
      const double tolerance=abundance_tolerance(regions);
      if(p.microscopic) {
        std::vector<double> mixing;
        if(!lagged_mixing.empty())mixing=lagged_mixing;
        else if(finite)mixing=finite_mixing_conductances(current,coupled);
        else {
          mixing=secular_mixing_diffusivities(current,coupled);
          for(std::size_t i=0;i<mixing.size();++i) {
            const double r=.5*(current.r(i)+current.r(i+1)),rho=.5*(current.rho(i)+current.rho(i+1));
            const double area_mass=4*M_PI*r*r*rho;
            mixing[i]*=area_mass*area_mass/(current.m[i+1]-current.m[i]);
          }
        }
        auto base_flux=[&](std::size_t i,const Composition& left,const Composition& right,bool derivatives) {
          return microscopic_face(*p.microscopic,i,current.m[i],current.m[i+1],
              current.y[i],left,current.y[i+1],right,derivatives).species;
        };
        // A local correction can reach its floating-point floor before the
        // integrated balance does. Keep conservation and outer coupling tight
        // without requiring this inner Newton solve to resolve that floor.
        SpeciesTransportOptions transport_options;transport_options.abundance_tolerance=10*tolerance;
        transport_options.integrated_balance_tolerance=tolerance*.1;
        transport_options.seed_present_species=p.microscopic->requires_positive_species_guess();
        if(metal) {
          // Reuse the last coupling iterate; previous still supplies every
          // storage term and conservation reference in the implicit solve.
          transport_options.initial_guess=current.comp;
          transport_options.evaluation_threads=options.relaxation.zone_threads;
          MetalCNFlux flux=[&](std::size_t i,const Composition& left,const Composition& right,bool derivatives) {
            const auto face=metal_species_face(*metal,i,current.m[i],current.m[i+1],
                current.y[i],left,current.y[i+1],right,derivatives);
            auto combined=common_metal_cn_flux(face,left,right,derivatives);
            metal->add_mixing_flux(combined,i,current.m[i],current.m[i+1],
                current.y[i],left,current.y[i+1],right,derivatives);
            return combined;
          };
          auto full=burn_metal_cn_and_diffuse(current,previous,*cn_network,regions,flux,dt,transport_options,mixing);
          Model updated=current;updated.comp=std::move(full.composition);
          auto redistribution=reconstruct_metal_fluxes(updated,previous,*cn_network,regions,full.boundary_fluxes,dt,options.relaxation.zone_threads);
          auto continuous=[&]() {
            for(const auto& cell:redistribution.cell_balances)for(double balance:cell)
              if(std::abs(balance)>transport_options.integrated_balance_tolerance)return false;
            return true;
          };
          if(!continuous()) {
            // A small Newton correction and a closed global inventory can
            // coexist with a larger regional flux residual. Refine species
            // from this answer before discarding the entire stellar step.
            // The independent continuity and energy requirements are unchanged.
            transport_options.initial_guess=updated.comp;
            transport_options.abundance_tolerance=std::min(
                transport_options.abundance_tolerance,transport_options.integrated_balance_tolerance);
            full=burn_metal_cn_and_diffuse(current,previous,*cn_network,regions,flux,dt,transport_options,mixing);
            updated.comp=std::move(full.composition);
            redistribution=reconstruct_metal_fluxes(updated,previous,*cn_network,regions,
                full.boundary_fluxes,dt,options.relaxation.zone_threads);
            if(!continuous())
              throw std::runtime_error("evolve_step: reconstructed metal species continuity exceeds tolerance after refinement");
          }
          return {std::move(updated.comp),{},std::move(redistribution.face_rates)};
        }
        SpeciesTransportResult species;
        if(const auto* network=dynamic_cast<const PPCNNetwork*>(p.nuclear)) {
          transport_options.evaluation_threads=options.relaxation.zone_threads;
          if(p.cn_microscopic==CNMicroscopicApproximation::unselected)
            throw std::invalid_argument("evolve_step: explicit CN microscopic approximation required");
          CNSpeciesFlux flux=[&](std::size_t i,const Composition& left,const Composition& right,bool derivatives) {
            return trace_cn_flux(base_flux(i,left,right,derivatives),left,right,p.cn_microscopic,derivatives);
          };
          auto full=burn_cn_and_diffuse(current,previous,*network,regions,flux,dt,transport_options,mixing);
          species.composition=std::move(full.composition);
          for(const auto& face:full.boundary_fluxes)
            species.boundary_fluxes.push_back({face.face,{face.rate[0],face.rate[1]}});
          // Fixed-GS98 material thermodynamics depends on the H/He lookup
          // coordinates. Its carried heat uses their conservative rates;
          // CN-dependent material enthalpy/collision feedback is not included
          // in this explicitly selected trace approximation.
        } else {
          SpeciesFlux flux=[&](std::size_t i,const Composition& left,const Composition& right,bool derivatives) {
            auto f=base_flux(i,left,right,derivatives);const double g=mixing[i];
            for(std::size_t k=0;k<2;++k) {
              f.rate[k]+=g*(left.X[k]-right.X[k]);
              if(derivatives){f.dleft[k][k]+=g;f.dright[k][k]-=g;}
            }
            return f;
          };
          species=burn_and_diffuse(current,previous,*p.nuclear,regions,flux,dt,transport_options);
        }
        if(!frozen)return {std::move(species.composition),{},{}};
        Model updated=current;updated.comp=std::move(species.composition);
        auto redistribution=reconstruct_species_fluxes(updated,previous,*p.nuclear,regions,species.boundary_fluxes,dt);
        for(const auto& cell:redistribution.cell_balances)for(double balance:cell)
          if(std::abs(balance)>transport_options.integrated_balance_tolerance)
            throw std::runtime_error("evolve_step: reconstructed species continuity exceeds tolerance");
        return {std::move(updated.comp),std::move(redistribution.face_rates),{}};
      }
      if(regions.size()==1 || (p.alpha_semiconvection==0 && p.alpha_thermohaline==0))
        return {burn_and_mix(current,previous,*p.nuclear,regions,dt,tolerance*.1),{},{}};
      return {burn_and_transport(current,previous,*p.nuclear,regions,secular_mixing_diffusivities(current,coupled),
          dt,tolerance*.1),{},{}};
    };
    MixingRegions pending_regions;
    CompositionUpdate pending;
    bool have_pending=false;
    std::string coupling_reason;
    double heat_residual_q=0,heat_residual_surface=0;
    std::vector<BurningResponse> burning_response;
    for(std::size_t iteration=0;iteration<options.max_coupling_iterations;++iteration) {
      // The convergence check below already evaluates the next composition
      // on this exact thermal state. Retain it when another coupling pass is
      // needed; no state or convection boundary changes between these calls.
      const auto regions=have_pending?std::move(pending_regions):mixing_regions();
      auto update=have_pending?std::move(pending):burning(regions);
      current.comp=std::move(update.composition);
      if(frozen) {frozen->rates=update.total_rates;frozen->metal_rates=update.total_metal_rates;}
      have_pending=false;
      double change=0;
      if(options.abundance_cap_after_mixing && !finite)
        change=abundance_change_after_mixing(previous,current.comp,regions);
      else for(std::size_t i=0;i<current.size();++i)
        change=std::max(change,composition_difference(current.comp[i],previous.comp[i]));
      if(change>options.max_abundance_change) throw std::runtime_error("evolve_step: abundance change exceeds step limit");
      if(options.linearized_burning && iteration==0)
        burning_response=local_burning_response(current,previous,*p.nuclear,regions,
            dt,.1*abundance_tolerance(regions));
      auto responsive=coupled;
      if(!burning_response.empty()) {
        for(const auto [begin,end]:regions)for(std::size_t i=begin;i<end;++i) {
          auto& r=burning_response[i];
          if(end!=begin+1)r.dEps_dlnT=r.dEps_dlnRho=0;
          r.lnT_ref=current.y[i].lnT;r.lnRho_ref=current.y[i].lnrho;
        }
        responsive.burning_response=burning_response;
      }
      const auto structure=relax(current,responsive,atmosphere,options.relaxation,dt,&previous);
      result.coupling_iterations=iteration+1;result.residual=structure.residual;result.correction=structure.correction;
      if(!structure.converged) throw std::runtime_error("evolve_step: "+structure.message);
      current=structure.model;
      auto next_regions=mixing_regions();
      auto next=burning(next_regions);
      double residual=0;
      for(std::size_t i=0;i<current.size();++i)
        residual=std::max(residual,composition_difference(current.comp[i],next.composition[i]));
      result.abundance_residual=residual;
      result.material_heat_residual=0;
      heat_residual_q=heat_residual_surface=0;
      if(frozen) {
        double luminosity_floor=0;
        for(const auto& point:current.y)luminosity_floor=std::max(luminosity_floor,
            minimum_luminosity_scale_fraction*std::abs(point.L));
        std::vector<double> face_residual(current.size()-1),absolute_heat_change(current.size()-1);
        detail::independent_evaluations(current.size()-1,options.relaxation.zone_threads,[&](std::size_t i) {
          const auto old_heat=microscopic_heat(*frozen,i,current.m[i],current.m[i+1],
              current.y[i],current.comp[i],current.y[i+1],current.comp[i+1],false);
          // The supplied rate enters linearly. Reuse its exact enthalpy
          // derivative instead of repeating the entire kinetic evaluation.
          double heat_change=0;
          for(std::size_t k=0;k<2;++k)
            heat_change+=old_heat.total_rate_enthalpy[k]*(metal?
                next.total_metal_rates[i][k]-update.total_metal_rates[i][k]:
                next.total_rates[i][k]-update.total_rates[i][k]);
          if(metal)heat_change+=old_heat.total_metal_rate_enthalpy*
              (next.total_metal_rates[i][2]-update.total_metal_rates[i][2]);
          const double scale=std::max({std::abs(current.y[i].L),std::abs(current.y[i+1].L),
              luminosity_floor,std::numeric_limits<double>::min()});
          face_residual[i]=std::abs(heat_change)/scale;
          absolute_heat_change[i]=std::abs(heat_change);
        });
        for(std::size_t i=0;i<face_residual.size();++i) {
          if(face_residual[i]>result.material_heat_residual) {
            result.material_heat_residual=face_residual[i];
            heat_residual_q=.5*(current.m[i]+current.m[i+1])/current.M;
          }
          heat_residual_surface=std::max(heat_residual_surface,absolute_heat_change[i]/
              std::max(std::abs(current.y.back().L),std::numeric_limits<double>::min()));
        }
      }
      bool deuterium_coupled=true;
      if(dynamic_cast<const PPDeuterium*>(p.nuclear) || cn_network) {
        for(std::size_t i=0;i<current.size();++i) {
          const double a=current.comp[i][Species::H2],b=next.composition[i][Species::H2];
          deuterium_coupled &= std::abs(a-b)<=1e-8*std::max(a,b)+2*std::numeric_limits<double>::denorm_min();
        }
      }
      const double stop=std::max(options.coupling_stop_tolerance,abundance_tolerance(next_regions));
      if(!deuterium_coupled || regions!=next_regions || residual>stop
          || result.material_heat_residual>options.material_heat_tolerance) {
        coupling_reason=regions!=next_regions?"convective boundary changes":
          (!deuterium_coupled?"deuterium coupling":(residual>stop?"composition coupling":"transported heat coupling"));
        pending_regions=std::move(next_regions);pending=std::move(next);
        have_pending=true;continue;
      }
      // The last species solve uses the final thermal state. Returning the
      // preceding composition can violate its integrated source balance even
      // when the outer iteration's abundance change is small. Retain this
      // solution and its matching heat fluxes only if the actual returned
      // structure and convective partition still satisfy convergence.
      const auto final_regions=finite?convective_mixing_regions(current,coupled,options.relaxation.zone_threads):next_regions;
      if(current.comp!=next.composition || update.total_rates!=next.total_rates
          || update.total_metal_rates!=next.total_metal_rates || !burning_response.empty()) {
        current.comp=std::move(next.composition);
        if(frozen) {frozen->rates=next.total_rates;frozen->metal_rates=next.total_metal_rates;}
        double final_change=0;
        if(options.abundance_cap_after_mixing && !finite)
          final_change=abundance_change_after_mixing(previous,current.comp,next_regions);
        else for(std::size_t i=0;i<current.size();++i)
          final_change=std::max(final_change,composition_difference(current.comp[i],previous.comp[i]));
        if(final_change>options.max_abundance_change)
          throw std::runtime_error("evolve_step: abundance change exceeds step limit");
        auto verification=options.relaxation;verification.max_iterations=0;
        if(options.verification_residual_tolerance>0)verification.residual_tolerance=options.verification_residual_tolerance;
        if(options.verification_correction_tolerance>0)verification.correction_tolerance=options.verification_correction_tolerance;
        const auto checked=relax(current,coupled,atmosphere,verification,dt,&previous);
        result.residual=checked.residual;result.correction=checked.correction;
        if(!checked.converged || convective_mixing_regions(current,coupled,options.relaxation.zone_threads)!=final_regions
            || (finite && mixing_regions()!=next_regions)) {
          coupling_reason=checked.converged?"final convective boundary changes":"final structure verification";
          continue;
        }
        update.total_rates=std::move(next.total_rates);
        update.total_metal_rates=std::move(next.total_metal_rates);
      }
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
        mass_release-=weights[i]*(p.nuclear->rest_energy_correction(current.comp[i])
          -p.nuclear->rest_energy_correction(previous.comp[i]))/dt;
      }
      result.luminosity_balance=(result.nuclear_luminosity+result.gravitational_luminosity
          -result.thermal_neutrino_luminosity)/current.y.back().L-1;
      const double release=result.nuclear_luminosity+result.neutrino_luminosity;
      result.nuclear_mass_balance=release>0 ? mass_release/release-1 : 0;
      const auto physical_regions=finite?convective_mixing_regions(current,coupled,options.relaxation.zone_threads):regions;
      for(auto [begin,end]:physical_regions) if(end>begin+1) {
        ++result.mixed_regions;
        for(std::size_t i=begin;i<end;++i) result.convective_mass_fraction+=weights[i]/current.M;
      }
      if(finite) {
        const auto faces=convective_mixing_faces(current,coupled);
        for(std::size_t i=0;i<faces.size();++i)if(faces[i].buoyancy_contrast>0) {
          const double ratio=std::abs(faces[i].composition_term)/faces[i].buoyancy_contrast;
          if(ratio>result.convection_composition_ratio) {
            result.convection_composition_ratio=ratio;result.convection_composition_face=i;
          }
        }
      }
      current.age=previous.age+dt;
      result.total_species_rates=std::move(update.total_rates);
      result.total_metal_species_rates=std::move(update.total_metal_rates);
      result.model=std::move(current);result.converged=true;result.message="converged";return result;
    }
    std::ostringstream reason;reason.precision(4);
    reason<<"evolve_step: coupling iteration limit ("<<coupling_reason
      <<"; abundance="<<result.abundance_residual<<", transported heat="<<result.material_heat_residual
      <<", heat q="<<heat_residual_q<<", max heat/surface L="<<heat_residual_surface
      <<", structure="<<result.residual<<", correction="<<result.correction<<")";
    result.message=reason.str();
  } catch(const std::exception& e) {result.message=e.what();}
  return result;
}
} // namespace ember

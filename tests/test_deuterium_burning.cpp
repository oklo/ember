#include "ember/deuterium_burning.hpp"
#include "ember/constants.hpp"
#include "../apps/evolution_checkpoint.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <numeric>

using namespace ember;
int main() {
  int failures=0;
  auto check=[&](bool ok,const char* label,double value=0.) {
    failures+=!ok;std::printf("[%s] %s %.4g\n",ok?"PASS":"FAIL",label,value);
  };
  constexpr auto d=static_cast<std::size_t>(Species::H2);
  auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  const double ions=c.mu_ions_inv(),electrons=c.mu_elec_inv();
  c.X[0]-=2e-5;c.X[d]=2e-5;
  check(std::abs(c.sum()-1)<1e-15 && std::abs(c.Z()-.02)<1e-15,"D is conserved hydrogen isotope mass, not metal mass");
  check(std::abs(c.mu_ions_inv()-ions+1e-5)<1e-15 && std::abs(c.mu_elec_inv()-electrons+1e-5)<1e-15,
        "GS98 charge and ion counts include explicit D");
  PPDeuterium network;PPChains pp(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn);
  const auto response=network.composition_response(1e6,1,c);
  auto without=c;without.X[d]=0;
  const auto capture=deuterium_capture(1e6,1,without,c.X[d]);const auto ordinary=pp.eval(1e6,1,c);
  check(std::abs((response.state.eps-ordinary.eps)/capture.source.eps-1)<1e-14,
        "explicit source adds only initial-D capture to existing pp heat");
  double derivative=0;
  for(std::size_t j=0;j<NSPEC;++j) {
    constexpr double h=1e-8;auto plus=c,minus=c;plus.X[j]+=h;
    const bool edge=c.X[j]<h;minus.X[j]+=edge?2*h:-h;
    const auto a=network.eval(1e6,1,plus),b=network.eval(1e6,1,minus);
    auto difference=[&](double x,double y,double base){return edge?(-3*base+4*x-y)/(2*h):(x-y)/(2*h);};
    derivative=std::max(derivative,std::abs(difference(a.eps,b.eps,response.state.eps)-response.deps_dX[j])/
        std::max(response.state.eps,std::abs(response.deps_dX[j])));
    for(auto i:{0UL,1UL,d})derivative=std::max(derivative,
        std::abs(difference(a.dXdt[i],b.dXdt[i],response.state.dXdt[i])-response.d_dXdt_dX[i][j])/
        std::max(std::abs(response.state.dXdt[1]),std::abs(response.d_dXdt_dX[i][j])));
  }
  check(derivative<1e-6,"all isotope/composition derivatives agree with independent perturbations",derivative);
  Model m;m.M=10;m.m={1,2,4,10};
  for(int i=0;i<4;++i) {
    auto q=c;q.X[d]+=i*2e-6;q.X[0]-=i*2e-6;
    m.comp.push_back(q);m.y.push_back({std::log(1.+i),std::log(1.+.1*i),std::log(1e6),1.});
  }
  const auto weights=nodal_mass_weights(m);auto cold=m;
  for(auto& y:cold.y)y.lnT=std::log(3000.);
  const auto mixed=burn_and_mix(cold,m,network,{{0,2},{2,4}},1e12);
  double conserved=0;
  for(std::size_t j=0;j<NSPEC;++j) {
    double delta=0;for(std::size_t i=0;i<m.size();++i)delta+=weights[i]*(mixed[i].X[j]-m.comp[i].X[j]);
    conserved=std::max(conserved,std::abs(delta)/m.M);
  }
  check(conserved<1e-15 && mixed[0].X==mixed[1].X && mixed[2].X==mixed[3].X && mixed[0].X!=mixed[2].X,
        "cold convection mixes every isotope and preserves separate regions",conserved);
  // Repeated convection in an already homogeneous cold star must not create
  // fictitious fuel or inert species through the weighted-mean arithmetic.
  auto uniform=cold;uniform.m.clear();uniform.y.clear();uniform.comp.clear();
  for(int i=0;i<512;++i) {
    uniform.m.push_back(uniform.M*std::pow((i+1)/512.,3));
    uniform.y.push_back({std::log(1.+i),0.,std::log(3000.),1.});
    uniform.comp.push_back(c);
  }
  const auto unchanged=burn_and_mix(uniform,uniform,network,{{0,512}},1.,1e-14);
  check(std::all_of(unchanged.begin(),unchanged.end(),[&](const auto& q){return q.X==c.X;}),
        "cold convection preserves a uniform isotope inventory exactly on an unequal mesh");
  double worst_residual=0,worst_heat=0;
  double initial_D=0;for(std::size_t i=0;i<m.size();++i)initial_D+=weights[i]*m.comp[i].X[d]/m.M;
  for(double dt:{1e8,1e11,1e14}) {
    const auto next=burn_and_mix(m,m,network,{{0,4}},dt,1e-14);
    double heat=0,rest=0;std::array<double,NSPEC> residual{};
    for(std::size_t i=0;i<m.size();++i) {
      const auto s=network.eval(m.T(i),m.rho(i),next[i]);heat+=weights[i]*(s.eps+s.eps_neutrino);
      for(std::size_t j=0;j<NSPEC;++j) {
        const double delta=next[i].X[j]-m.comp[i].X[j];
        residual[j]+=weights[i]*(delta-dt*s.dXdt[j])/m.M;
        rest-=weights[i]*(nuclides[j].A/mass_numbers[j]-1)*delta*constants::c*constants::c/dt;
      }
      check(next[i].X[d]>=0 && next[i].X[d]<initial_D && std::abs(next[i].sum()-1)<2e-15,
            "implicit burning keeps isotope fractions positive and normalized");
    }
    for(double v:residual)worst_residual=std::max(worst_residual,std::abs(v));
    worst_heat=std::max(worst_heat,std::abs(rest/heat-1));
  }
  check(worst_residual<1e-13,"integrated mixed-region source equations close for stiff capture",worst_residual);
  check(worst_heat<1e-7,"D-inclusive mass defect closes the integrated heat budget",worst_heat);

  const MixingRegions separate{{0,1},{1,2},{2,3},{3,4}};
  const auto finite_cold=burn_and_transport(cold,cold,network,separate,{.01,.02,.01},2,1e-14);
  double flux_residual=0,moderate_flux_residual=0,flux_backward_error=0;
  double finite_conservation=0,finite_heat_error=0;
  auto check_flux=[&](const Model& state,const std::vector<Composition>& result,
                      const std::vector<double>& diffusivity,double dt) {
    std::array<double,NSPEC> total{};double heat=0,rest=0;
    for(std::size_t i=0;i<state.size();++i) {
      const auto rate=network.eval(state.T(i),state.rho(i),result[i]);
      heat+=weights[i]*(rate.eps+rate.eps_neutrino);
      for(std::size_t k=0;k<NSPEC;++k) {
        const double delta=result[i].X[k]-state.comp[i].X[k];
        double residual=weights[i]*(delta-dt*rate.dXdt[k]);total[k]+=residual;
        double residual_scale=weights[i]*(std::abs(result[i].X[k])+std::abs(state.comp[i].X[k])
                                         +dt*std::abs(rate.dXdt[k]));
        rest-=weights[i]*(nuclides[k].A/mass_numbers[k]-1)*delta*constants::c*constants::c/dt;
        for(int j:{static_cast<int>(i)-1,static_cast<int>(i)+1})if(j>=0 && j<static_cast<int>(state.size())) {
          const double r=.5*(state.r(i)+state.r(j)),rho=.5*(state.rho(i)+state.rho(j));
          const double area_mass=4*std::acos(-1.)*r*r*rho;
          const double g=dt*area_mass*area_mass*diffusivity[std::min(i,static_cast<std::size_t>(j))]
              /std::abs(state.m[i]-state.m[j]);
          residual+=g*(result[i].X[k]-result[j].X[k]);
          residual_scale+=g*(std::abs(result[i].X[k])+std::abs(result[j].X[k]));
        }
        flux_residual=std::max(flux_residual,std::abs(residual)/state.M);
        if(dt<=1e8)moderate_flux_residual=std::max(moderate_flux_residual,std::abs(residual)/state.M);
        // Strong diffusion multiplies differences of rounded abundances by g.
        // Measure its local equation error relative to the uncancelled terms;
        // the homogeneous limit is tested separately against another solver.
        if(residual_scale>0)flux_backward_error=std::max(flux_backward_error,std::abs(residual)/residual_scale);
      }
    }
    for(double v:total)finite_conservation=std::max(finite_conservation,std::abs(v)/state.M);
    if(heat>0)finite_heat_error=std::max(finite_heat_error,std::abs(rest/heat-1));
  };
  check_flux(cold,finite_cold,{.01,.02,.01},2);
  check(finite_cold.front()[Species::H2]>cold.comp.front()[Species::H2]
        && finite_cold.back()[Species::H2]<cold.comp.back()[Species::H2],
        "finite mixing smooths the initial D gradient");
  const auto zero=burn_and_transport(m,m,network,{{0,2},{2,4}},{0,0,0},1e8,1e-14);
  check(zero==burn_and_mix(m,m,network,{{0,2},{2,4}},1e8,1e-14),
        "zero finite mixing preserves the original D burning solution exactly");
  const auto barrier=burn_and_transport(cold,cold,network,separate,{.01,0,.01},2,1e-14);
  double barrier_error=0;
  for(auto [begin,end]:MixingRegions{{0,2},{2,4}})for(std::size_t k=0;k<NSPEC;++k) {
    double sum=0;for(std::size_t i=begin;i<end;++i)sum+=weights[i]*(barrier[i].X[k]-cold.comp[i].X[k]);
    barrier_error=std::max(barrier_error,std::abs(sum)/m.M);
  }
  check(barrier_error<1e-14,"zero-flux face preserves each separate isotope inventory",barrier_error);
  const auto joined=burn_and_transport(cold,cold,network,{{0,2},{2,4}},{.01,.02,.01},2,1e-14);
  check(joined[0]==joined[1] && joined[2]==joined[3] && joined[0]!=joined[2],
        "finite D transport coexists with instantaneous convection inside each region");
  auto varied=m;varied.y[0].lnT=std::log(7e5);varied.y[2].lnT=std::log(1.5e6);
  double uniform_error=0;
  for(double dt:{1e8,1e11,1e14}) {
    const std::vector<double> finite{1e-13,2e-13,1e-13};
    const auto result=burn_and_transport(varied,varied,network,separate,finite,dt,1e-14);
    check_flux(varied,result,finite,dt);
    const auto stiff=burn_and_transport(varied,varied,network,separate,{1e30,1e30,1e30},dt,1e-14);
    const auto uniform=burn_and_mix(varied,varied,network,{{0,4}},dt,1e-14);
    for(std::size_t i=0;i<m.size();++i)for(std::size_t k=0;k<NSPEC;++k)
      uniform_error=std::max(uniform_error,std::abs(stiff[i].X[k]-uniform[i].X[k]));
  }
  check(moderate_flux_residual<2e-13,"moderate finite D mixing closes each local isotope equation",moderate_flux_residual);
  check(flux_backward_error<2e-14,"finite D continuity has small backward error including stiff diffusion",flux_backward_error);
  std::printf("Maximum unscaled local abundance residual across all stiffnesses %.4g\n",flux_residual);
  check(finite_conservation<2e-13,"finite D burning and mixing conserve the integrated reaction inventory",finite_conservation);
  check(finite_heat_error<2e-7,"finite D reaction/transport mass defect agrees with nuclear heat",finite_heat_error);
  check(uniform_error<2e-13,"arbitrarily stiff D mixing agrees with independent homogeneous burning",uniform_error);
  bool metals_rejected=false;auto bad_metals=m;
  bad_metals.comp.back().X[3]+=1e-6;bad_metals.comp.back().X[2]-=1e-6;
  try {burn_and_transport(bad_metals,bad_metals,network,separate,{.01,.02,.01},2);}
  catch(const std::invalid_argument&) {metals_rejected=true;}
  check(metals_rejected,"D-only finite transport refuses an unrepresented metal gradient");

  // Exhausting trace initial D must not freeze at the absolute Newton
  // tolerance while the rest of the composition is already well converged.
  double trace_fuel_error=0;
  for(double xd:{1e-18,1e-30}) {
    auto trace=m;auto q=c;q.X[0]+=q.X[d]-xd;q.X[d]=xd;
    for(auto& comp:trace.comp)comp=q;
    constexpr double step=1e11;
    const auto burned=burn_and_mix(trace,trace,network,{{0,trace.size()}},step,1e-16);
    double rate=0;
    for(std::size_t i=0;i<trace.size();++i) {
      auto other=burned[i];other.X[d]=0;
      rate+=weights[i]/trace.M*deuterium_capture(trace.T(i),trace.rho(i),other,burned[i].X[d]).dX_deuterium_dt;
    }
    const double change=burned.front().X[d]-xd;
    check(change<0,"trace initial D continues burning below absolute abundance tolerance",xd);
    trace_fuel_error=std::max(trace_fuel_error,std::abs(change/(step*rate)-1));
  }
  check(trace_fuel_error<2e-10,"trace initial-D consumption matches the implicit capture source",trace_fuel_error);

  const auto tag=std::chrono::steady_clock::now().time_since_epoch().count();
  auto path=std::filesystem::temp_directory_path()/("ember-d-checkpoint-"+std::to_string(tag));
  driver::Checkpoint state;state.model=m;state.next_dt=1e8;
  driver::Selections selected{"test","D","ideal","grey","none"};driver::Identities ids{{"executable","test"}};
  driver::write_checkpoint(path,state,selected,1e-12,ids);
  const auto restored=driver::read_checkpoint(path,m.size(),m.M,c,selected,1e-12,ids);
  bool exact=true;for(std::size_t i=0;i<m.size();++i)exact&=restored.model.comp[i].X==m.comp[i].X;
  check(exact,"extended checkpoint roundtrip preserves every D abundance exactly");
  std::filesystem::remove(path);
  for(auto& q:state.model.comp){q.X[0]+=q.X[d];q.X[d]=0;}
  driver::write_checkpoint(path,state,selected,1e-12,ids);
  std::ifstream legacy(path);std::string first;std::getline(legacy,first);
  const auto old=driver::read_checkpoint(path,m.size(),m.M,state.model.comp.front(),selected,1e-12,ids);
  check(first=="EMBER_EVOLUTION_CHECKPOINT 1" && old.model.comp==state.model.comp,
        "D-free checkpoints retain the original eight-species format");
  std::filesystem::remove(path);
  bool rejected=false;try{burn_and_mix(m,m,pp,{{0,4}},1e8);}catch(const std::invalid_argument&){rejected=true;}
  check(rejected,"reduced two-isotope solver cannot silently lose D");
  return failures?1:0;
}

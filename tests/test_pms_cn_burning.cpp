#include "ember/metal_cn_transport.hpp"
#include "ember/deuterium_burning.hpp"
#include "ember/constants.hpp"
#include "../src/metal_cn_source.hpp"
#include "../apps/evolution_checkpoint.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace ember;
namespace {
Composition fuel(double d,double he3=2e-5) {
  auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
  c.metal_inventory=MetalInventory::gs98;c.X[0]-=d;c[Species::H2]=d;
  c.X[1]=he3;c.X[2]-=he3;c.cn_molality=initial_gs98_cn(c);
  return explicit_cn_material(c);
}
Model model(double d,double temperature,double density=10) {
  Model m;m.M=1;m.m={.02,.08,.3,.7,1};
  for(std::size_t i=0;i<m.m.size();++i) {
    m.comp.push_back(fuel(d));
    m.y.push_back({std::log(1.+static_cast<double>(i)),std::log(density),std::log(temperature),1.});
  }
  return m;
}
}
int main() {
  int failures=0,checks=0;
  auto check=[&](bool ok,const char* label,double error=0) {
    ++checks;if(!ok)++failures;
    std::printf("%s %s %.4g\n",ok?"PASS":"FAIL",label,error);
  };
  auto rejects=[](auto f){try{f();}catch(const std::exception&){return true;}return false;};
  PPCNNetwork nuclear(PPRates::solar_fusion_iii,PPScreening::salpeter_van_horn);
  constexpr auto D=static_cast<std::size_t>(Species::H2);
  try {
    double maximum_derivative=0,maximum_energy=0,maximum_balance=0,maximum_pp_difference=0;
    for(double T:{1e6,1.4e7}) {
      const auto c=fuel(2e-5);constexpr double rho=100;
      const auto total=nuclear.composition_response(T,rho,c);
      auto other=c;other[Species::H2]=0;
      const auto capture=deuterium_capture(T,rho,other,c[Species::H2]);
      const auto pp=nuclear.pp().eval(T,rho,c);
      const auto cn=nuclear.cn().response(T,rho,c,*c.cn_molality).physical.state;
      check(std::abs(total.state.eps/(pp.eps+capture.source.eps+cn.eps)-1)<1e-13,
          "heat equals independent pp, initial-D and CN sources");
      check(total.state.eps_neutrino==pp.eps_neutrino+cn.eps_neutrino,
          "D capture adds no neutrino loss or duplicate pp energy");
      const auto source=detail::metal_cn_source(T,rho,c,nuclear.light(),nuclear.cn());
      const auto initial=metal_cn_abundances(c);
      for(std::size_t col=0;col<METAL_CN_SIZE;++col) {
        const double h=1e-5*std::max(initial[col],1e-5);
        auto up=initial,dn=initial;up[col]+=h;
        const bool edge=initial[col]<h;dn[col]+=edge?2*h:-h;
        const auto a=detail::metal_cn_source(T,rho,metal_cn_composition(c,up),nuclear.light(),nuclear.cn());
        const auto b=detail::metal_cn_source(T,rho,metal_cn_composition(c,dn),nuclear.light(),nuclear.cn());
        for(std::size_t row=0;row<METAL_CN_SIZE;++row) {
          const double difference=edge?(-3*source.source[row]+4*a.source[row]-b.source[row])/(2*h)
              :(a.source[row]-b.source[row])/(2*h);
          const double scale=std::max({std::abs(source.jacobian[row][col]),std::abs(source.source[row]),1e-30});
          maximum_derivative=std::max(maximum_derivative,std::abs(difference-source.jacobian[row][col])/scale);
        }
      }
    }
    check(maximum_derivative<2e-4,"all seven independent composition derivatives",maximum_derivative);

    for(double T:{1e6,1.4e7})for(double dt:{1e6,1e10})for(bool mixed:{false,true}) {
      auto m=model(2e-5,T);const auto weights=nodal_mass_weights(m);
      MixingRegions regions;
      if(mixed)regions={{0,m.size()}};
      else for(std::size_t i=0;i<m.size();++i)regions.push_back({i,i+1});
      const auto next=burn_and_mix(m,m,nuclear,regions,dt,1e-14);
      long double rest=0,heat=0;MetalCNVector balance{};
      double catalysts=0,inert=0;
      for(std::size_t i=0;i<m.size();++i) {
        const auto& c=next[i];const auto old=metal_cn_abundances(m.comp[i]),now=metal_cn_abundances(c);
        const auto light=nuclear.light().eval(T,10,c);
        const auto cn=nuclear.cn().response(T,10,c,*c.cn_molality).physical.state;
        const MetalCNVector source{light.dXdt[0]+cn.dXdt[0],light.dXdt[1],
          cn.dXdt[3],cn.dXdt[4],cn.dXdt[5],0,light.dXdt[D]};
        for(std::size_t k=0;k<METAL_CN_SIZE;++k)balance[k]+=weights[i]*(now[k]-old[k]-dt*source[k]);
        catalysts+=weights[i]*((now[2]-old[2])/12+(now[3]-old[3])/13+(now[4]-old[4])/14);
        inert+=weights[i]*(now[5]-old[5]);
        const auto power=nuclear.eval(T,10,c);heat+=weights[i]*(power.eps+power.eps_neutrino);
        for(std::size_t k=0;k<NSPEC;++k)
          rest-=static_cast<long double>(weights[i])*(nuclides[k].A/mass_numbers[k]-1)
            *(static_cast<long double>(c.X[k])-m.comp[i].X[k])*constants::c*constants::c/dt;
        rest-=static_cast<long double>(weights[i])*(static_cast<long double>(nuclear.rest_energy_correction(c))
            -nuclear.rest_energy_correction(m.comp[i]))/dt;
        check(c.X[D]>=0 && c.X[D]<m.comp[i].X[D] && std::abs(c.sum()-1)<2e-15,
            "combined solve consumes D with positive normalized abundances");
      }
      for(double b:balance)maximum_balance=std::max(maximum_balance,std::abs(b));
      maximum_energy=std::max(maximum_energy,std::abs(static_cast<double>(rest/heat-1)));
      check(std::abs(catalysts)<1e-14 && std::abs(inert)<1e-14,
          "CN catalyst number and inert metal mass are conserved separately");
      if(T==1e6) {
        auto d_only=m;
        for(auto& c:d_only.comp){c.cn_molality.reset();c.cn_mass_convention=CNMassConvention::fixed_metal_proxy;}
        PPDeuterium pp;
        const auto reference=burn_and_mix(d_only,d_only,pp,regions,dt,1e-14);
        for(std::size_t i=0;i<m.size();++i)for(auto k:{0UL,1UL,2UL,D})
          maximum_pp_difference=std::max(maximum_pp_difference,std::abs(next[i].X[k]-reference[i].X[k]));
      }
    }
    check(maximum_balance<5e-14,"independently integrated isotope equations",maximum_balance);
    check(maximum_energy<2e-5,"combined mass defect equals nuclear heat plus neutrinos",maximum_energy);
    check(maximum_pp_difference<5e-14,"cold combined network agrees with independent PMS solver",maximum_pp_difference);

    // Cold finite mixing moves every physical isotope, including D and CN,
    // while keeping the global inventories fixed. The closed barrier is exact.
    auto cold=model(2e-5,3000,1);const auto weights=nodal_mass_weights(cold);
    MixingRegions separate;
    for(std::size_t i=0;i<cold.size();++i) {
      auto v=metal_cn_abundances(cold.comp[i]);v[METAL_CN_D]*=1.+static_cast<double>(i);v[2]*=1.+.1*static_cast<double>(i);
      cold.comp[i]=metal_cn_composition(cold.comp[i],v);separate.push_back({i,i+1});
    }
    const auto mixed=burn_and_transport(cold,cold,nuclear,separate,{.01,.01,0.,.01},2.,1e-14);
    double conservation=0;
    for(std::size_t k=0;k<METAL_CN_SIZE;++k) {
      double a=0,b=0;
      for(std::size_t i=0;i<cold.size();++i) {
        const double change=weights[i]*(metal_cn_abundances(mixed[i])[k]-metal_cn_abundances(cold.comp[i])[k]);
        (i<3?a:b)+=change;
      }
      conservation=std::max({conservation,std::abs(a),std::abs(b)});
    }
    check(conservation<1e-14 && mixed.front()[Species::H2]>cold.comp.front()[Species::H2],
        "finite mixing conserves D and metals on each side of a closed barrier",conservation);
    check(rejects([&]{common_metal_cn_flux({},cold.comp.front(),cold.comp.back(),true);}),
        "three-mass microscopic provider cannot silently drop D transport");

    // Compare the trace-D update with its scalar implicit destruction law,
    // including values whose unscaled nuclear rate underflows.
    double trace_error=0;
    for(double d:{1e-15,1e-100,1e-250,1e-310,1e-320,std::numeric_limits<double>::denorm_min()}) {
      auto m=model(d,1e6);constexpr double dt=1e12;
      const auto next=burn_and_mix(m,m,nuclear,{{0,m.size()}},dt,1e-14);
      const auto response=nuclear.light().composition_response(1e6,10,next.front());
      const double expected=d/(1-dt*response.d_dXdt_dX[D][D]);
      const double error=std::abs(next.front()[Species::H2]-expected);
      trace_error=std::max(trace_error,error/std::max(d,1e-300));
      check(error<=1e-9*d+2*std::numeric_limits<double>::denorm_min(),
          "trace D follows the same implicit depletion through underflow",error);
    }
    check(trace_error<1e-9,"maximum scaled trace-D error",trace_error);

    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    const auto path=std::filesystem::temp_directory_path()/("ember-pms-cn-"+std::to_string(stamp)+".checkpoint");
    driver::Selections selections{"D+pp+CN","SFIII","control","control","mixing"};
    driver::Identities identities{{"executable","combined inventories test"}};
    for(double d:{2e-5,1e-310,1e-320,std::numeric_limits<double>::denorm_min(),0.}) {
      const auto c=fuel(d);const auto m=model(d,1e6);
      driver::write_checkpoint(path,{m,1e6,0,0},selections,1e-14,identities);
      const auto saved=driver::read_checkpoint(path,m.size(),m.M,c,selections,1e-14,identities);
      check(saved.model.comp==m.comp,"checkpoint preserves D and actual CN together exactly");
      std::filesystem::remove(path);
    }
  } catch(const std::exception& e){++failures;std::printf("FAIL exception: %s\n",e.what());}
  std::printf("checks %d failures %d\n",checks,failures);
  return failures?1:0;
}

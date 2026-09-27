#include "ember/cn_transport.hpp"
#include "ember/constants.hpp"
#include "../src/cn_source.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ember;
namespace {
using Six=std::array<double,6>;
Six physical(const Composition& c) {
  const auto p=cn_physical_ledger(c,*c.cn_molality);const auto x=cn_transport_abundances(c);
  return {x[0],x[1],p.helium4,x[2],x[3],x[4]};
}
Composition composition(double h,double conversion,double enrichment) {
  auto c=solar_scaled(h,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=.001;c.X[2]-=.001;const auto y=initial_gs98_cn(c);
  c.cn_molality=CNAbundances{enrichment*y[0]*(1-conversion),enrichment*.03*y[0],
    enrichment*(y[2]+(conversion-.03)*y[0])};return c;
}
SpeciesFaceResponse hhe_flux(const Composition& a,const Composition& b,double factor) {
  SpeciesFaceResponse f;f.rate={factor*(.1+a.X[0]-b.X[0]),factor*.3*(a.X[1]-b.X[1])};
  f.dleft[0][0]=factor;f.dright[0][0]=-factor;f.dleft[1][1]=.3*factor;f.dright[1][1]=-.3*factor;return f;
}
}
int main() {
  int checks=0,failures=0;double max_derivative=0,max_balance=0,max_energy=0,max_region=0;
  auto check=[&](bool ok,const char* label,double error=0.){++checks;if(!ok){++failures;std::printf("FAIL %s %.8g\n",label,error);}};
  auto rejects=[](auto f){try{f();}catch(const std::exception&){return true;}return false;};
  PPCNNetwork nuclear;auto left=composition(.2,.1,.7),right=composition(.55,.8,1.3);
  for(double sign:{-1.,1.}) {
    const auto baseline=hhe_flux(left,right,sign);
    const auto f=trace_cn_flux(baseline,left,right,CNMicroscopicApproximation::helium_velocity,true);
    const double group_flux=-(baseline.rate[0]+baseline.rate[1]);const auto& donor=group_flux>0?left:right;
    const auto x=physical(donor);const double group=x[2]+x[3]+x[4]+x[5];double sum=0;for(double rate:f.rate)sum+=rate;
    check(std::abs(-sum-group_flux*x[2]/group)<1e-15,"helium plus CN preserves exact baryonic closure");
    for(std::size_t side=0;side<2;++side)for(std::size_t col=0;col<5;++col) {
      auto a=cn_transport_abundances(side?right:left);const double delta=1e-6*std::max(a[col],1e-5);
      auto up=a,dn=a;up[col]+=delta;dn[col]-=delta;
      auto au=left,ad=left,bu=right,bd=right;
      if(side){bu=cn_transport_composition(right,up);bd=cn_transport_composition(right,dn);}
      else{au=cn_transport_composition(left,up);ad=cn_transport_composition(left,dn);}
      const auto plus=trace_cn_flux(hhe_flux(au,bu,sign),au,bu,CNMicroscopicApproximation::helium_velocity,false);
      const auto minus=trace_cn_flux(hhe_flux(ad,bd,sign),ad,bd,CNMicroscopicApproximation::helium_velocity,false);
      for(std::size_t row=0;row<5;++row) {
        const double expected=side?f.dright[row][col]:f.dleft[row][col];
        const double error=std::abs((plus.rate[row]-minus.rate[row])/(2*delta)-expected)/std::max(1e-8,std::abs(expected));
        max_derivative=std::max(max_derivative,error);check(error<2e-6,"analytic five-variable trace flux derivative",error);
      }
    }
  }
  const auto src=detail::cn_full_source(1.4e7,300,left,nuclear.pp(),nuclear.cn());
  for(std::size_t col=0;col<5;++col) {
    const auto x=cn_transport_abundances(left);auto up=x,dn=x;const double h=1e-5*std::max(x[col],1e-5);up[col]+=h;dn[col]-=h;
    const auto a=detail::cn_full_source(1.4e7,300,cn_transport_composition(left,up),nuclear.pp(),nuclear.cn());
    const auto b=detail::cn_full_source(1.4e7,300,cn_transport_composition(left,dn),nuclear.pp(),nuclear.cn());
    for(std::size_t row=0;row<5;++row) {
      const double derivative=(a.source[row]-b.source[row])/(2*h),expected=src.jacobian[row][col];
      check(std::abs(derivative-expected)<2e-6*std::max(1e-25,std::abs(expected)),"full independent isotope source derivatives",derivative-expected);
    }
  }
  Model m;m.M=10;m.m={1,4,10};m.comp={left,composition(.35,.5,1.),right};m.y.assign(3,{0,0,0,1});
  const auto w=nodal_mass_weights(m);
  for(bool warm:{false,true})for(int mode=0;mode<5;++mode) {
    for(std::size_t i=0;i<m.size();++i){m.y[i].lnT=std::log(warm?1.2e7+1e6*static_cast<double>(i):1e4);m.y[i].lnrho=std::log(warm?300.:1.);}
    const double dt=warm?1e13:2.;const MixingRegions regions=mode==3?MixingRegions{{0,2},{2,3}}:MixingRegions{{0,1},{1,2},{2,3}};
    CNSpeciesFlux flux=[&](std::size_t,const Composition& a,const Composition& b,bool derivatives) {
      CNSpeciesFaceResponse f;
      if(mode>=2)f=trace_cn_flux(hhe_flux(a,b,1/dt),a,b,CNMicroscopicApproximation::helium_velocity,derivatives);
      return f;
    };
    SpeciesTransportOptions options;options.abundance_tolerance=1e-14;
    const std::vector<double> common(2,mode>0 && mode<4?.3/dt:0.);
    const auto result=burn_cn_and_diffuse(m,m,nuclear,regions,flux,dt,options,common);
    if(warm && mode==0) {
      auto separate=options;separate.abundance_tolerance=1e-6;
      const auto loose=burn_cn_and_diffuse(m,m,nuclear,regions,flux,dt,separate,common);
      double loose_balance=0;for(double x:loose.integrated_balance)loose_balance=std::max(loose_balance,std::abs(x));
      check(loose_balance>1e-14,"independent-balance control exercises a loose residual",loose_balance);
      separate.integrated_balance_tolerance=1e-14;
      const auto bounded=burn_cn_and_diffuse(m,m,nuclear,regions,flux,dt,separate,common);
      for(double x:bounded.integrated_balance)check(std::abs(x)<=1e-14,
          "independent integrated CN balance remains enforced",x);
      check(bounded.abundance_correction<=separate.abundance_tolerance,"independent local CN correction bound");
    }
    Six global{};double heat=0,rest=0;std::vector<Six> residual;
    for(std::size_t i=0;i<m.size();++i) {
      const auto a=physical(m.comp[i]),b=physical(result.composition[i]);
      const auto pp=nuclear.pp().eval(m.T(i),m.rho(i),result.composition[i]);
      const auto cn=nuclear.cn().response(m.T(i),m.rho(i),result.composition[i],*result.composition[i].cn_molality).physical.state;
      Six e{};
      for(std::size_t k=0;k<6;++k) {
        check(b[k]>=0,"positive physical species");e[k]=w[i]/m.M*(b[k]-a[k]-dt*(pp.dXdt[k]+cn.dXdt[k]));global[k]+=e[k];
        rest-=w[i]/dt*(b[k]-a[k])*(nuclides[k].A/mass_numbers[k]-1)*constants::c*constants::c;
      }
      residual.push_back(e);heat+=w[i]*(pp.eps+pp.eps_neutrino+cn.eps+cn.eps_neutrino);
    }
    for(double e:global){max_balance=std::max(max_balance,std::abs(e));check(std::abs(e)<2e-14,"global physical fuel and catalyst conservation",e);}
    if(warm){const double e=std::abs(rest/heat-1);max_energy=std::max(max_energy,e);check(e<2e-6,"physical nuclear mass-energy balance",e);}
    for(std::size_t region=0;region<regions.size();++region) {
      const auto [begin,end]=regions[region];Six e{};
      for(std::size_t i=begin;i<end;++i)for(std::size_t k=0;k<6;++k)e[k]+=residual[i][k];
      for(int side:{-1,1}) {
        if((side<0 && region==0) || (side>0 && region+1==regions.size()))continue;
        const auto& f=result.boundary_fluxes[side<0?region-1:region].rate;
        Six all{f[0],f[1],-(f[0]+f[1]+f[2]+f[3]+f[4]),f[2],f[3],f[4]};
        for(std::size_t k=0;k<6;++k)e[k]+=side*dt/m.M*all[k];
      }
      for(double v:e){max_region=std::max(max_region,std::abs(v));check(std::abs(v)<2e-14,"independent regional physical flux equation",v);}
      for(std::size_t i=begin+1;i<end;++i)check(result.composition[i]==result.composition[begin],"mixed region homogeneous including total CN");
    }
    if(mode)check(result.composition.front().cn_molality!=m.comp.front().cn_molality,"local catalyst inventory actually transported");
  }
  check(rejects([&]{trace_cn_flux({},left,right,CNMicroscopicApproximation::unselected,true);}),"missing approximation rejected");
  const auto stopped=trace_cn_flux(hhe_flux(left,right,1),left,right,CNMicroscopicApproximation::zero_catalyst_drift,true);
  check(stopped.rate[2]==0 && stopped.rate[3]==0 && stopped.rate[4]==0,"explicit zero-catalyst-drift control");
  // Independent evaluations must preserve complete solutions and the first
  // domain error. A larger partition exercises all requested worker counts.
  Model parallel;parallel.M=160;
  MixingRegions partition;
  for(std::size_t i=0;i<160;++i) {
    parallel.m.push_back(static_cast<double>(i+1));
    parallel.y.push_back({0,0,std::log(1e4),1});
    parallel.comp.push_back(composition(.2+.3*static_cast<double>(i)/160,.1+.4*static_cast<double>(i)/160,1.));
    partition.push_back({i,i+1});
  }
  const CNSpeciesFlux parallel_flux=[](std::size_t,const Composition& a,const Composition& b,bool derivatives) {
    return trace_cn_flux(hhe_flux(a,b,.1),a,b,CNMicroscopicApproximation::helium_velocity,derivatives);
  };
  SpeciesTransportOptions parallel_options;parallel_options.abundance_tolerance=1e-13;
  const auto serial=burn_cn_and_diffuse(parallel,parallel,nuclear,partition,parallel_flux,2.,parallel_options);
  for(std::size_t threads:{2,4}) {
    parallel_options.evaluation_threads=threads;
    const auto result=burn_cn_and_diffuse(parallel,parallel,nuclear,partition,parallel_flux,2.,parallel_options);
    check(result.composition==serial.composition,"parallel CN composition identical");
    check(result.integrated_balance==serial.integrated_balance && result.residual_history==serial.residual_history
      && result.abundance_correction==serial.abundance_correction && result.iterations==serial.iterations,
      "parallel CN conservation and Newton history identical");
    for(std::size_t i=0;i<serial.boundary_fluxes.size();++i)
      check(result.boundary_fluxes[i].face==serial.boundary_fluxes[i].face
        && result.boundary_fluxes[i].rate==serial.boundary_fluxes[i].rate,"parallel CN boundary flux identical");
  }
  const CNSpeciesFlux failing_flux=[](std::size_t face,const Composition&,const Composition&,bool) {
    if(face==3 || face==90)throw std::domain_error("deliberate face "+std::to_string(face));
    return CNSpeciesFaceResponse{};
  };
  for(std::size_t threads:{1,2,4}) {
    parallel_options.evaluation_threads=threads;
    std::string message;
    try{(void)burn_cn_and_diffuse(parallel,parallel,nuclear,partition,failing_flux,2.,parallel_options);}
    catch(const std::domain_error& e){message=e.what();}
    check(message=="deliberate face 3","parallel CN reports first failing face");
  }
  parallel_options.evaluation_threads=0;
  check(rejects([&]{burn_cn_and_diffuse(parallel,parallel,nuclear,partition,parallel_flux,2.,parallel_options);}),
    "zero CN evaluation threads rejected");
  parallel_options.evaluation_threads=65;
  check(rejects([&]{burn_cn_and_diffuse(parallel,parallel,nuclear,partition,parallel_flux,2.,parallel_options);}),
    "excess CN evaluation threads rejected");
  std::printf("%d checks, %d failures; max derivative %.4g region %.4g global %.4g mass-energy %.4g\n",
    checks,failures,max_derivative,max_region,max_balance,max_energy);
  return failures?1:0;
}

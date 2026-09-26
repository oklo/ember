#include "ember/eos_deuterium.hpp"
#include "ember/eos_composite.hpp"
#include "ember/eos_variable_metal.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ember;
int main(int argc,char** argv) {
  try {
    CompositeEos direct;DeuteriumApproxEos mapped(direct);
    double discrepancy=0,derivative_error=0,inversion_error=0;int failures=0;
    const auto check=[&](bool pass,const char* label,double value=0.) {
      failures+=!pass;std::printf("[%s] %s %.12g\n",pass?"PASS":"FAIL",label,value);
    };
    for(bool cn:{false,true})for(double D:{0.,2e-5,1e-4})for(auto [T,rho]:{std::pair{3000.,1e-8},
        std::pair{1e6,10.},std::pair{4e6,100.},std::pair{1e7,1e5}}) {
      auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
      c.metal_inventory=MetalInventory::gs98;c.X[0]-=D;c[Species::H2]=D;
      if(cn){c.cn_molality=initial_gs98_cn(c);c=explicit_cn_material(c);}
      const auto exact=direct.eval_with_derivatives(T,rho,c);
      const auto out=mapped.eval_with_derivatives(T,rho,c);
      for(auto [a,b]:{std::pair{out.state.P,exact.state.P},std::pair{out.state.E,exact.state.E},
          std::pair{out.state.S,exact.state.S},std::pair{out.state.cp,exact.state.cp},
          std::pair{out.state.cv,exact.state.cv},std::pair{out.state.grad_ad,exact.state.grad_ad},
          std::pair{out.state.mu,exact.state.mu},std::pair{out.state.free_e,exact.state.free_e}})
        discrepancy=std::max(discrepancy,std::abs(a/b-1));
      const double h=1e-5;
      for(int j:{0,1}) {
        const auto p=mapped.eval(T*std::exp(j==0?h:0),rho*std::exp(j==1?h:0),c);
        const auto m=mapped.eval(T*std::exp(j==0?-h:0),rho*std::exp(j==1?-h:0),c);
        const double de=(p.E-m.E)/(2*h),analytic=j==0?out.state.cv*T:out.dE_dlnRho;
        derivative_error=std::max(derivative_error,std::abs(de-analytic)/out.state.E);
        const double dcp=(p.cp-m.cp)/(2*h),acp=j==0?out.dcp_dlnT:out.dcp_dlnRho;
        derivative_error=std::max(derivative_error,std::abs(dcp-acp)/out.state.cp);
      }
      inversion_error=std::max(inversion_error,std::abs(mapped.rho_from_PT(T,out.state.P,c,rho*.9)/rho-1));
    }
    check(discrepancy<2e-12,"number mapping and isotope entropy recover direct ion/electron/radiation EOS",discrepancy);
    check(derivative_error<2e-8,"mapped caloric derivatives agree with finite differences",derivative_error);
    check(inversion_error<2e-8,"mapped PT inversion recovers density",inversion_error);
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
    c.X[0]-=.001;c[Species::H2]=.001;bool rejected=false;
    try{mapped.eval(1e6,10,c);}catch(const std::domain_error&){rejected=true;}
    check(rejected,"non-trace D rejected without capping or dropping fuel");
    struct RecordingEos final:Eos {
      CompositeEos source;mutable Composition last;mutable double density{};
      EosState eval(double T,double rho,const Composition& q)const override {
        last=q;density=rho;(void)cn_physical_ledger(q,*q.cn_molality);return source.eval(T,rho,q);
      }
      const char* name()const override{return "record material composition mapping";}
    } recording;
    c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[0]-=2e-5;c[Species::H2]=2e-5;c.cn_molality=initial_gs98_cn(c);c=explicit_cn_material(c);
    const auto unchanged=c;
    DeuteriumApproxEos with_cn(recording);(void)with_cn.eval(1e6,10,c);
    const auto physical=cn_physical_ledger(c,*c.cn_molality);
    const auto proxy=cn_physical_ledger(recording.last,*recording.last.cn_molality);
    const double f=1-c[Species::H2]/2;double number_error=0;
    for(std::size_t i=0;i<3;++i)number_error=std::max(number_error,
        std::abs(recording.density*(*recording.last.cn_molality)[i]-10*(*c.cn_molality)[i]));
    check(number_error<1e-17 && std::abs(f*proxy.metal_fraction-physical.metal_fraction)<1e-17,
          "D EOS mapping preserves CN number density and physical metal mass",number_error);
    check(c==unchanged,"material mapping does not alter the nuclear fuel inventory");
    struct RecordingAtmosphere final:Atmosphere {
      mutable Composition last;
      AtmosphereState eval(double T,double,const Composition& q)const override {
        last=q;(void)cn_physical_ledger(q,*q.cn_molality);
        AtmosphereState a;a.T=T;a.P=1e5;a.rho=1e-7;a.dlnT_dlnTeff=1;return a;
      }
      const char* name()const override{return "record atmosphere composition mapping";}
    } atmosphere_source;
    TraceDeuteriumAtmosphere atmosphere(mapped,atmosphere_source);
    const auto boundary=atmosphere.eval(3000,1e5,c);
    check(atmosphere_source.last.cn_molality==c.cn_molality
        && atmosphere_source.last.cn_mass_convention==c.cn_mass_convention
        && atmosphere_source.last[Species::H2]==0
        && atmosphere_source.last.X[0]==c.X[0]+c[Species::H2],
        "trace-D atmosphere keeps CN and metal mass while looking up total H mass");
    check(std::abs(mapped.eval(boundary.T,boundary.rho,c).P/boundary.P-1)<2e-8,
        "CN/D atmosphere density uses the actual interior isotope mixture");
    auto invalid=c;(*invalid.cn_molality)[0]=1;rejected=false;
    try{with_cn.eval(1e6,10,invalid);}catch(const std::domain_error&){rejected=true;}
    check(rejected,"D mapping rejects catalysts exceeding available metal mass");
    if(argc==2) {
      VariableMetalHelmholtzEos source(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
      DeuteriumApproxEos eos(source);double response_error=0,maxwell_error=0,trace_effect=0;
      auto mixture=solar_scaled(.7,.02);mixture.basis=AbundanceBasis::baryon_mass;
      mixture.metal_inventory=MetalInventory::gs98;mixture[Species::H2]=2e-5;mixture.X[0]-=2e-5;
      mixture.X[1]=.003;mixture.X[2]-=.003;
      for(auto [T,rho]:{std::pair{4000.,1e-7},std::pair{32000.,1e-4},
          std::pair{5e5,10.},std::pair{3e6,100.}}) {
        const auto r=eos.eval_with_derivatives(T,rho,mixture);const auto& s=r.state;
        const auto cr=eos.composition_response(T,rho,mixture);const double h=1e-5;
        for(int j:{0,1}) {
          auto plus=mixture,minus=mixture;plus.X[j]+=h;plus.X[2]-=h;minus.X[j]-=h;minus.X[2]+=h;
          const auto p=eos.eval(T,rho,plus),m=eos.eval(T,rho,minus);
          response_error=std::max({response_error,std::abs((p.P-m.P)/(2*h)-cr.dP[j])/s.P,
            std::abs((p.E-m.E)/(2*h)-cr.dE[j])/(T*s.cv)});
        }
        const auto a=eos.eval(T,rho*std::exp(h),mixture),b=eos.eval(T,rho*std::exp(-h),mixture);
        maxwell_error=std::max({maxwell_error,std::abs((a.E-b.E)/(2*h)-s.P/rho*(1-s.chiT))/(T*s.cv),
          std::abs((a.S-b.S)/(2*h)+s.P*s.chiT/(rho*T))/s.cv});
        auto nominal=mixture;nominal.X[0]+=nominal[Species::H2];nominal[Species::H2]=0;
        const auto ordinary=source.eval(T,rho,nominal);
        trace_effect=std::max({trace_effect,std::abs(s.P/ordinary.P-1),std::abs(s.cp/ordinary.cp-1),
          std::abs(s.grad_ad/ordinary.grad_ad-1)});
        check(std::abs(eos.rho_from_PT(T,s.P,mixture,rho)/rho-1)<1e-9,"material-table PT inversion");
      }
      check(response_error<1e-5,"table H/He3 composition derivatives at fixed D",response_error);
      check(maxwell_error<1e-5,"number-mapped table preserves first-law and Maxwell identities",maxwell_error);
      check(trace_effect<1e-3,"trace material effect on P/cp/adiabatic gradient in sampled states",trace_effect);
    }
    return failures?1:0;
  }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}
}

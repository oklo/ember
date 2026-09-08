#include "ember/eos_helmholtz.hpp"
#include "ember/opacity_aesopus.hpp"
#include "ember/opacity_tops.hpp"
#include "ember/boundary.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>

using namespace ember;
static int failures=0;
static void check(bool ok,const char* label,double value=0) {
  if (!ok) ++failures;
  std::printf("[%s] %s (%.9g)\n",ok?"PASS":"FAIL",label,value);
}
template<class F> bool throws(F f) { try {f();} catch(const std::exception&) {return true;} return false; }

int main() {
  const std::string data=EMBER_DATA_DIR;
  const auto file=data+"/eos/freeeos300_hhe_x070_potential.dat";
  HelmholtzTableEos eos(file,HelmholtzTableEos::Mixture::allow_documented_proxy);
  const auto comp=solar_scaled(.7,.02);
  check(throws([&]{HelmholtzTableEos strict(file);}),"metals-as-helium proxy requires explicit selection");
  auto changed=comp; changed[Species::He3]=.001; changed[Species::He4]-=.001;
  check(throws([&]{eos.eval(1e6,1,changed);}),"fixed table rejects He3/composition changes");
  check(throws([&]{eos.eval(1000,1e-5,comp);}) && throws([&]{eos.eval(3000,1e3,comp);})
     && throws([&]{eos.rho_from_PT(3000,1e30,comp);}),"temperature/density/pressure extrapolation is rejected");
  check(eos.has_internal_energy(),"caloric energy is provided by the same potential");
  const std::array<std::pair<double,double>,10> states{{
    {3.371,-6.773},{3.731,-4.881},{4.173,-5.19},{4.461,-3.17},
    {5.173,-1.177},{6.119,.17},{6.713,2.31},{4.42549,-.931112},
    {3.9572579388,-2.1821411875},{5.2066765669,.6246495827}}};
  constexpr double h=1e-5;
  double identity=0,transport=0,inversion=0;
  for (auto [lt,lr]:states) {
    const double T=std::pow(10.,lt),rho=std::pow(10.,lr);
    const auto a=eos.eval_with_derivatives(T,rho,comp);const auto& e=a.state;
    inversion=std::max(inversion,std::abs(eos.rho_from_PT(T,e.P,comp)/rho-1));
    for (int v=0;v<2;++v) {
      const auto p=eos.eval(T*std::exp(v==0?h:0),rho*std::exp(v==1?h:0),comp);
      const auto m=eos.eval(T*std::exp(v==0?-h:0),rho*std::exp(v==1?-h:0),comp);
      const double dE=(p.E-m.E)/(2*h),dS=(p.S-m.S)/(2*h);
      const double targetE=v==0?T*e.cv:e.P/rho*(1-e.chiT);
      const double targetS=v==0?e.cv:-e.P*e.chiT/(rho*T);
      identity=std::max({identity,std::abs(dE-targetE)/std::max(T*e.cv,std::abs(targetE)),
                         std::abs(dS-targetS)/std::max(e.cv,std::abs(targetS))});
      auto diff=[&](double plus,double minus,double analytic,double scale) {
        const double numerical=(plus-minus)/(2*h);
        transport=std::max(transport,std::abs(numerical-analytic)/std::max({scale,std::abs(numerical),std::abs(analytic)}));
      };
      diff(std::log(p.P),std::log(m.P),v==0?e.chiT:e.chiRho,1);
      diff(p.cp,m.cp,v==0?a.dcp_dlnT:a.dcp_dlnRho,e.cp);
      diff(p.delta,m.delta,v==0?a.ddelta_dlnT:a.ddelta_dlnRho,1);
      diff(p.grad_ad,m.grad_ad,v==0?a.dgrad_ad_dlnT:a.dgrad_ad_dlnRho,1);
    }
  }
  check(identity<3e-6,"independent finite differences of E and S satisfy both first-law identities",identity);
  check(transport<3e-5,"pressure and all transport responses differentiate the actual potential",transport);
  check(inversion<1e-10,"safeguarded PT inversion recovers density",inversion);
  // Source-data regressions use off-grid, independently queried FreeEOS
  // states, rather than values reconstructed by the table importer.
  std::ifstream in(std::string(EMBER_TEST_DATA_DIR)+"/freeeos300_hhe_x070_reference.dat");
  double T,rho,P,E,S,cv,cp,ad;int count=0;double value_error=0,response_error=0;
  while(in>>T>>rho>>P>>E>>S>>cv>>cp>>ad) {
    const auto e=eos.eval(T,rho,comp);++count;
    value_error=std::max({value_error,std::abs(e.P/P-1),std::abs(e.E/E-1),std::abs(e.S/S-1)});
    response_error=std::max({response_error,std::abs(e.cv/cv-1),std::abs(e.cp/cp-1),std::abs(e.grad_ad/ad-1)});
  }
  check(count==10 && value_error<.001,"independent FreeEOS thermodynamic values agree within 0.1%",value_error);
  check(response_error<.02,"source-fit join effects on heat capacities/adiabat stay below 2% at audit points",response_error);
  // C2 continuity: P, E, S and heat capacities stay continuous across a
  // table knot even though their next derivatives can jump.
  double continuity=0;
  for (double lt:{3.7,4.425,5.2,6.5}) {
    const double r=std::pow(10.,.913+1.5*(lt-6)),t=std::pow(10.,lt);
    const auto p=eos.eval(t*std::exp(1e-10),r,comp),m=eos.eval(t*std::exp(-1e-10),r,comp);
    continuity=std::max({continuity,std::abs(p.P/m.P-1),std::abs(p.E/m.E-1),std::abs(p.S/m.S-1),std::abs(p.cp/m.cp-1)});
  }
  check(continuity<1e-6,"thermodynamic values and heat capacities are continuous across knots",continuity);
  AesopusOpacity low(data+"/opacity/aesopus21_gs98_z020.dat");TopsOpacity high(data+"/opacity");
  BlendedOpacity opacity(low,high,4.4,4.5);PPChains nuclear;
  Physics physics{&eos,&opacity,&nuclear,1.9};
  double jacobian=0;
  for (auto [lt,lr]:{std::pair{3.731,-4.881},std::pair{4.42549,-.931112},std::pair{6.319,1.079}}) {
    Model m;m.M=.1*constants::Msun;m.m={.4*m.M,.41*m.M};m.comp={comp,comp};
    m.y={{std::log(.08*constants::Rsun),lr*std::log(10.),lt*std::log(10.),.001*constants::Lsun},
         {std::log(.083*constants::Rsun),lr*std::log(10.)-.1,lt*std::log(10.)-.04,.0011*constants::Lsun}};
    Model prev=m;for(auto& p:prev.y) {p.lnT-=.004;p.lnrho-=.006;}
    for (double dt:{0.,1e12}) {
      const auto a=zone_residual(m,0,physics,dt,&prev),b=zone_residual_numerical(m,0,physics,dt,&prev,1e-6);
      const double dm=m.m[1]-m.m[0];
      for(int k=0;k<4;++k) for(int v=0;v<4;++v) {
        const double scale=dm*(v==3?m.y[0].L:1)/(k==2?m.y[0].L:1);
        for(auto pair:{std::pair{a.dfdy_lo[k][v],b.dfdy_lo[k][v]},std::pair{a.dfdy_hi[k][v],b.dfdy_hi[k][v]}})
          jacobian=std::max(jacobian,std::abs((pair.first-pair.second)*scale)/std::max(1.,std::abs(pair.second*scale)));
      }
    }
    const auto center=central_residual(m,physics,1e12,&prev);
    check(std::isfinite(center.f[1]),"fixed-composition central energy boundary accepts validated caloric state");
  }
  check(jacobian<5e-5,"static AND time-dependent energy/MLT zone Jacobians match independent differences",jacobian);
  return failures?1:0;
}

#include "ember/evolution.hpp"
#include "ember/convection.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* label,double error=0) {
  failures+=!ok;std::printf("[%s] %s (%.8g)\n",ok?"PASS":"FAIL",label,error);
}
class Reaction final:public Nuclear {
public:
  double rate{};
  NuclearState eval(double,double,const Composition& c) const override {
    NuclearState s{};s.dXdt[0]=-rate*c.X[0]*c.X[0];s.dXdt[2]=-s.dXdt[0];return s;
  }
  NuclearResponse composition_response(double T,double rho,const Composition& c) const override {
    NuclearResponse r{};r.state=eval(T,rho,c);r.d_dXdt_dX[0][0]=-2*rate*c.X[0];r.d_dXdt_dX[2][0]=2*rate*c.X[0];return r;
  }
  const char* name() const override {return "test quadratic conversion";}
};
class Gas final:public Eos {
public:
  EosState eval(double T,double rho,const Composition& c) const override {
    const double R=constants::R_gas*(c.mu_ions_inv()+c.mu_elec_inv());
    EosState s{};s.P=R*rho*T;s.cv=1.5*R;s.cp=2.5*R;s.delta=s.chiT=s.chiRho=1;s.grad_ad=.4;return s;
  }
  const char* name() const override {return "test ideal mixture";}
};
class ConstantOpacity final:public Opacity {
public:
  OpacityState eval(double,double,const Composition&) const override {return {1,0,0};}
  const char* name() const override {return "test constant";}
};
Composition comp(double x,double y=.01) {
  auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.X[1]=y;c.X[2]-=y;return c;
}
}
int main() {
  Reaction nuclear;Model m;m.M=10;m.m={1,4,10};m.comp={comp(.3),comp(.5),comp(.7)};
  m.y.assign(3,{0,0,std::log(1e6),0});const auto w=nodal_mass_weights(m);
  MixingRegions regions{{0,1},{1,2},{2,3}};std::vector<double> D{.01,.01};
  const auto mixed=burn_and_transport(m,m,nuclear,regions,D,2);
  double before=0,after=0;for(int i=0;i<3;++i){before+=w[i]*m.comp[i].X[0];after+=w[i]*mixed[i].X[0];}
  check(std::abs(after-before)<1e-13,"finite-volume diffusion conserves total hydrogen",after-before);
  check(mixed[0].X[0]>.3 && mixed[2].X[0]<.7,"diffusion smooths a gradient without overshoot");
  double residual=0;
  for(int i=0;i<3;++i) {
    double f=w[i]*(mixed[i].X[0]-m.comp[i].X[0]);
    for(int j:{i-1,i+1})if(j>=0 && j<3) {
      const double g=2*16*M_PI*M_PI*D[std::min(i,j)]/std::abs(m.m[i]-m.m[j]);
      f+=g*(mixed[i].X[0]-mixed[j].X[0]);
    }
    residual=std::max(residual,std::abs(f));
  }
  check(residual<1e-13,"solution satisfies independently assembled backward-Euler flux equations",residual);
  const auto stopped=burn_and_transport(m,m,nuclear,regions,{0,0},2);
  const auto original=burn_and_mix(m,m,nuclear,regions,2);
  check(stopped[0].X==original[0].X && stopped[2].X==original[2].X,
        "zero coefficients reproduce existing local burning exactly");
  const auto stiff=burn_and_transport(m,m,nuclear,regions,{1e30,1e30},2);
  check(std::abs(stiff[0].X[0]-before/m.M)<2e-15 && stiff[0].X==stiff[2].X,
        "extremely stiff diffusion recovers mass-weighted uniform limit",stiff[0].X[0]-before/m.M);
  const auto barrier=burn_and_transport(m,m,nuclear,regions,{1e30,0},2);
  check(std::abs(barrier[0].X[0]-(w[0]*.3+w[1]*.5)/(w[0]+w[1]))<1e-14 && std::abs(barrier[2].X[0]-.7)<1e-15,
        "zero-flux face preserves a composition barrier");
  const auto collapsed=burn_and_transport(m,m,nuclear,{{0,2},{2,3}},{0,.01},2);
  check(collapsed[0].X==collapsed[1].X,"instantaneous convection and finite diffusion coexist");
  auto inconsistent=m;inconsistent.comp.back().metal_inventory=MetalInventory::gs98;
  bool rejected=false;
  try {burn_and_transport(inconsistent,inconsistent,nuclear,regions,D,2);}
  catch(const std::invalid_argument&) {rejected=true;}
  check(rejected,"diffusion rejects different metal inventories across separate regions");
  nuclear.rate=.2;
  const auto reactive=burn_and_transport(m,m,nuclear,regions,{1e30,1e30},2);
  const double mean=before/m.M,exact=2*mean/(1+std::sqrt(1+4*2*nuclear.rate*mean));
  check(std::abs(reactive[0].X[0]-exact)<1e-13,"stiff diffusion plus nonlinear burning matches exact implicit uniform solution",reactive[0].X[0]-exact);
  const auto local=burn_and_transport(m,m,nuclear,regions,D,2);
  double conservation=0;
  for(int i=0;i<3;++i)conservation+=w[i]*(local[i].X[0]-m.comp[i].X[0]+2*nuclear.rate*local[i].X[0]*local[i].X[0]);
  check(std::abs(conservation)<1e-12,"burning and finite diffusion jointly conserve integrated nuclear conversion",conservation);
  PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);
  for(int i=0;i<3;++i){m.y[i].lnT=std::log(4.5e6+.5e6*i);m.y[i].lnrho=std::log(300.);}
  const auto pp_diffuse=burn_and_transport(m,m,pp,regions,{1e30,1e30},1e16);
  const auto pp_uniform=burn_and_mix(m,m,pp,{{0,3}},1e16);
  double pp_error=0;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)pp_error=std::max(pp_error,std::abs(pp_diffuse[i].X[j]-pp_uniform[i].X[j]));
  check(pp_error<3e-13,"stiff two-species pp transport matches independent instantaneous burn/mix solver",pp_error);

  Gas eos;ConstantOpacity opacity;Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux,.1,1};
  m.M=1e32;m.m={.4e32,1e32};m.comp={comp(.3),comp(.7)};m.y.resize(2);
  const double T=5e6,P=1e17,contrast=-1;
  m.y[0]={std::log(3e9),std::log(eos.rho_from_PT(T,P,m.comp[0])),std::log(T),0};
  m.y[1]={std::log(4e9),std::log(eos.rho_from_PT(T,P*std::exp(contrast),m.comp[1])),std::log(T),0};
  auto luminosity=[&](double rad) {
    const double Pb=.5*(eos.eval(T,m.rho(0),m.comp[0]).P+eos.eval(T,m.rho(1),m.comp[1]).P);
    const double L=rad*16*M_PI*constants::a_rad*constants::c*constants::G*.5*(m.m[0]+m.m[1])*std::pow(T,4)/(3*Pb);
    m.y[0].L=m.y[1].L=L;
  };
  luminosity(.5);auto coefficients=secular_mixing_diffusivities(m,physics);
  check(coefficients[0]>0 && convective_mixing_regions(m,physics).size()==2,"semiconvection acts between Schwarzschild and Ledoux thresholds");
  luminosity(2);check(secular_mixing_diffusivities(m,physics)[0]==0,"MLT-unstable face excludes secular transport");
  std::swap(m.comp[0],m.comp[1]);
  m.y[0].lnrho=std::log(eos.rho_from_PT(T,P,m.comp[0]));m.y[1].lnrho=std::log(eos.rho_from_PT(T,P*std::exp(contrast),m.comp[1]));
  luminosity(-.1);coefficients=secular_mixing_diffusivities(m,physics);
  check(coefficients[0]>0,"thermohaline mixing acts for thermally stable inverted composition");
  physics.alpha_semiconvection=physics.alpha_thermohaline=0;
  check(secular_mixing_diffusivities(m,physics)[0]==0,"zero efficiencies retain explicit no-secular-mixing control");
  return failures?1:0;
}

#include "microscopic_test_transport.hpp"
#include "ember/structure.hpp"
#include "ember/convection.hpp"
#include "ember/energy_grid.hpp"
#include "ember/evolution.hpp"
#include "ember/constants.hpp"
#include "ember/conduction.hpp"
#include "ember/evolution_proxies.hpp"
#include "ember/opacity_blend.hpp"
#include "../src/thermal_transport.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace ember;
namespace {
int checks=0;double max_jacobian=0,max_flux=0;
void require(bool p,const char* label){++checks;if(!p)throw std::runtime_error(label);}
struct Ideal final:Eos {
  EosState eval(double T,double rho,const Composition&)const override {
    EosState e;e.P=constants::R_gas*rho*T/.6;e.E=1.5*constants::R_gas*T/.6;
    e.cv=1.5*constants::R_gas/.6;e.cp=2.5*constants::R_gas/.6;
    e.chiT=e.chiRho=e.delta=1;e.grad_ad=.4;return e;
  }
  EosResponse eval_with_derivatives(double T,double rho,const Composition& c)const override {
    EosResponse e;e.state=eval(T,rho,c);return e;
  }
  const char* name()const override{return "ideal structure test";}
};
struct Radiation final:Opacity {
  OpacityState eval(double T,double rho,const Composition&)const override {
    return {20*std::pow(T/1e7,-.8)*std::pow(rho/1e3,.3),-.8,.3};
  }
  const char* name()const override{return "analytic radiative opacity";}
};
struct NoBurn final:Nuclear {
  NuclearState eval(double,double,const Composition&)const override{return {};}
  const char* name()const override{return "no burning";}
};
Model fixture() {
  Model m;m.M=2e32;m.m={1e32,2e32};
  m.y={{std::log(1e9),std::log(1e3),std::log(1e7),1e31},
       {std::log(1.1e9),std::log(9e2),std::log(9e6),1e31}};
  auto c=solar_scaled(.15,.02);c.basis=AbundanceBasis::baryon_mass;c.X[1]=.002;c.X[2]-=.002;
  m.comp={c,c};m.comp[1].X[0]+=.003;m.comp[1].X[2]-=.003;return m;
}
void jacobian(const Model& m,const Physics& p) {
  const auto a=zone_residual(m,0,p,0),b=zone_residual_numerical(m,0,p,0,nullptr,2e-5);
  const double Lunit=std::max({std::abs(m.y[0].L),std::abs(m.y[1].L),1e20});
  for(std::size_t row=0;row<4;++row) {
    double error=0,scale=0;
    for(std::size_t end=0;end<2;++end)for(std::size_t col=0;col<4;++col) {
      const double unit=col==3?Lunit:1.;
      const double actual=(end?a.dfdy_hi:a.dfdy_lo)[row][col]*unit;
      const double expected=(end?b.dfdy_hi:b.dfdy_lo)[row][col]*unit;
      error=std::max(error,std::abs(actual-expected));scale=std::max({scale,std::abs(actual),std::abs(expected)});
    }
    max_jacobian=std::max(max_jacobian,error/scale);require(error/scale<2e-6,"thermal structure Jacobian");
  }
}
}
int main() {
  try {
    Ideal eos;Radiation opacity;NoBurn nuclear;AnalyticMicroscopicTransport micro;
    Physics p{&eos,&opacity,&nuclear,1.9};p.microscopic=&micro;
    for(auto grid:{LuminosityGrid::mass_nodes,LuminosityGrid::volume_faces})
    for(double contrast:{-.2,-1e-4,-1e-10,0.,1e-10,.1})for(double target:{-.1,.2,.39,.41,2.,1e4}) {
      auto m=fixture();m.luminosity_grid=grid;m.y[1].lnT=m.y[0].lnT+contrast;
      const auto a=eos.eval(m.T(0),m.rho(0),m.comp[0]),b=eos.eval(m.T(1),m.rho(1),m.comp[1]);
      const double T=.5*(m.T(0)+m.T(1)),rho=.5*(m.rho(0)+m.rho(1)),P=.5*(a.P+b.P),mass=.5*(m.m[0]+m.m[1]);
      const double k=.5*(opacity.eval(m.T(0),m.rho(0),m.comp[0]).kappa+opacity.eval(m.T(1),m.rho(1),m.comp[1]).kappa);
      const auto f=microscopic_face(micro,0,m.m[0],m.m[1],m.y[0],m.comp[0],m.y[1],m.comp[1],true);
      const double thermal=4*constants::a_rad*constants::c*T*T*T/(3*rho*k)+f.conductivity;
      const double logmean=contrast==0?m.T(0):m.T(0)*std::expm1(contrast)/contrast;
      const double luminosity=f.carried_luminosity+target*4*M_PI*constants::G*mass*rho*thermal*logmean/P;
      m.y[0].L=m.y[1].L=luminosity;
      // An unrelated outer-face luminosity must not change this face's
      // gradient, convective classification or associated Jacobian.
      if(grid==LuminosityGrid::volume_faces)m.y[1].L=-3*luminosity;
      jacobian(m,p);
      const auto regions=convective_mixing_regions(m,p);
      require((regions.size()==1)==(target>.4),"convection partition disagrees with thermal equation");
      const auto finite=convective_mixing_faces(m,p).front();
      require((finite.diffusivity>0)==(target>.4),"finite mixing disagrees with thermal stability");
      require(finite.diffusivity==finite.length*finite.velocity/3,"finite convection diffusion/velocity relation");
      if(target>.4) {
        const auto equation=zone_equations(m,0,p,0);
        const double grad=(m.y[1].lnT-m.y[0].lnT-(m.m[1]-m.m[0])*equation[3])/(std::log(b.P)-std::log(a.P));
        EosState mid{};mid.P=P;mid.cp=.5*(a.cp+b.cp);mid.delta=.5*(a.delta+b.delta);
        const auto tr=detail::thermal_transport<0>(T,rho,P,mass,k,luminosity,true,f.carried_luminosity,
            f.conductivity,m.y[0].lnT,m.y[1].lnT);
        const double r=.5*(m.r(0)+m.r(1)),gravity=constants::G*mass/(r*r);
        const double U=mixing_length_U(T,rho,tr.opacity.value,gravity,mid,p.alpha_mlt);
        // Independent flux relation: Fconv/Fgradient = 9 q^3/(8 U).
        const double recovered_q=std::cbrt(8*U*(target-grad)/9);
        require(std::abs(recovered_q*recovered_q/finite.buoyancy_contrast-1)<2e-8,
            "mixing velocity does not correspond to the thermal equation's convective flux");
      }
      auto baseline=p;baseline.microscopic=nullptr;
      const auto full=zone_residual(m,0,p,0),plain=zone_residual(m,0,baseline,0);
      for(std::size_t row=0;row<3;++row) {
        require(full.f[row]==plain.f[row],"microscopic heat altered a nonthermal equation");
        require(full.dfdy_lo[row]==plain.dfdy_lo[row] && full.dfdy_hi[row]==plain.dfdy_hi[row],"nonthermal derivative changed");
      }
      if(target<.4) {
        const double dP=std::log(b.P)-std::log(a.P);
        const double recovered=(m.y[1].lnT-m.y[0].lnT-(m.m[1]-m.m[0])*full.f[3])/dP;
        max_flux=std::max(max_flux,std::abs(recovered-target));
        require(std::abs(recovered-target)<2e-12,"radiative gradient does not reproduce prescribed flux");
      }
    }
    // Composition derivatives include changing mobility, not just driving-force derivatives.
    const auto m=fixture();const auto f=microscopic_face(micro,0,m.m[0],m.m[1],m.y[0],m.comp[0],m.y[1],m.comp[1],true);
    for(std::size_t end=0;end<2;++end)for(std::size_t col=0;col<2;++col) {
      auto a=m.comp[0],b=m.comp[1];auto& c=end?b:a;const double h=1e-6;c.X[col]+=h;c.X[2]-=h;
      const auto plus=micro.eval(0,m.m[0],m.m[1],m.y[0],a,m.y[1],b,false);
      c.X[col]-=2*h;c.X[2]+=2*h;
      const auto minus=micro.eval(0,m.m[0],m.m[1],m.y[0],a,m.y[1],b,false);
      for(std::size_t row=0;row<2;++row) {
        const double derivative=(end?f.species.dright:f.species.dleft)[row][col];
        require(std::abs((plus.species.rate[row]-minus.species.rate[row])/(2*h)-derivative)
            <1e-8*std::max(std::abs(derivative),1.),"changing species mobility derivative");
      }
    }
    struct Cond final:Conduction {
      OpacityState eval(double,double,const Composition&)const override{return {1,0,0};}
      const char* name()const override{return "test conductor";}
    };
    auto rad=std::make_shared<Radiation>();CombinedOpacity combined(rad,std::make_shared<Cond>());
    BlendedOpacity wrapped(*rad,combined);NominalAbundanceOpacity proxy(wrapped);
    for(const Opacity* op:{static_cast<const Opacity*>(&combined),static_cast<const Opacity*>(&proxy)}) {
      auto bad=p;bad.opacity=op;bool rejected=false;
      try{zone_residual(m,0,bad,0);}catch(const std::invalid_argument&){rejected=true;}
      require(rejected,"combined conductivity was counted twice");
    }
    std::cout<<checks<<" microscopic structure checks passed; maximum Jacobian error "<<max_jacobian
             <<", radiative gradient error "<<max_flux<<'\n';return 0;
  }catch(const std::exception& e){std::cerr<<"after "<<checks<<" checks: "<<e.what()<<'\n';return 1;}
}

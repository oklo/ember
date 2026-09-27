#include "ember/envelope_transport.hpp"
#include "ember/convective_evolution_checks.hpp"
#include "ember/constants.hpp"
#include <iostream>
using namespace ember;
namespace {
class Source final:public MetalMicroscopicTransport {
public:
  explicit Source(double s,bool hot):scale(s),hot_(hot){}
  double scale;
  mutable int heat_calls{};
  bool requires_positive_species_guess()const override{return hot_;}
  MetalMicroscopicFaceResponse metal_eval(std::size_t,double,double,const Point& p,
      const Composition& a,const Point& q,const Composition& b,bool)const override {
    if(!hot_ || p.lnT<std::log(2e6) || q.lnT<std::log(2e6))
      throw std::domain_error("test microscopic species domain");
    MetalMicroscopicFaceResponse r;
    r.species.rate={a.X[0]-b.X[0],a.X[1]-b.X[1],a.Z()-b.Z()};return r;
  }
  MicroscopicHeatResponse heat(std::size_t i,double a,double b,const Point& p,
      const Composition& x,const Point& q,const Composition& y,bool d)const override {
    return heat_with_total_metal_rate(i,a,b,p,x,q,y,{},d);
  }
  MicroscopicHeatResponse heat_with_total_metal_rate(std::size_t,double,double,const Point& p,
      const Composition&,const Point& q,const Composition&,const MetalSpeciesVector& rate,bool d)const override {
    ++heat_calls;
    if(hot_ && (p.lnT<std::log(2e6) || q.lnT<std::log(2e6)))throw std::domain_error("test microscopic heat domain");
    MicroscopicHeatResponse r;
    const double h=scale*std::exp(.3*p.lnT+.2*q.lnT+.1*p.lnrho);
    r.total_rate_enthalpy={h,2*h};r.total_metal_rate_enthalpy=3*h;
    r.carried_luminosity=h*(rate[0]+2*rate[1]+3*rate[2])+(hot_?7.:0.);
    r.conductivity=scale*std::exp(.6*p.lnT+.4*q.lnT-.2*q.lnrho);
    if(d) {
      const double variable=r.carried_luminosity-(hot_?7.:0.);
      r.dcarried_lo[1]=.1*variable;r.dcarried_lo[2]=.3*variable;r.dcarried_hi[2]=.2*variable;
      r.dconductivity_lo[2]=.6*r.conductivity;r.dconductivity_hi[2]=.4*r.conductivity;
      r.dconductivity_hi[1]=-.2*r.conductivity;
    }
    return r;
  }
  const char* name()const override{return "analytic transport test";}
private:bool hot_;
};
template<class F> bool rejects(F f){try{f();}catch(const std::exception&){return true;}return false;}
class ConstantOpacity final:public Opacity {
public:
  OpacityState eval(double,double,const Composition&)const override{return {1,0,0};}
  const char* name()const override{return "constant test opacity";}
};
}
int main() {
  int failures=0;auto check=[&](bool ok,const char* text){if(!ok){++failures;std::cerr<<"FAIL "<<text<<'\n';}};
  Source cold(1,false),hot(5,true);EnvelopeTransport join(cold,hot,2e6,3e6);
  auto x=solar_scaled(.6,.02),y=solar_scaled(.5,.03);x.X[1]=y.X[1]=0;
  const MetalSpeciesVector rate{3,-2,.4};double worst=0;
  for(double T:{1e6,2e6,2.05e6,2.5e6,2.95e6,3e6,4e6}) {
    Point p{std::log(1e9),std::log(300.),std::log(T),1e30},q=p;q.lnT+=.01;
    auto eval=[&](const Point& a,const Point& b,bool d){return join.heat_with_total_metal_rate(0,1,2,a,x,b,y,rate,d);};
    hot.heat_calls=0;cold.heat_calls=0;
    const auto r=eval(p,q,true);
    if(T<=2e6)check(hot.heat_calls==0,"cool envelope does not query unsupported microscopic heat");
    if(T>=3e6)check(cold.heat_calls==0,"hot endpoint exactly selects microscopic heat");
    if(T<=2e6 || T>=3e6) {
      const auto exact=(T<=2e6?cold:hot).heat_with_total_metal_rate(0,1,2,p,x,q,y,rate,true);
      check(r.carried_luminosity==exact.carried_luminosity && r.conductivity==exact.conductivity
          && r.dcarried_lo==exact.dcarried_lo && r.dconductivity_hi==exact.dconductivity_hi,"exact endpoint values and derivatives");
    }
    for(std::size_t side=0;side<2;++side)for(std::size_t k=0;k<NVAR;++k) {
      auto a=p,b=q,c=p,e=q;const double eps=2e-6;
      (side?b:a)[static_cast<Var>(k)]+=eps;(side?e:c)[static_cast<Var>(k)]-=eps;
      const auto up=eval(a,b,false),dn=eval(c,e,false);
      const double dQ=(up.carried_luminosity-dn.carried_luminosity)/(2*eps),dK=(up.conductivity-dn.conductivity)/(2*eps);
      worst=std::max({worst,std::abs(dQ-(side?r.dcarried_hi[k]:r.dcarried_lo[k]))/std::max(1.,std::abs(r.carried_luminosity)),
        std::abs(dK-(side?r.dconductivity_hi[k]:r.dconductivity_lo[k]))/r.conductivity});
    }
    for(std::size_t k=0;k<3;++k) {
      auto up=rate,dn=rate;up[k]+=.001;dn[k]-=.001;
      const auto a=join.heat_with_total_metal_rate(0,1,2,p,x,q,y,up,false),b=join.heat_with_total_metal_rate(0,1,2,p,x,q,y,dn,false);
      const double exact=k<2?r.total_rate_enthalpy[k]:r.total_metal_rate_enthalpy;
      check(std::abs((a.carried_luminosity-b.carried_luminosity)/.002/exact-1)<1e-10,"total isotope and metal enthalpy rates retained");
    }
    if(T<2e6)check(rejects([&]{join.metal_eval(0,1,2,p,x,q,y,true);}),"cool radiative species exchange rejects");
    else check(join.metal_eval(0,1,2,p,x,q,y,true).species.rate==hot.metal_eval(0,1,2,p,x,q,y,true).species.rate,
        "microscopic species rate is not attenuated by heat overlap");
    const std::vector<MetalSpeciesVector> rates{rate};join.diagnostic_rates=rates;
    const auto diagnostic=join.heat(0,1,2,p,x,q,y,true);join.diagnostic_rates={};
    check(diagnostic.carried_luminosity==r.carried_luminosity,"accepted-model assessment uses total mixing and diffusion rates");
  }
  check(worst<1e-7,"heat Jacobian includes both changing temperature weights");
  check(join.requires_positive_species_guess(),"implicit species guess capability retained");
  check(rejects([&]{EnvelopeTransport bad(cold,hot,3e6,2e6);}),"invalid overlap rejected");
  // A radiative hot core attached to a cool mixed envelope. Set each face
  // luminosity from a specified diffusion gradient, independently of the
  // assessment's convection partition. This tests the newly supported phase.
  IdealEos eos;ConstantOpacity opacity;Physics physics{&eos,&opacity,nullptr,1.9,ConvectiveCriterion::ledoux};
  physics.microscopic=&join;
  Model m;m.M=1e32;m.m={2e31,5e31,1e32};m.luminosity_grid=LuminosityGrid::volume_faces;
  x=solar_scaled(.5,.02);x.basis=AbundanceBasis::baryon_mass;x.metal_inventory=MetalInventory::gs98;
  y=solar_scaled(.6,.02);y.basis=x.basis;y.metal_inventory=x.metal_inventory;
  m.comp={x,y,y};m.y={{std::log(1e9),std::log(200.),std::log(5e6),0},
      {std::log(2e9),std::log(100.),std::log(3e6),0},{std::log(3e9),std::log(50.),std::log(1e6),0}};
  std::vector<MetalSpeciesVector> rates(2);join.diagnostic_rates=rates;
  for(std::size_t i=0;i<2;++i) {
    const double T=.5*(m.T(i)+m.T(i+1)),rho=.5*(m.rho(i)+m.rho(i+1)),mass=.5*(m.m[i]+m.m[i+1]);
    const double P=.5*(eos.eval(m.T(i),m.rho(i),m.comp[i]).P+eos.eval(m.T(i+1),m.rho(i+1),m.comp[i+1]).P);
    const auto h=join.heat(i,m.m[i],m.m[i+1],m.y[i],m.comp[i],m.y[i+1],m.comp[i+1],false);
    const double K=4*constants::a_rad*constants::c*T*T*T/(3*rho)+h.conductivity;
    const double Tlog=(m.T(i+1)-m.T(i))/(m.y[i+1].lnT-m.y[i].lnT);
    m.y[i].L=h.carried_luminosity+(i?.8:.2)*4*M_PI*constants::G*mass*rho*K*Tlog/P;
  }
  m.y.back().L=m.y[1].L;
  check(convective_mixing_regions(m,physics)==MixingRegions{{0,1},{1,3}},"control has a radiative core and mixed cool envelope");
  const auto assessed=driver::check_envelope_transport(m,physics,join,rates,2e6,.01);
  const auto weights=nodal_mass_weights(m);
  check(assessed.radiative_boundaries==1 && std::abs(assessed.convective_mass_fraction-(weights[1]+weights[2])/m.M)<1e-15,
      "assessment accepts a hot radiative boundary and measures the convective mass");
  auto cool_core=m;cool_core.y[0].lnT=std::log(1.9e6);cool_core.y[0].L*=.001;
  check(convective_mixing_regions(cool_core,physics)==MixingRegions{{0,1},{1,3}},
      "cold-boundary negative control retains its radiative core");
  bool domain_rejected=false;
  try{driver::check_envelope_transport(cool_core,physics,join,rates,2e6,.01);}
  catch(const std::domain_error& e){domain_rejected=std::string(e.what()).find("microscopic species domain")!=std::string::npos;}
  check(domain_rejected,
      "unsupported cold radiative boundary cannot silently stop settling");
  EvolutionOptions finite;finite.convective_mixing=ConvectiveMixing::finite_implicit;
  finite.instantaneous_mixing_below_T=2e6;
  check(instantaneous_mixing_regions(m,physics,finite)==MixingRegions{{0,1},{1,3}},
      "cool convection is collapsed independently of a hot radiative interface");
  check(instantaneous_mixing_regions(cool_core,physics,finite)==MixingRegions{{0,1},{1,3}},
      "cool radiative interface must remain exposed in finite mode");
  check(rejects([&]{driver::check_envelope_transport(cool_core,physics,join,rates,2e6,.01,finite);}),
      "finite mode must also reject unsupported cool radiative exchange");
  finite.instantaneous_mixing_below_T=m.T(2);
  check(instantaneous_mixing_regions(m,physics,finite)==MixingRegions{{0,1},{1,2},{2,3}},
      "threshold equality leaves a face finite");
  finite.instantaneous_mixing_below_T=std::nextafter(m.T(2),INFINITY);
  check(instantaneous_mixing_regions(m,physics,finite)==MixingRegions{{0,1},{1,3}},
      "the current thermal state must determine the instantaneous partition");
  auto unmixed=m;unmixed.comp[2].X[0]-=1e-8;unmixed.comp[2].X[2]+=1e-8;
  check(rejects([&]{driver::check_envelope_transport(unmixed,physics,join,rates,2e6,.01);}),
      "instantaneous mixing cannot accept a residual envelope gradient");
  check(rejects([&]{driver::check_envelope_transport(m,physics,join,{},2e6,.01);}),
      "missing total composition rates are rejected");
  join.diagnostic_rates={};
  std::cout<<"maximum scaled derivative error "<<worst<<", failures "<<failures<<'\n';
  return failures?1:0;
}

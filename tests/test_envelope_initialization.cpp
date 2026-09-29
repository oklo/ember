#include "ember/stellar_seed.hpp"
#include "ember/envelope_atmosphere.hpp"
#include "ember/eos.hpp"
#include "ember/constants.hpp"
#include <iostream>

using namespace ember;
namespace {
struct Heating final : Nuclear {
  NuclearState eval(double,double,const Composition&) const override { NuclearState s{};s.eps=1;return s; }
  const char* name()const override{return "constant heating";}
};
struct Boundary final : Atmosphere {
  AtmosphereState eval(double,double,const Composition&) const override {
    AtmosphereState s{};s.T=3e4;s.rho=1e-5;s.P=1e8;return s;
  }
  const char* name()const override{return "fixed test boundary";}
};
struct Gas final : Eos {
  mutable Composition last;
  EosState eval(double T,double rho,const Composition& c)const override {
    last=c;const double R=constants::R_gas*(c.mu_ions_inv()+c.mu_elec_inv());
    EosState s{};s.P=rho*R*T;s.cp=2.5*R;s.cv=1.5*R;s.delta=s.chiRho=s.chiT=1;
    s.grad_ad=.4;s.E=s.cv*T;return s;
  }
  const char* name()const override{return "ideal gas for envelope limit";}
};
struct PureGas final : Eos {
  EosState eval(double T,double rho,const Composition& c)const override {return Gas{}.eval(T,rho,c);}
  const char* name()const override{return "thread-safe ideal gas";}
};
struct Transparent final : Opacity {
  OpacityState eval(double,double,const Composition&)const override {OpacityState s{};s.kappa=1e-20;return s;}
  const char* name()const override{return "radiatively isothermal limit";}
};
struct NarrowBoundary final : Atmosphere {
  mutable int refused{};
  AtmosphereState eval(double Teff,double g,const Composition& c)const override {
    if(std::abs(Teff/3000-1)>.002){++refused;throw std::domain_error("temperature outside test atmosphere");}
    return Boundary{}.eval(Teff,g,c);
  }
  const char* name()const override{return "narrow atmosphere for restart test";}
};
}
int main() {
  int failed=0;
  auto check=[&](bool ok,const char* name,double value=0.) {
    std::cout<<(ok?"PASS ":"FAIL ")<<name<<" "<<value<<'\n';failed+=!ok;
  };
  try {
    Heating nuclear;Boundary boundary;auto c=solar_scaled(.7,.02);
    for(double fraction:{0.,.0015,.01}) {
      const double M=.1*constants::Msun;
      const auto m=example::stellar_seed(128,M,1.6*constants::Rsun,c,nuclear,boundary,
          1.5,2857,LuminosityGrid::volume_faces,fraction*M);
      check(m.valid_outer_mass(),"initial mesh and outer reservoir have the declared mass");
      double total=0,hydrogen=0,largest=0;
      for(std::size_t i=0;i<m.size();++i) {
        const double lo=i==0?0:.5*(m.m[i-1]+m.m[i]);
        const double hi=i+1==m.size()?M:.5*(m.m[i]+m.m[i+1]);
        total+=hi-lo;hydrogen+=(hi-lo)*m.comp[i].h1();
        if(i) {
          const double rb=.5*(m.r(i)+m.r(i-1)),rho=.5*(m.rho(i)+m.rho(i-1));
          const double integral=4*M_PI*rb*rb*rb*rho*(m.y[i].lnr-m.y[i-1].lnr);
          largest=std::max(largest,std::abs(integral-(m.m[i]-m.m[i-1]))/M);
        }
      }
      check(std::abs(total/M-1)<1e-14 && std::abs(hydrogen/M-c.h1())<1e-14,
            "whole-star mass and hydrogen include the outer reservoir");
      check(largest<1e-14,"seed still satisfies the discrete mass equation",largest);
    }
    bool rejected=false;
    try{example::stellar_seed(128,1e32,1e11,c,nuclear,boundary,1.5,3000,LuminosityGrid::mass_nodes,1e29);}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"reservoir requires conservative face luminosities");

    Gas gas;Transparent opacity;
    c.basis=AbundanceBasis::baryon_mass;c.X[0]-=2e-5;c[Species::H2]=2e-5;
    const double M=.1*constants::Msun,R=.1*constants::Rsun,Teff=3000;
    EnvelopeAtmosphere envelope(boundary,opacity,gas,1.9,M,.0015*M,20);
    const auto forward=envelope.from_photosphere(Teff,R,c);
    check(gas.last==c,"native layer uses the actual interior composition, including D and metals");
    check(std::abs(forward.T_base/3e4-1)<1e-9 && forward.P_base>forward.P100
          && forward.r_base<R,"transparent layer is isothermal and hydrostatically ordered");
    const double pseudo=Teff*std::sqrt(R/forward.r_base),g=constants::G*M/(forward.r_base*forward.r_base);
    const auto inverse=envelope.photosphere(pseudo,g,c);
    check(std::abs(inverse.R/R-1)<1e-9 && std::abs(inverse.Teff/Teff-1)<1e-9,
          "initial photosphere and base-boundary evaluations recover the same physical surface");
    NarrowBoundary narrow;
    EnvelopeAtmosphere restarted(narrow,opacity,gas,1.9,M,.0015*M,20);
    const auto cold_start=restarted.photosphere(pseudo,g,c);
    check(narrow.refused>0 && std::abs(cold_start.R/R-1)<1e-9,
          "restart finds a supported surface even when the base-radius trial is outside the atmosphere");
    const auto state=envelope.eval(pseudo,g,c);
    const auto value_start=envelope.integrations();
    const auto values=envelope.eval_value(pseudo,g,c);
    const auto value_count=envelope.integrations()-value_start;
    const auto full_start=envelope.integrations();
    envelope.eval(pseudo,g,c);
    check(envelope.integrations()-full_start>value_count,"value-only envelope avoids sensitivity integrations");
    check(std::abs(values.P/state.P-1)<1e-9 && std::abs(values.T/state.T-1)<1e-9,
          "value-only envelope preserves the physical boundary");
    PureGas shared_gas;
    EnvelopeAtmosphere parallel(boundary,opacity,shared_gas,1.9,M,.0015*M,20);
    parallel.from_photosphere(Teff,R,c);
    const auto serial=parallel.eval(pseudo,g,c);
    for(std::size_t threads:{2,4}) {
      parallel.evaluation_threads(threads);
      const auto other=parallel.eval(pseudo,g,c);
      check(std::abs(other.P/serial.P-1)<1e-9 && std::abs(other.T/serial.T-1)<1e-9
        && std::abs(other.dlnP_dlng-serial.dlnP_dlng)<1e-8
        && std::abs(other.dlnP_dlnTeff-serial.dlnP_dlnTeff)<1e-8
        && std::abs(other.dlnT_dlng-serial.dlnT_dlng)<1e-8
        && std::abs(other.dlnT_dlnTeff-serial.dlnT_dlnTeff)<1e-8,
        "parallel boundary values and derivatives match serial");
    }
    constexpr double difference=5e-5;
    const auto plus_g=envelope.photosphere(pseudo,g*std::exp(difference),c);
    const auto minus_g=envelope.photosphere(pseudo,g*std::exp(-difference),c);
    const double slope=(std::log(plus_g.P_base)-std::log(minus_g.P_base))/(2*difference);
    check(std::abs(slope-state.dlnP_dlng)<1e-5,"native envelope pressure derivative agrees at a second increment",slope-state.dlnP_dlng);

    parallel.jacobian_reuse(.001);
    const auto anchor=parallel.eval(pseudo,g,c);
    const auto computations=parallel.jacobian_computed();
    const auto nearby=parallel.eval(pseudo*std::exp(2e-5),g*std::exp(3e-5),c);
    const auto exact=envelope.eval(pseudo*std::exp(2e-5),g*std::exp(3e-5),c);
    check(parallel.jacobian_computed()==computations && parallel.jacobian_reused()==1,
          "nearby Newton iteration reuses the boundary derivatives");
    check(std::abs(nearby.P/exact.P-1)<1e-9 && std::abs(nearby.T/exact.T-1)<1e-9
          && std::abs(nearby.P/anchor.P-1)>1e-7,
          "reused Jacobian does not reuse boundary values");
    for(int i=0;i<7;++i)parallel.eval(pseudo,g,c);
    check(parallel.jacobian_computed()==computations+1,"boundary Jacobian refreshes periodically");
    const auto refreshed=parallel.jacobian_computed();
    parallel.eval(pseudo*std::exp(.002),g,c);
    check(parallel.jacobian_computed()==refreshed+1,"structural displacement refreshes boundary Jacobian");
    auto changed=c;changed.X[0]-=2e-5;changed.X[2]+=2e-5;
    parallel.eval(pseudo*std::exp(.002),g,changed);
    check(parallel.jacobian_computed()==refreshed+2,"composition displacement refreshes boundary Jacobian");
    restarted.jacobian_reuse(.001);
    restarted.eval(pseudo,g,c);
    rejected=false;
    try{restarted.eval(pseudo*1.1,g,c);}catch(const std::domain_error&){rejected=true;}
    check(rejected,"Jacobian reuse never bypasses boundary-value domain checks");

  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
  return failed?1:0;
}

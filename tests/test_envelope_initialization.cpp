#include "ember/stellar_seed.hpp"
#include "ember/envelope_atmosphere.hpp"
#include "ember/eos.hpp"
#include "ember/constants.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iomanip>

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

    // An independent ideal-gas limit: exact source nodes in composition and
    // log-linear P(T,rho), so interpolated H/He has an analytic mixture answer.
    const auto file=std::filesystem::temp_directory_path()/"ember-envelope-ideal-source.txt";
    struct Remove {std::filesystem::path p;~Remove(){std::filesystem::remove(p);}} remove{file};
    const std::array<double,3> xs{.6,.7,.8};const std::array<double,2> ys{0,.02};
    std::array<double,6> lt{},lr{};
    for(int i=0;i<6;++i){lt[i]=std::log(1e2)+i*std::log(10.);lr[i]=std::log(1e-10)+i*std::log(100.);}
    {std::ofstream out(file);out<<std::setprecision(17)<<"EMBER_ENVELOPE_SOURCE_V1\n3 2 6 6\n";
      for(auto axis:{std::vector<double>(xs.begin(),xs.end()),std::vector<double>(ys.begin(),ys.end()),
                     std::vector<double>(lt.begin(),lt.end()),std::vector<double>(lr.begin(),lr.end())}) {
        for(double v:axis)out<<v<<' ';out<<'\n';
      }
      for(double X:xs)for(double Y:ys) {
        const double gas_constant=constants::R_gas*(2*X+Y+.75*(1-X-Y));
        for(int q=0;q<5;++q)for(double t:lt)for(double r:lr)
          out<<(q==0?t+r+std::log(gas_constant):q==1?.4:q==2?2.5*gas_constant:1)<<'\n';
      }
    }
    EnvelopeSource source(file.string());
    Composition mixture{};mixture.basis=AbundanceBasis::baryon_mass;
    mixture[Species::H1]=.98*.7;mixture[Species::He3]=.98*.02;
    mixture[Species::He4]=.98*.28;mixture[Species::O16]=.02;
    const double T=1e4,P=1e10,Rh=constants::R_gas*(2*.7+.02+.75*.28);
    for(auto mode:{EnvelopeMetals::neutral,EnvelopeMetals::ionized}) {
      double guess=std::log(.01);const auto s=source.at_pressure(std::log(T),std::log(P),mixture,guess,mode);
      const double gas_constant=.98*Rh+.02*constants::R_gas*(mode==EnvelopeMetals::neutral?1:9)/16;
      check(std::abs(s.rho/(P/(gas_constant*T))-1)<1e-7 && std::abs(s.cp/(2.5*gas_constant)-1)<1e-7,
            "metal mixture recovers independent ideal-gas density and heat capacity");
      check(std::abs(s.grad_ad-.4)<1e-7 && std::abs(s.chiRho-1)<1e-7,
            "metal mixture recovers ideal-gas convection and compressibility");
      constexpr double h=1e-5;double gp=guess,gm=guess;
      const auto plus=source.at_pressure(std::log(T)+h,std::log(P),mixture,gp,mode);
      const auto minus=source.at_pressure(std::log(T)-h,std::log(P),mixture,gm,mode);
      const double fd=-(std::log(plus.rho)-std::log(minus.rho))/(2*h);
      check(std::abs(fd-s.delta)<1e-7,"thermal expansion agrees with density differences",fd-s.delta);
      gp=gm=guess;
      const auto pp=source.at_pressure(std::log(T),std::log(P)+h,mixture,gp,mode);
      const auto pm=source.at_pressure(std::log(T),std::log(P)-h,mixture,gm,mode);
      const double fdP=(std::log(pp.rho)-std::log(pm.rho))/(2*h);
      check(std::abs(fdP-1/s.chiRho)<1e-7,"compressibility agrees with density differences",fdP-1/s.chiRho);
    }
    for(auto mode:{EnvelopeMetals::neutral,EnvelopeMetals::ionized}) {
      EnvelopeAtmosphere layer(boundary,opacity,source,1.9,M,.00015*M,20,mode);
      const double layer_radius=1.6*constants::Rsun;
      const auto surface=layer.from_photosphere(Teff,layer_radius,mixture);
      const double tb=Teff*std::sqrt(layer_radius/surface.r_base),gb=constants::G*M/(surface.r_base*surface.r_base);
      layer.evaluation_threads(1);const auto one=layer.eval(tb,gb,mixture);
      layer.evaluation_threads(2);const auto two=layer.eval(tb,gb,mixture);
      check(std::abs(one.T/two.T-1)<1e-12 && std::abs(one.dlnT_dlnTeff-two.dlnT_dlnTeff)<1e-10,
            "parallel boundary derivatives retain the selected metal approximation");
    }
    double guess=std::log(.01);rejected=false;
    try{source.at_pressure(std::log(T),std::log(P),mixture,guess,EnvelopeMetals::reject);}
    catch(const std::domain_error&){rejected=true;}
    check(rejected,"metal approximation requires explicit selection");
    // A missing source node must not be smoothed into a supported query.
    const auto masked_file=std::filesystem::path(file.string()+".masked");
    Remove remove_masked{masked_file};
    std::ifstream input(file);std::vector<std::string> lines;std::string line;
    while(std::getline(input,line))lines.push_back(line);
    lines[6+3*5*36+2*36+2*6+3]="0";
    {std::ofstream output(masked_file);for(const auto& value:lines)output<<value<<'\n';}
    EnvelopeSource masked(masked_file.string());rejected=false;
    try{masked.eval(std::log(1e4),std::log(.001),.7,.02);}
    catch(const std::domain_error&){rejected=true;}
    check(rejected,"missing envelope source response is rejected across its interpolation stencil");

  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
  return failed?1:0;
}

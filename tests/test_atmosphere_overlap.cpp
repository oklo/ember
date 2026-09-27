#include "ember/atmosphere_overlap.hpp"
#include <iomanip>
#include <iostream>
using namespace ember;
namespace {
class Gas final:public Eos {
public:
  EosState eval(double t,double rho,const Composition& c)const override {
    EosState s;const double pg=constants::R_gas*c.mu_ions_inv()*rho*t,pr=constants::a_rad*std::pow(t,4)/3;
    s.P=pg+pr;s.chiRho=pg/s.P;s.chiT=(pg+4*pr)/s.P;return s;
  }
  const char* name()const override{return "test gas";}
};
class Source final:public Atmosphere {
public:
  Source(const Eos& eos,bool main):eos_(eos),main_(main){}
  double tau=100;
  AtmosphereState eval(double t,double g,const Composition& c)const override {
    const double lg=std::log10(g);
    if(main_?lg<4.9:(lg>5.5 || c.h1()<.695))throw std::domain_error("test source unsupported");
    AtmosphereState s;const double offset=main_?.05:0;
    s.T=t*std::exp(offset)*std::pow(g/1e5,.03);
    s.Pgas=10*std::exp(2*offset)*std::pow(t/3000.,-.4)*std::pow(g/1e5,.7);
    const double pr=constants::a_rad*std::pow(s.T,4)/3;s.P=s.Pgas+pr;s.tau=tau;
    s.dlnT_dlnTeff=1;s.dlnT_dlng=.03;
    s.dlnP_dlnTeff=(-.4*s.Pgas+4*pr)/s.P;s.dlnP_dlng=(.7*s.Pgas+.12*pr)/s.P;
    s.rho=eos_.rho_from_PT(s.T,s.P,c);return s;
  }
  const char* name()const override{return "synthetic overlapping atmosphere";}
private:
  const Eos& eos_;bool main_;
};
template<class F>bool rejects(F f){try{f();}catch(const std::exception&){return true;}return false;}
}
int main() {
  int failed=0;const auto check=[&](bool ok,const char* s){std::cout<<(ok?"PASS ":"FAIL ")<<s<<'\n';failed+=!ok;};
  try {
    Gas eos;Source a(eos,false),b(eos,true);AtmosphereOverlap blend(eos,a,b,{4.95,5.10,.6985,.6995});
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;
    for(double lg:{4.,4.95,5.10,6.}) {
      const double g=std::pow(10.,lg);const auto s=blend.eval(2800,g,c),ref=(lg<5.?a:b).eval(2800,g,c);
      check(s.T==ref.T && s.P==ref.P && s.rho==ref.rho && s.dlnP_dlng==ref.dlnP_dlng,
          "outside transition exactly preserves the selected source");
    }
    double error=0;
    for(double h:{.7,.6991,.6985})for(double lg:{4.95,4.98,5.025,5.08,5.1}) {
      c=solar_scaled(h,.02);c.basis=AbundanceBasis::baryon_mass;const auto before=c;
      const double t=2800,g=std::pow(10.,lg),eps=1e-6;const auto s=blend.eval(t,g,c);
      const auto gp=blend.eval(t,g*std::exp(eps),c),gm=blend.eval(t,g*std::exp(-eps),c);
      const auto tp=blend.eval(t*std::exp(eps),g,c),tm=blend.eval(t*std::exp(-eps),g,c);
      error=std::max({error,std::abs(std::log(gp.T/gm.T)/(2*eps)-s.dlnT_dlng),
        std::abs(std::log(gp.P/gm.P)/(2*eps)-s.dlnP_dlng),
        std::abs(std::log(tp.T/tm.T)/(2*eps)-s.dlnT_dlnTeff),
        std::abs(std::log(tp.P/tm.P)/(2*eps)-s.dlnP_dlnTeff)});
      check(c==before && std::abs(eos.eval(s.T,s.rho,c).P/s.P-1)<2e-10,"inventory and actual EOS density retained");
    }
    std::cout<<"maximum derivative error "<<std::setprecision(5)<<error<<'\n';
    check(error<2e-8,"thermal and gravity derivatives include changing blend weight and radiation");
    c=solar_scaled(.6,.02);c.basis=AbundanceBasis::baryon_mass;
    const auto s=blend.eval(2800,1e5,c),m=b.eval(2800,1e5,c);
    check(s.T==m.T && s.P==m.P,"evolved composition does not return to the contraction table during expansion");
    check(rejects([&]{blend.eval(2800,1e4,c);}),"no fallback when selected source lacks coverage");
    c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;b.tau=25;
    check(rejects([&]{blend.eval(2800,std::pow(10.,5.025),c);}),"different matching depths cannot be blended");
    check(rejects([&]{AtmosphereOverlap bad(eos,a,b,{5.1,4.9,.6,.7});}),"invalid overlap rejected");
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
  return failed?1:0;
}

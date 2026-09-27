#include "ember/atmosphere_hydrogen_envelope.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <chrono>
#include "ember/runtime_identity.hpp"

using namespace ember;
namespace {
unsigned checks=0;
void require(bool ok,const char* message) {++checks;if(!ok)throw std::runtime_error(message);}
template<class F> bool rejects(F f) {try {f();return false;}catch(const std::exception&){return true;}}
bool close(double a,double b,double tol=3e-10) {
  return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});
}
Composition composition(double helium,double z,double he3=0) {
  auto c=solar_scaled(1-helium-z,z);c.X[1]=he3;c.X[2]=helium-he3;
  c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;return c;
}
class Gas final:public Eos {
public:
  EosState eval(double t,double rho,const Composition& c)const override {
    EosState s{};const double pg=constants::R_gas*c.mu_ions_inv()*rho*t;
    const double pr=constants::a_rad*std::pow(t,4)/3;
    s.P=pg+pr;s.chiT=(pg+4*pr)/s.P;s.chiRho=pg/s.P;return s;
  }
  const char* name()const override{return "analytic test gas";}
};
class Boundary final:public Atmosphere {
public:
  Boundary(const Eos& eos,bool upper):eos_(eos),upper_(upper){}
  double tau=100;
  AtmosphereState eval(double t,double g,const Composition& c)const override {
    if((upper_ && g<1e5) || (!upper_ && g>1e6))throw std::domain_error("source gravity");
    AtmosphereState s{};s.T=5000*std::pow(t/4700.,.9)*std::pow(g/1e5,.04);
    s.Pgas=10*std::pow(t/4700.,-.3)*std::pow(g/1e5,.8);
    const double pr=constants::a_rad*std::pow(s.T,4)/3;
    s.P=s.Pgas+pr;s.tau=tau;s.dlnT_dlnTeff=.9;s.dlnT_dlng=.04;
    s.dlnP_dlnTeff=(-.3*s.Pgas+3.6*pr)/s.P;
    s.dlnP_dlng=(.8*s.Pgas+.16*pr)/s.P;
    s.rho=eos_.rho_from_PT(s.T,s.P,c);return s;
  }
  const char* name()const override{return "analytic power-law boundary";}
private:
  const Eos& eos_;bool upper_;
};
double derivative_error(const Atmosphere& source,double t,double g,const Composition& c) {
  const auto s=source.eval(t,g,c);double maximum=0;
  for(int axis=0;axis<2;++axis) {
    constexpr double step=1e-5;double tv[4],pv[4];int n=0;
    for(int sign:{-2,-1,1,2}) {
      const auto q=source.eval(t*std::exp(axis==0?sign*step:0),g*std::exp(axis==1?sign*step:0),c);
      tv[n]=std::log(q.T);pv[n++]=std::log(q.P);
    }
    for(int field=0;field<2;++field) {
      const double* v=field?pv:tv;const double fd=(v[0]-8*v[1]+8*v[2]-v[3])/(12*step);
      const double exact=field?(axis?s.dlnP_dlng:s.dlnP_dlnTeff):(axis?s.dlnT_dlng:s.dlnT_dlnTeff);
      maximum=std::max(maximum,std::abs(fd-exact));
    }
  }
  return maximum;
}
std::string table() {
  const auto metals=solar_scaled(.7,.02);
  std::ostringstream out;out<<std::setprecision(17);
  out<<"EMBER_HYDROGEN_DOMINATED_ATMOSPHERE 1\nsource \"analytic test grid\"\n"
       "approximation \"trace helium neglected\"\nbasis baryon_mass\ntau 100\n"
       "maximum_helium .001\nsource_helium 0\nmetal_pattern";
  for(std::size_t k=0;k<NMETALS;++k)out<<' '<<metals.X[METAL_BEGIN+k]/metals.Z();
  out<<"\nmetallicity 2 0 .001\nlog_teff 2 3.5 4\nlog_g 2 5 6\ndata\n";
  for(double z:{0.,.001})for(double lt:{3.5,4.})for(double lg:{5.,6.})
    out<<.9*lt+.04*lg+.25+3*z<<' '<<-.3*lt+.8*lg-2+20*z<<'\n';
  return out.str();
}
}
int main() {
  try {
    Gas eos;Boundary a(eos,false),b(eos,true);
    GravityIntervalAtmosphere gravity(eos,a,b,5.2,5.8);
    const auto c=composition(.0005,1e-5,.0001);
    for(double lg:{5.1,5.2,5.4,5.8,5.9}) {
      const double g=std::pow(10.,lg);const auto s=gravity.eval(4700,g,c),ref=a.eval(4700,g,c);
      require(close(s.T,ref.T) && close(s.P,ref.P),"power law not reproduced by gravity interpolation");
      require(close(eos.eval(s.T,s.rho,c).P,s.P),"actual-composition EOS pressure mismatch");
      require(derivative_error(gravity,4700,g,c)<2e-9,"gravity/thermal derivative omits radiation or endpoint slope");
    }
    require(rejects([&]{gravity.eval(4700,0,c);}),"zero gravity accepted");
    require(rejects([&]{GravityIntervalAtmosphere invalid(eos,a,b,6,5);}),"reversed interval accepted");
    b.tau=25;
    require(rejects([&]{gravity.eval(4700,3e5,c);}),"different atmosphere depths joined");
    b.tau=100;
    std::istringstream input(table());
    HydrogenDominatedAtmosphereGrid grid(eos,input,
        HydrogenDominatedAtmosphereGrid::Approximation::neglect_trace_atmospheric_helium,.001);
    const auto s=grid.eval(4700,3e5,c);
    require(close(std::log10(s.T),.9*std::log10(4700)+.04*std::log10(3e5)+.25+3*c.Z()),
        "hydrogen table interpolation differs from analytic source");
    require(derivative_error(grid,4700,3e5,c)<2e-9,"pure-H interpolation derivative");
    auto different_isotope=c;different_isotope.X[1]=0;different_isotope.X[2]=.0005;
    const auto h=grid.eval(4700,3e5,different_isotope);
    require(s.T==h.T && s.P==h.P && std::abs(s.rho/h.rho-1)>1e-8,"actual isotope lost in density inversion");
    auto deuterium=c;deuterium[Species::H2]=1e-5;deuterium[Species::H1]-=1e-5;
    require(!grid.covers(4700,3e5,deuterium),"deuterium read as a metal or silently discarded");
    require(!grid.covers(4700,3e5,composition(.002,1e-5)),"helium approximation exceeded");
    require(!grid.covers(4700,3e5,composition(.0005,.002)),"metals extrapolated");
    require(!grid.covers(4700,1e7,c),"gravity extrapolated");

    // Use retained physical sources to check the actual selection, including
    // composition joins and both gravity intervals, without a stellar solve.
    const auto path=std::filesystem::path(EMBER_SOURCE_DIR)/"data/atmosphere/lifetime_hydrogen_envelope/envelope.dat";
    HydrogenEnvelopeAtmosphere envelope(eos,a,path);
    for(double he:{.0003,.0005,.0007,.00089,.00095})
      for(double z:{2e-9,5e-8,8e-7,1e-6,3e-6,1e-5,2e-5}) {
        auto q=composition(he,z,std::min(he/4.,.0001));
        const auto v=envelope.eval(4731,std::pow(10.,5.43),q);
        require(close(eos.eval(v.T,v.rho,q).P,v.P),"composition join does not preserve EOS");
        require(derivative_error(envelope,4731,std::pow(10.,5.43),q)<2e-8,"composition join derivative");
      }
    for(double lg:{5.49,5.54,5.61,5.68,5.77,5.83,5.87,5.93,6.03}) {
      require(derivative_error(envelope,4731,std::pow(10.,lg),composition(.0001,1e-30))<2e-8,
          "retained source gravity join derivative");
    }
    for(double he:{.0004,.0006,.0009}) {
      const auto below=envelope.eval(4731,3e5,composition(he-1e-10,1.4e-6));
      const auto above=envelope.eval(4731,3e5,composition(he+1e-10,1.4e-6));
      require(std::abs(std::log(below.P/above.P))<1e-4,"discontinuous helium join");
    }
    require(rejects([&]{envelope.eval(4500,1e6,composition(.0001,1e-30));}),"cool unsupported state extrapolated");
    require(rejects([&]{envelope.eval(4731,2e6,composition(.0001,1e-30));}),"high-gravity unsupported state extrapolated");
    require(rejects([&]{envelope.eval(4731,1e6,composition(.006,1e-30));}),"large helium fraction accepted");
    const auto stamp=std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto temporary=std::filesystem::temp_directory_path()/("ember-mixed-helium-"+std::to_string(stamp));
    std::filesystem::create_directory(temporary);
    const auto mixed_path=temporary/"mixed.dat";
    {
      std::ofstream out(mixed_path);out<<std::setprecision(17)
        <<"EMBER_COMPOSITION_ATMOSPHERE 4\nsource \"analytic mixed helium\"\napproximation \"test only\"\nbasis baryon_mass\ntau 100\nmetals 0 0 0 0 0\nmetal_tolerance 0\n"
        <<"hydrogen 2 .78 .9995\nhelium3 1 0\nlog_teff 2 3.5 4\nlog_g 2 4.5 6.3\ndata\n";
      for(double hydrogen:{.78,.9995})for(double t:{3.5,4.})for(double g:{4.5,6.3})
        out<<"1 "<<.9*t+.04*g+.25+.02*(1-hydrogen)<<' '<<-.3*t+.8*g+4+.1*(1-hydrogen)<<'\n';
      out<<"cells\n1\n";
    }
    auto spec=HydrogenEnvelopeAtmosphere::read(path);
    spec.mixed_helium=mixed_path;spec.mixed_maximum_helium3=.02;spec.mixed_maximum_metals=1e-10;
    spec.mixed_helium_join={.0005,.0009};spec.mixed_metal_join={1e-12,1e-10};
    HydrogenEnvelopeAtmosphere mixed(eos,a,spec);
    for(double he:{.0001,.0005,.0007,.0009,.01,.1}) {
      const auto q=composition(he,1e-30,std::min(he/5,.009));
      const auto v=mixed.eval(4901,std::pow(10.,6.03),q);
      require(close(eos.eval(v.T,v.rho,q).P,v.P),"mixed-helium envelope lost actual-composition density");
      require(derivative_error(mixed,4901,std::pow(10.,6.03),q)<2e-8,"mixed-helium connection derivative");
      if(he<=.0005) {
        const auto old=envelope.eval(4901,std::pow(10.,6.03),q);
        require(old.T==v.T && old.P==v.P && old.rho==v.rho,"optional mixed source changed preceding models");
      }
    }
    for(double edge:{.0005,.0009}) {
      const auto l=mixed.eval(4901,1e6,composition(edge-1e-12,1e-30));
      const auto r=mixed.eval(4901,1e6,composition(edge+1e-12,1e-30));
      require(std::abs(std::log(l.P/r.P))<1e-6,"mixed-helium connection is discontinuous");
    }
    require(rejects([&]{mixed.eval(4901,1e6,composition(.03,1e-30,.021));}),"mixed envelope ignored its isotope bound");
    const auto manifest=temporary/"envelope.dat";
    {
      std::ifstream in(path);std::ofstream out(manifest);std::string line;std::getline(in,line);
      out<<"EMBER_HYDROGEN_ENVELOPE_ATMOSPHERE 2\n";
      for(int i=0;std::getline(in,line);++i) {
        if(i<4) {std::istringstream row(line);std::string key,file;row>>key>>std::quoted(file);out<<key<<' '<<std::quoted((path.parent_path()/file).string())<<'\n';}
        else out<<line<<'\n';
      }
      out<<"mixed_helium \"mixed.dat\"\nmixed_maximum_helium3 .02\nmixed_maximum_metals 1e-10\nmixed_helium_join .0005 .0009\nmixed_metal_join 1e-12 1e-10\n";
    }
    HydrogenEnvelopeAtmosphere parsed(eos,a,manifest);
    const auto q=composition(.1,1e-30,.009);
    require(parsed.eval(4901,1e6,q).P==mixed.eval(4901,1e6,q).P,"mixed-helium manifest changed its selection");
    driver::RuntimeIdentity identity;identity.hydrogen_envelope("boundary",manifest);
    require(identity.values.contains("boundary.mixed_helium")
        && identity.values.contains("boundary.mixed_helium_join.low"),"mixed source absent from restart identity");
    std::filesystem::remove_all(temporary);
    std::cout<<checks<<" atmosphere checks passed\n";
  }catch(const std::exception& e){std::cerr<<"after "<<checks<<" checks: "<<e.what()<<'\n';return 1;}
}

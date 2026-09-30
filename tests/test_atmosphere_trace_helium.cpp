#include "ember/atmosphere_trace_helium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
using namespace ember;
namespace {
int checks=0;
void require(bool ok,const char* message){++checks;if(!ok)throw std::runtime_error(message);}
bool close(double a,double b,double tol=2e-11){return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});}
template<class F>bool rejects(F f){try{f();return false;}catch(const std::exception&){return true;}}
class Gas final:public Eos{public:
 EosState eval(double t,double rho,const Composition& c)const override{
  EosState s{};double pg=constants::R_gas*c.mu_ions_inv()*rho*t,pr=constants::a_rad*std::pow(t,4)/3;
  s.P=pg+pr;s.chiT=(pg+4*pr)/s.P;s.chiRho=pg/s.P;return s;
 }
 const char* name()const override{return "analytic test gas";}
};
double lt(double h,double z,double t,double g){return .2+.95*t+.02*g+.1*h+2e4*z+3e3*h*z*t*g;}
double lp(double h,double z,double t,double g){return 5.2+.2*t+.3*g-.4*h+4e4*z-1e3*h*z*t*g;}
Composition comp(double h,double y,double z){auto c=solar_scaled(h,z);c.X[1]=y;c.X[2]-=y;c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;return c;}
std::string fixture(int missing=-1,double max_y=.0003){
 auto c=comp(.7,0,.02);std::ostringstream s;s<<std::setprecision(17);
 s<<"EMBER_TRACE_HELIUM_ATMOSPHERE 1\nsource \"synthetic analytic test; not a physical table\"\napproximation \"test of explicit atmospheric isotope approximation\"\nbasis baryon_mass\ntau 100\nmaximum_helium3 "<<max_y<<"\nmetal_pattern";
 for(size_t j=METAL_BEGIN;j<METAL_END;++j)s<<' '<<c.X[j]/c.Z();
 s<<"\nhydrogen 3 .999 .9995 .99999\nmetallicity 3 1e-8 1e-7 1e-6\nlog_teff 3";
 for(double t:{4600.,4800.,5200.})s<<' '<<std::log10(t);
 s<<"\nlog_g 2 5.4 5.65\ndata\n";int k=0;
 for(double h:{.999,.9995,.99999})for(double z:{1e-8,1e-7,1e-6})for(double t:{4600.,4800.,5200.})for(double g:{5.4,5.65}){
  if(k++==missing)s<<"0\n";else s<<"1 "<<lt(h,z,std::log10(t),g)<<' '<<lp(h,z,std::log10(t),g)<<'\n';
 }
 return s.str();
}
double share_lt(double u,double z,double t,double g){return t+.2*u+50*z-.01*g;}
double share_lp(double u,double z,double t,double g){return 7+.4*u-100*z-2*(t-3.66)+.7*(g-5.5);}
std::string share_fixture(){
 std::ostringstream s;s<<std::setprecision(17);
 s<<"EMBER_TRACE_HELIUM_ATMOSPHERE 2\nsource \"analytic coordinate test\"\napproximation \"test only\"\nbasis baryon_mass\ntau 100\nmaximum_helium3 .0003\nmetal_pattern";
 auto c=comp(.7,0,.02);for(std::size_t j=METAL_BEGIN;j<METAL_END;++j)s<<' '<<c.X[j]/c.Z();
 s<<"\nhydrogen_share 2 .98 .999999\nmetallicity 2 0 .0001\nlog_teff 2 "<<std::log10(4500.)<<' '<<std::log10(4800.)<<"\nlog_g 2 5.3 5.7\ndata\n";
 for(double u:{.98,.999999})for(double z:{0.,.0001})for(double t:{4500.,4800.})for(double g:{5.3,5.7})
  s<<"1 "<<share_lt(u,z,std::log10(t),g)<<' '<<share_lp(u,z,std::log10(t),g)<<'\n';
 return s.str();
}

}
int main(){try{
 Gas eos;const auto approx=TraceHeliumAtmosphereGrid::Approximation::neglect_atmospheric_helium3;
 std::istringstream in(fixture());TraceHeliumAtmosphereGrid grid(eos,in,approx,.0003);
 auto c=comp(.9992,.00004,5e-8);double teff=4700,g=std::pow(10.,5.52);auto s=grid.eval(teff,g,c);
 require(close(std::log10(s.T),lt(c.h1(),c.Z(),std::log10(teff),5.52)),"interior T differs from independent multiaffine expression");
 require(close(std::log10(s.Pgas),lp(c.h1(),c.Z(),std::log10(teff),5.52)),"interior P differs from independent multiaffine expression");
 require(close(eos.eval(s.T,s.rho,c).P,s.P),"pressure inversion does not use actual composition");
 auto zc=comp(c.h1(),0,c.Z());auto zero=grid.eval(teff,g,zc);
 require(zero.T==s.T&&zero.Pgas==s.Pgas,"documented isotope approximation must affect source boundary only");
 require(std::abs(zero.rho/s.rho-1)>1e-9,"EOS density silently projected helium3 away");
 double max_derivative_error=0;
 for(int axis=0;axis<2;++axis){double h=1e-5;double valuesT[4],valuesP[4];int j=0;
  for(int sign:{-2,-1,1,2}){auto q=grid.eval(teff*(axis==0?std::exp(sign*h):1),g*(axis==1?std::exp(sign*h):1),c);valuesT[j]=std::log(q.T);valuesP[j++]=std::log(q.P);}
  for(int field=0;field<2;++field){auto v=field==0?valuesT:valuesP;double fd=(v[0]-8*v[1]+8*v[2]-v[3])/(12*h);double analytic=field==0?(axis==0?s.dlnT_dlnTeff:s.dlnT_dlng):(axis==0?s.dlnP_dlnTeff:s.dlnP_dlng);double error=std::abs(fd-analytic)/(1+std::abs(analytic));max_derivative_error=std::max(max_derivative_error,error);require(error<2e-9,"temperature/gravity derivative fails independent differences");}
 }
 auto responses=grid.composition_response(teff,g,c);
 for(int axis=0;axis<2;++axis){double h=axis==0?1e-7:1e-11;double valuesT[4],valuesP[4];int j=0;
  for(int sign:{-2,-1,1,2}){auto q=grid.eval(teff,g,comp(c.h1()+(axis==0?sign*h:0),c.X[1],c.Z()+(axis==1?sign*h:0)));valuesT[j]=std::log(q.T);valuesP[j++]=std::log(q.P);}
  for(int field=0;field<2;++field){auto v=field==0?valuesT:valuesP;double fd=(v[0]-8*v[1]+8*v[2]-v[3])/(12*h);double analytic=field==0?(axis==0?responses.dlnT_dXH:responses.dlnT_dZ):(axis==0?responses.dlnP_dXH:responses.dlnP_dZ);double error=std::abs(fd-analytic)/(1+std::abs(analytic));max_derivative_error=std::max(max_derivative_error,error);require(error<3e-7,"composition derivative fails independent differences");}
 }
 for(double h:{.999,.9995,.99999})for(double z:{1e-8,1e-7,1e-6})for(double t:{4600.,4800.,5200.})for(double lg:{5.4,5.65}){auto q=grid.eval(t,std::pow(10.,lg),comp(h,0,z));require(close(std::log10(q.T),lt(h,z,std::log10(t),lg))&&close(std::log10(q.Pgas),lp(h,z,std::log10(t),lg)),"source node not recovered");}
 for(auto bad:{comp(.9989,0,5e-8),comp(.9992,.00031,5e-8),comp(.9992,0,9e-9),comp(.9992,0,1.1e-6)})require(rejects([&]{grid.eval(teff,g,bad);}),"out-of-domain composition accepted");
 require(rejects([&]{grid.eval(4599,g,c);})&&rejects([&]{grid.eval(5201,g,c);})&&rejects([&]{grid.eval(teff,std::pow(10.,5.39),c);})&&rejects([&]{grid.eval(teff,std::pow(10.,5.651),c);}),"out-of-domain T or gravity accepted");
 auto bad=c;bad.basis=AbundanceBasis::atomic_mass;require(rejects([&]{grid.eval(teff,g,bad);}),"wrong mass convention accepted");
 bad=c;bad.X[3]+=1e-10;bad.X[4]-=1e-10;require(rejects([&]{grid.eval(teff,g,bad);}),"unsupported metal pattern accepted");
 std::istringstream missing(fixture(0));TraceHeliumAtmosphereGrid holes(eos,missing,approx,.0003);
 require(!holes.covers(teff,g,c)&&rejects([&]{holes.eval(teff,g,c);}),"missing value/derivative corner accepted");
 std::istringstream edge(fixture(36));TraceHeliumAtmosphereGrid edges(eos,edge,approx,.0003);
 require(edges.covers(teff,g,comp(.9995,0,5e-8)),"complete lower stencil at exact knot was lost");
 require(!edges.covers(teff,g,comp(.9998,0,5e-8)),"missing upper stencil silently extrapolated");
 require(rejects([&]{std::istringstream too(fixture(-1,.0004));TraceHeliumAtmosphereGrid q(eos,too,approx,.0003);}),"table exceeded caller-selected isotope bound");
 {
  std::istringstream source(share_fixture());TraceHeliumAtmosphereGrid share(eos,source,approx,.0003);
  for(double u:{.98,.99,.9999,.999999})for(double z:{0.,.000025,.0001}){
   auto mixture=comp(u*(1-z),0,z);auto q=share.eval(4660,std::pow(10.,5.45),mixture);
   require(close(std::log10(q.T),share_lt(u,z,std::log10(4660.),5.45))&&close(std::log10(q.Pgas),share_lp(u,z,std::log10(4660.),5.45)),"hydrogen-share values differ from analytic atmosphere");
  }
  auto mixture=comp(.9999*(1-5e-5),0,5e-5);auto d=share.composition_response(4660,std::pow(10.,5.45),mixture);
  for(int axis=0;axis<2;++axis){
   const double h=axis==0?1e-7:1e-8;
   auto lo=share.eval(4660,std::pow(10.,5.45),comp(mixture.h1()-(axis==0?h:0),0,mixture.Z()-(axis==1?h:0)));
   auto hi=share.eval(4660,std::pow(10.,5.45),comp(mixture.h1()+(axis==0?h:0),0,mixture.Z()+(axis==1?h:0)));
   const double a=axis==0?d.dlnT_dXH:d.dlnT_dZ,b=axis==0?d.dlnP_dXH:d.dlnP_dZ;
   const double error=std::max(std::abs(std::log(hi.T/lo.T)/(2*h)-a)/(1+std::abs(a)),std::abs(std::log(hi.P/lo.P)/(2*h)-b)/(1+std::abs(b)));
   max_derivative_error=std::max(max_derivative_error,error);
   require(error<1e-6,"hydrogen-share composition derivative fails independent differences");
  }
  require(!share.covers(4660,std::pow(10.,5.45),comp(.99995,0,.0001)),"negative helium accepted");
  require(!share.covers(4660,std::pow(10.,5.45),comp(.9798,0,.0001)),"outside hydrogen-share domain accepted");
  const auto support=share.support();
  require(close(support.hydrogen[0],.98*(1-.0001))&&close(support.hydrogen[1],.999999),"support did not return physical hydrogen bounds");
 }
 std::cout<<std::setprecision(17)<<"{\"checks\":"<<checks<<",\"maximum_derivative_error\":"<<max_derivative_error<<",\"actual_isotope_density_fractional_difference\":"<<s.rho/zero.rho-1<<"}\n";
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

#include "ember/atmosphere_metal_interval.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
using namespace ember;
unsigned checks=0;
void require(bool x,const char* m){++checks;if(!x)throw std::runtime_error(m);}
bool close(double a,double b,double tol=2e-11){return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});}
template<class F>bool rejects(F f){try{f();return false;}catch(const std::exception&){return true;}}
Composition comp(double h,double y,double z){auto c=solar_scaled(h,z);c.X[1]=y;c.X[2]-=y;c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;return c;}
class Gas final:public Eos{public:EosState eval(double t,double rho,const Composition& c)const override {EosState s{};double pg=constants::R_gas*c.mu_ions_inv()*rho*t,pr=constants::a_rad*std::pow(t,4)/3;s.P=pg+pr;s.chiT=(pg+4*pr)/s.P;s.chiRho=pg/s.P;return s;}const char* name()const override{return "analytic gas";}};
class Boundary final:public Atmosphere{const Eos& eos;bool low;public:Boundary(const Eos& e,bool l):eos(e),low(l){}AtmosphereState eval(double t,double g,const Composition& c)const override {const double eps=8*std::numeric_limits<double>::epsilon();if(!(t>0&&g>0)||!std::isfinite(t+g)||c.h1()>(low?.99999:.9995)||(low?c.Z()>1e-6*(1+eps):c.Z()<1e-5*(1-eps)))throw std::domain_error("synthetic source domain");AtmosphereState s{};s.T=(low?5000:6000)*std::pow(t/4700.,low?.9:1.1)*std::pow(g/3e5,low?.02:.07);s.Pgas=(low?2e7:3e7)*std::pow(t/4700.,low?.3:.6)*std::pow(g/3e5,low?.8:1.2);double pr=constants::a_rad*std::pow(s.T,4)/3;s.P=s.Pgas+pr;s.tau=100;s.dlnT_dlnTeff=low?.9:1.1;s.dlnT_dlng=low?.02:.07;s.dlnP_dlnTeff=(s.Pgas*(low?.3:.6)+4*pr*s.dlnT_dlnTeff)/s.P;s.dlnP_dlng=(s.Pgas*(low?.8:1.2)+4*pr*s.dlnT_dlng)/s.P;s.rho=eos.rho_from_PT(s.T,s.P,c);return s;}const char* name()const override{return "synthetic analytic atmosphere";}};
int main(){try{Gas eos;Boundary low(eos,true),high(eos,false);MetalIntervalAtmosphere bridge(eos,low,high,1e-6,1e-5);double maxerr=0;
for(double z:{1e-8,1e-6,3e-6,8e-6,1e-5,1e-4}){auto c=comp(.9992,.00004,z);auto s=bridge.eval(4700,3e5,c);double f=std::clamp((z-1e-6)/9e-6,0.,1.);require(close(std::log(s.T),(1-f)*std::log(5000)+f*std::log(6000)),"independent temperature interpolation");require(close(std::log(s.Pgas),(1-f)*std::log(2e7)+f*std::log(3e7)),"independent gas pressure interpolation");require(close(eos.eval(s.T,s.rho,c).P,s.P),"EOS inversion at actual isotope composition");
for(int axis=0;axis<2;++axis){double h=1e-5,tv[4],pv[4];int n=0;for(int sign:{-2,-1,1,2}){auto q=bridge.eval(4700*std::exp(axis==0?sign*h:0),3e5*std::exp(axis==1?sign*h:0),c);tv[n]=std::log(q.T);pv[n++]=std::log(q.P);}for(int field=0;field<2;++field){auto v=field?pv:tv;double fd=(v[0]-8*v[1]+8*v[2]-v[3])/(12*h);double deriv=field?(axis?s.dlnP_dlng:s.dlnP_dlnTeff):(axis?s.dlnT_dlng:s.dlnT_dlnTeff);double err=std::abs(fd-deriv)/(1+std::abs(deriv));maxerr=std::max(maxerr,err);require(err<2e-9,"physical T/g derivative");}}}
for(double factor:{1.,1+std::numeric_limits<double>::epsilon()}){auto c=comp(.99999,0,1e-6*factor);auto a=bridge.eval(4700,3e5,c),b=low.eval(4700,3e5,c);require(a.T==b.T&&a.P==b.P&&a.rho==b.rho,"lower endpoint incorrectly queries unavailable upper H");}
require(rejects([&]{bridge.eval(4700,3e5,comp(.9998,0,5e-6));}),"unavailable upper source extrapolated");require(rejects([&]{bridge.eval(4700,3e5,comp(.99999,.000009,5e-6));}),"negative helium anchor accepted");
auto c=comp(.9992,.00004,5e-6);auto zero=comp(.9992,0,5e-6);auto a=bridge.eval(4700,3e5,c),b=bridge.eval(4700,3e5,zero);require(a.T==b.T&&a.P==b.P&&std::abs(a.rho/b.rho-1)>1e-9,"trace isotope discarded from EOS");c.basis=AbundanceBasis::atomic_mass;require(rejects([&]{bridge.eval(4700,3e5,c);}),"wrong composition basis");
std::cout<<std::setprecision(17)<<"{\"checks\":"<<checks<<",\"maximum_derivative_error\":"<<maxerr<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

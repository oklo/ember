// Version-3 (irregular Delaunay) TraceHeliumAtmosphereGrid tests. Fixture: Kuhn triangulation of a 4D box in raw
// coordinates (h=X/(1-Z), Z, log10 Teff, log10 g): 16 vertices, 24 simplices sharing the main diagonal.
#include "ember/atmosphere_trace_helium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
using namespace ember;
namespace {
int checks=0;
void require(bool ok,const char* m){++checks;if(!ok)throw std::runtime_error(m);}
bool close(double a,double b,double tol){return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});}
template<class F>bool rejects(F f){try{f();return false;}catch(const std::exception&){return true;}}
class Gas final:public Eos{public:
 EosState eval(double t,double rho,const Composition& c)const override{EosState s{};double pg=constants::R_gas*c.mu_ions_inv()*rho*t,pr=constants::a_rad*std::pow(t,4)/3;s.P=pg+pr;s.chiT=(pg+4*pr)/s.P;s.chiRho=pg/s.P;return s;}
 const char* name()const override{return "analytic test gas";}};
Composition comp(double h,double y,double z){auto c=solar_scaled(h,z);c.X[1]=y;c.X[2]-=y;c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;return c;}
const std::array<double,4> lo{.98,0.,std::log10(4500.),5.3},hi{.999999,1e-4,std::log10(4800.),5.7};
const std::array<double,4> A{.2,50.,1.,-.01},B{.4,-100.,-2.,.7};
double affine(const std::array<double,4>& a,double c,const std::array<double,4>& q){return c+a[0]*q[0]+a[1]*q[1]+a[2]*q[2]+a[3]*q[3];}
std::array<double,4> vertex(unsigned bits){std::array<double,4> q{};for(int k=0;k<4;++k)q[k]=bits&(1u<<k)?hi[k]:lo[k];return q;}
// field: 0 affine; 1 affine plus a vertex-only perturbation (piecewise-linear, not affine)
std::string fixture(int field,std::string mutate=""){
 std::ostringstream s;s<<std::setprecision(17);auto c=comp(.7,0,.02);
 s<<"EMBER_TRACE_HELIUM_ATMOSPHERE 3\nsource \"irregular analytic test\"\napproximation \"test only\"\nbasis baryon_mass\ntau 100\nmaximum_helium3 .0003\nmetal_pattern";
 for(std::size_t j=METAL_BEGIN;j<METAL_END;++j)s<<' '<<c.X[j]/c.Z();
 if(mutate=="offset")s<<"\nirregular_scale 1e999 0 3.66 5.5 .01 1e-5 .01 .1\n";
 else s<<(mutate=="scale"?"\nirregular_scale .99 0 3.66 5.5 .01 0 .01 .1\n":"\nirregular_scale .99 0 3.66 5.5 .01 1e-5 .01 .1\n");
 s<<"vertices 16\n";
 for(unsigned b=0;b<16;++b){auto q=vertex(b);double pT=field?.01*std::sin(3.*b):0,pP=field?.02*std::cos(5.*b):0;
  s<<q[0]<<' '<<q[1]<<' '<<q[2]<<' '<<q[3]<<' '<<(mutate=="overflow"?400.:affine(A,3.6,q)+pT)<<' '<<affine(B,7,q)+pP<<'\n';}
 std::array<int,4> perm{0,1,2,3};std::ostringstream simp;int m=0;
 do{unsigned b=0;simp<<b;for(int k:perm){b|=1u<<k;simp<<' '<<b;}simp<<'\n';++m;}while(std::next_permutation(perm.begin(),perm.end()));
 if(mutate=="degenerate"){simp<<"0 1 3 7 7\n";++m;}
 if(mutate=="index"){simp<<"0 1 3 7 16\n";++m;}
 if(mutate=="single")s<<"simplices 1\n0 1 3 7 15\n";
 else s<<"simplices "<<m<<'\n'<<simp.str();
 if(mutate=="trailing")s<<"extra\n";
 return s.str();
}
}
int main(){try{
 Gas eos;const auto approx=TraceHeliumAtmosphereGrid::Approximation::neglect_atmospheric_helium3;
 std::istringstream in(fixture(0));TraceHeliumAtmosphereGrid grid(eos,in,approx,.0003);
 std::mt19937 rng(7);std::uniform_real_distribution<double> U(.02,.98);double worst=0;
 // 1. affine exactness and analytic derivatives at random interior points
 for(int n=0;n<200;++n){std::array<double,4> q{};for(int k=0;k<4;++k)q[k]=lo[k]+U(rng)*(hi[k]-lo[k]);
  auto c=comp(q[0]*(1-q[1]),0,q[1]);auto s=grid.eval(std::pow(10.,q[2]),std::pow(10.,q[3]),c);
  require(close(std::log10(s.T),affine(A,3.6,q),1e-13)&&close(std::log10(s.Pgas),affine(B,7,q),1e-13),"affine field not reproduced");
  require(close(s.dlnT_dlnTeff,A[2],1e-11)&&close(s.dlnT_dlng,A[3],1e-11),"affine Teff/g derivative");
  auto r=grid.composition_response(std::pow(10.,q[2]),std::pow(10.,q[3]),c);
  const double dx=std::log(10.)*A[0]/(1-q[1]),dz=std::log(10.)*(A[1]+A[0]*q[0]/(1-q[1]));
  worst=std::max({worst,std::abs(r.dlnT_dXH-dx)/std::abs(dx),std::abs(r.dlnT_dZ-dz)/std::abs(dz)});
  require(worst<1e-9,"affine composition derivative");}
 // 2. non-affine piecewise field: face continuity, one-sided derivatives, finite-offset validation
 std::istringstream in2(fixture(1));TraceHeliumAtmosphereGrid pw(eos,in2,approx,.0003);
 double jump=0,fdworst=0;
 for(int n=0;n<50;++n){std::array<double,4> u{};for(auto& x:u)x=U(rng);u[1]=u[0];          // on the face u0=u1
  auto at=[&](double e0){std::array<double,4> q{};for(int k=0;k<4;++k)q[k]=lo[k]+(u[k]+(k==0?e0:0))*(hi[k]-lo[k]);return q;};
  auto ev=[&](const std::array<double,4>& q){return pw.eval(std::pow(10.,q[2]),std::pow(10.,q[3]),comp(q[0]*(1-q[1]),0,q[1]));};
  const double e=1e-9;auto a=ev(at(-e)),b=ev(at(e)),f=ev(at(0));
  jump=std::max({jump,std::abs(std::log(a.T/b.T)),std::abs(std::log(f.T/a.T))});
  // one-sided derivative on each side against a finite offset kept on that side (piecewise linear: exact)
  for(double side:{-1.,1.}){auto q0=at(side*1e-4);auto s0=ev(q0);const double h=2e-6;
   auto qp=q0,qm=q0;qp[2]+=h;qm[2]-=h;auto sp=ev(qp),sm=ev(qm);
   const double fd=std::log(sp.T/sm.T)/(2*h*std::log(10.));fdworst=std::max(fdworst,std::abs(fd-s0.dlnT_dlnTeff));}
 }
 require(jump<1e-7,"values discontinuous across a shared face");
 require(fdworst<1e-7,"one-sided derivative fails finite-offset check");
 // 3. outside hull and composition validation
 auto c=comp(.99*(1-5e-5),0,5e-5);
 require(grid.covers(4650,std::pow(10.,5.5),c),"interior not covered");
 require(!grid.covers(4400,std::pow(10.,5.5),c)&&rejects([&]{grid.eval(4400,std::pow(10.,5.5),c);}),"outside hull accepted");
 require(!grid.covers(4650,std::pow(10.,5.8),c),"outside gravity accepted");
 require(!grid.covers(4650,std::pow(10.,5.5),comp(.99*(1-5e-5),.00031,5e-5)),"He3 bound ignored");
 auto bad=c;bad.X[3]+=1e-10;bad.X[4]-=1e-10;require(!grid.covers(4650,std::pow(10.,5.5),bad),"metal pattern ignored");
 // 4/5. malformed and degenerate files
 for(const char* m:{"degenerate","index","trailing","scale","offset","overflow"}){std::istringstream s(fixture(0,m));
  require(rejects([&]{TraceHeliumAtmosphereGrid q(eos,s,approx,.0003);}),m);}
 // A bounding box does not imply coverage between every set of source vertices.
 std::istringstream single(fixture(0,"single"));TraceHeliumAtmosphereGrid hull(eos,single,approx,.0003);
 auto hv=[&](std::array<double,4> u){std::array<double,4> q{};for(int k=0;k<4;++k)q[k]=lo[k]+u[k]*(hi[k]-lo[k]);return q;};
 auto inside=hv({.8,.6,.4,.2}),outside=hv({.2,.4,.6,.8});
 require(hull.has_missing_states(),"irregular hull advertised as complete rectangle");
 require(hull.covers(std::pow(10.,inside[2]),std::pow(10.,inside[3]),comp(inside[0]*(1-inside[1]),0,inside[1])),"simplex interior rejected");
 require(!hull.covers(std::pow(10.,outside[2]),std::pow(10.,outside[3]),comp(outside[0]*(1-outside[1]),0,outside[1])),"bounding-box hole accepted");
 std::cout<<std::setprecision(17)<<"{\"checks\":"<<checks<<",\"affine_composition_derivative_rel\":"<<worst<<",\"face_jump\":"<<jump<<",\"one_sided_fd_error\":"<<fdworst<<"}\n";
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

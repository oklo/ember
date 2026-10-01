#include "ember/cold_helium_base.hpp"
#include "ember/detail/ion_ocp_components.hpp"
#include "ember/constants.hpp"
#include "ember/gs98_mixture.hpp"
#include "ember/detail/ion_mixture.hpp"
#include "ember/ion_quantum.hpp"
#include "ember/ion_phase.hpp"
#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace ember {
namespace {
using K=detail::Taylor3<2>;
namespace io=detail::ioffe;
constexpr double pi=std::numbers::pi;
constexpr double me=9.1093837015e-28,clight=2.99792458e10,hbar=1.054571817e-27;

[[noreturn]] void fail(const std::string& m){throw std::domain_error("cold He base: "+m);}

// Every domain limit of the base, in scalar arithmetic (cheap; shared by evaluation and validation).
void check_domain(double T,double rho,const Composition& c,std::size_t channels,const ColdHeliumOptions& o) {
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0)||(channels!=1&&channels!=4&&channels!=10))
    fail("invalid state or derivative request");
  if(!(o.liquid_continuation_gamma==0. || (std::isfinite(o.liquid_continuation_gamma)
      && o.liquid_continuation_gamma>=175. && o.liquid_continuation_gamma<=300.)))
    fail("metal-liquid continuation must be zero or between Gamma 175 and 300");
  if(c.basis!=AbundanceBasis::baryon_mass||c.metal_inventory!=MetalInventory::gs98||c[Species::H2]!=0)
    fail("requires baryonic GS98 material after D mapping");
  for(double fraction:c.X)
    if(!std::isfinite(fraction)||fraction<0)fail("invalid mass fraction");
  if(std::abs(c.sum()-1)>1e-8)fail("mass fractions do not sum to one");
  if(!(c.X[0]<=o.max_mixture_hydrogen&&c.Z()<=o.max_mixture_metals))
    fail("liquid mixture limited to He-dominated material (X and Z limits)");
  if(!(rho>=o.minimum_density))fail("density below the assessed He pressure-ionization floor");
  const double X=c.X[0],Y3=c.X[1],Zm=c.Z(),X4=1-X-Y3-Zm;
  const double logne=std::log(rho)+std::log(X+(2./3)*Y3+.5*X4+gs98_ion_moment(1)*Zm)-std::log(constants::amu);
  const double x=(hbar/(me*clight))*std::cbrt(3*pi*pi)*std::exp(logne/3);
  const double eF=me*clight*clight*(std::sqrt(1+x*x)-1);
  if(!(constants::kB*T/eF<=o.maximum_electron_temperature_ratio))fail("electrons not degenerate enough for the Sommerfeld form");
  for(const double A:{3.,4.}){const auto p=io::plasma(std::log(T),logne,A,2.);
    if(!(p.rsi>=500&&p.rsi<=1.2e5&&p.tpt<=30))fail("He isotope outside the Baiko-Chugunov fitted quantum domain");}
  if((X>0||channels>1)
      &&!(io::plasma(std::log(T),logne,1.,1.).tpt<=(X<=o.trace_hydrogen?o.maximum_trace_hydrogen_quantum:o.maximum_hydrogen_quantum)))
    fail("hydrogen quantum parameter beyond the production limit");
  {   // helium host coupling; the trace-metal mean coupling is not a phase criterion (docs/DENSE_EOS.md)
    const double rs=std::exp((std::log(3/(4*pi))-logne)/3)/io::bohr;
    const double gamma=io::hartree_k/(rs*T)*std::pow(2.,5./3);
    if(o.mixture_phase) {
      if(!(gamma<=o.max_phase_helium_gamma))fail("helium coupling beyond the fitted liquid and solid branches");
    } else if(!(gamma<=o.max_helium_gamma))fail("helium coupling beyond the assessed liquid range; phase unsupported");
  }
}
std::array<K,2> electron_parts(const K& lt,const K& ln) {

    const K n=exp(ln),x=(hbar/(me*clight))*std::cbrt(3*pi*pi)*exp(ln/3),root=sqrt(1+x*x);
    const double Kc=std::pow(me,4)*std::pow(clight,5)/(8*pi*pi*std::pow(hbar,3));
    const K e_kin=Kc*(x*(2*x*x+1)*root-log(x+root))-n*(me*clight*clight);
    const K g=(x*me*clight)*(me*clight*clight)*root/(pi*pi*std::pow(hbar,3)*clight*clight);
    const K kT=constants::kB*exp(lt);
    const K rs=exp((std::log(3/(4*pi))-ln)/3)/io::bohr;
    return {constants::R_gas*(e_kin-(pi*pi/6)*g*kT*kT)/(kT*n),
            constants::R_gas*io::excor7(rs,io::hartree_k/rs*exp(-lt))};
}
void check_join(const ColdHeliumOptions& o) {
  if(!(o.join_cold>0&&o.join_hot>o.join_cold*1.05&&std::isfinite(o.join_hot)))fail("invalid join window");
  if(!(o.hydrogen_join_full>=0&&o.hydrogen_join_zero>o.hydrogen_join_full&&o.hydrogen_join_zero<=o.max_mixture_hydrogen))
    fail("invalid hydrogen join interval");
}
// Quintic C2 smoothstep, 1 at lo and 0 at hi, with its first three derivatives in x.
std::array<double,4> falling(double x,double lo,double hi) {
  if(x<=lo)return {1,0,0,0};
  if(x>=hi)return {0,0,0,0};
  const double h=hi-lo,u=(x-lo)/h;
  return {1-u*u*u*(10+u*(-15+6*u)),-30*u*u*(1-u)*(1-u)/h,-60*u*(1-u)*(1-2*u)/(h*h),-60*(1-6*u+6*u*u)/(h*h*h)};
}
ColdHeliumOptions transition_options(const ColdHeliumOptions& source) {
  auto o=source;
  o.dense_transition=false;o.mixture_phase=false;
  o.max_mixture_hydrogen=1.;o.hydrogen_join_full=.99;o.hydrogen_join_zero=1.;
  o.minimum_density=300.;o.join_cold=2e5;o.join_hot=3e5;
  o.maximum_hydrogen_quantum=2.5;o.maximum_trace_hydrogen_quantum=2.5;
  return o;
}
} // namespace

double dense_hhe_transition_weight(double T,double rho) {
  if(!(T>0&&rho>0&&std::isfinite(T)&&std::isfinite(rho)))fail("invalid H/He transition state");
  return falling(std::log(T),std::log(2e5),std::log(3e5))[0]
      *(1-falling(std::log(rho),std::log(300.),std::log(600.))[0]);
}
void dense_hhe_transition_validate(double T,double rho,const Composition& c,const ColdHeliumOptions& source) {
  if(!(T>=1.2e5&&rho<=4e3&&c.Z()<=1e-8))fail("outside assessed dense H/He transition range");
  const auto o=transition_options(source);
  for(double t:{T,o.join_cold,o.join_hot})cold_helium_validate(t,rho,c,o);
}
std::array<HelmholtzJet,10> dense_hhe_transition_jets(const ColdHeliumTable& table,double T,double rho,
    const Composition& c,std::size_t channels,const ColdHeliumOptions& source) {
  const double w=dense_hhe_transition_weight(T,rho);
  if(w==0)return table(T,rho,c,channels);
  dense_hhe_transition_validate(T,rho,c,source);
  const auto o=transition_options(source);
  auto base=cold_helium_material_jets(T,rho,c,channels,o);
  const auto alignment=cold_helium_alignment_jets(table,T,rho,c,channels,o);
  for(std::size_t k=0;k<channels;++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)
    base[k][i][j]+=alignment[k][i][j];
  if(w==1)return base;
  const auto old=table(T,rho,c,channels);
  const auto wt=falling(std::log(T),std::log(o.join_cold),std::log(o.join_hot));
  auto wr=falling(std::log(rho),std::log(300.),std::log(600.));
  wr[0]=1-wr[0];for(unsigned k=1;k<4;++k)wr[k]=-wr[k];
  constexpr unsigned binomial[4][4]={{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}};
  auto result=old;
  for(std::size_t k=0;k<channels;++k) {
    const unsigned order=k==0?0:(k<4?1:2);
    for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j+order<=3;++j)
      for(unsigned a=0;a<=i;++a)for(unsigned b=0;b<=j;++b)
        result[k][i][j]+=binomial[i][a]*binomial[j][b]*wt[a]*wr[b]
            *(base[k][i-a][j-b]-old[k][i-a][j-b]);
  }
  return result;
}

std::array<HelmholtzJet,10> cold_helium_material_jets(double T,double rho,const Composition& c,std::size_t channels,
    const ColdHeliumOptions& o) {
  check_domain(T,rho,c,channels,o);
  const double Ye=c.X[0]+(2./3)*c.X[1]+.5*c.X[2]+gs98_ion_moment(1)*c.Z();
  const K lt=K::variable(std::log(T),0),ln=K::variable(std::log(rho*Ye/constants::amu),1);
  const auto electrons=electron_parts(lt,ln);
  const K E=electrons[0]+electrons[1];
  auto out=detail::common_density_ion_jets(T,rho,c,channels,[&](double,double,double A,double Z) {
    const K f=io::classical_liquid(io::plasma(lt,ln,A,Z),o.liquid_continuation_gamma);
    const K total=(constants::R_gas/A)*f+(Z/A)*E;
    // Quantum term: the production Baiko-Chugunov implementation itself (He isotopes; H and metals by default).
    const HelmholtzJet q=detail::bc22_liquid_quantum_per_mass(lt.value(),ln.value(),A,Z);
    HelmholtzJet out{};
    for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)out[i][j]=total.derivative({i,j})+q[i][j];
    return out;
  });
  if(o.mixture_phase) {
    MixturePhaseOptions p;p.liquid_continuation_gamma=o.liquid_continuation_gamma;p.width=o.phase_width;p.minimum_solid_gamma=o.minimum_solid_gamma;p.trace_hydrogen=o.trace_hydrogen;
    const auto d=ion_mixture_phase_difference_jets(T,rho,c,channels,p);
    for(std::size_t k=0;k<channels;++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)out[k][i][j]+=d[k][i][j];
  }
  return out;
}
HelmholtzJet cold_helium_component_jet(int component,double T,double rho,const Composition& c,const ColdHeliumOptions& o) {
  check_domain(T,rho,c,1,o);
  if(component<0||component>1)fail("electron component must be 0 (ideal) or 1 (exchange-correlation)");
  const double Ye=c.X[0]+(2./3)*c.X[1]+.5*c.X[2]+gs98_ion_moment(1)*c.Z();
  const auto parts=electron_parts(K::variable(std::log(T),0),K::variable(std::log(rho*Ye/constants::amu),1));
  HelmholtzJet out{};
  for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)out[i][j]=Ye*parts[component].derivative({i,j});
  return out;
}
void cold_helium_validate(double T,double rho,const Composition& c,const ColdHeliumOptions& o) {
  check_join(o);check_domain(T,rho,c,10,o);
}
double cold_helium_join_weight(double T,const Composition& c,const ColdHeliumOptions& o) {
  check_join(o);
  return falling(std::log(T),std::log(o.join_cold),std::log(o.join_hot))[0]*falling(c.X[0],o.hydrogen_join_full,o.hydrogen_join_zero)[0];
}
std::array<double,6> cold_helium_solid_response(double T,double rho,const Composition& c,const ColdHeliumOptions& o) {
  check_join(o);
  if(!o.mixture_phase)return {};
  const auto wt=falling(std::log(T),std::log(o.join_cold),std::log(o.join_hot));
  const auto wx=falling(c.X[0],o.hydrogen_join_full,o.hydrogen_join_zero);
  if(wt[0]*wx[0]==0)return {};
  check_domain(T,rho,c,10,o);
  MixturePhaseOptions p;p.liquid_continuation_gamma=o.liquid_continuation_gamma;
  p.width=o.phase_width;p.minimum_solid_gamma=o.minimum_solid_gamma;p.trace_hydrogen=o.trace_hydrogen;
  const auto raw=ion_mixture_phase_weight_response(T,rho,c,p);
  auto result=raw;
  for(auto& v:result)v*=wt[0]*wx[0];
  result[1]+=raw[0]*wt[1]*wx[0];
  result[3]+=raw[0]*wt[0]*wx[1];
  return result;
}

std::array<HelmholtzJet,10> cold_helium_alignment_jets(const ColdHeliumTable& table,double T,double rho,
    const Composition& c,std::size_t channels,const ColdHeliumOptions& o) {
  check_join(o);
  const double T1=o.join_cold,T2=o.join_hot;
  std::array<HelmholtzJet,10> t1,t2;
  try{t1=table(T1,rho,c,channels);t2=table(T2,rho,c,channels);}
  catch(const std::exception& e){fail(std::string("join anchor outside table support: ")+e.what());}
  const auto b1=cold_helium_material_jets(T1,rho,c,channels,o),b2=cold_helium_material_jets(T2,rho,c,channels,o);
  std::array<HelmholtzJet,10> g{};
  for(std::size_t k=0;k<channels;++k){
    const unsigned order=k==0?0:(k<4?1:2);
    for(unsigned j=0;j+order<=3;++j){
      const double d1=t1[k][0][j]-b1[k][0][j],d2=t2[k][0][j]-b2[k][0][j];
      const double b=(d1-d2)/(1/T1-1/T2),a=d1-b/T1;
      for(unsigned i=0;i+j+order<=3;++i)g[k][i][j]=i==0?a+b/T:b*((i%2)?-1.:1.)/T;
    }
  }
  return g;
}

std::array<HelmholtzJet,10> cold_helium_joined_jets(const ColdHeliumTable& table,double T,double rho,
    const Composition& c,std::size_t channels,const ColdHeliumOptions& o) {
  check_join(o);
  if(!(std::isfinite(T)&&T>0))fail("invalid state or derivative request");
  const auto wt=falling(std::log(T),std::log(o.join_cold),std::log(o.join_hot));
  const auto wx=falling(c.X[0],o.hydrogen_join_full,o.hydrogen_join_zero);
  if(wt[0]==0||wx[0]==0)return table(T,rho,c,channels);
  auto base=cold_helium_material_jets(T,rho,c,channels,o);
  const auto g=cold_helium_alignment_jets(table,T,rho,c,channels,o);
  for(std::size_t k=0;k<channels;++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)base[k][i][j]+=g[k][i][j];
  if(wt[0]==1&&wx[0]==1)return base;
  std::array<HelmholtzJet,10> old;
  try{old=table(T,rho,c,channels);}
  catch(const std::exception& e){fail(std::string("fractional join weight needs the production table: ")+e.what());}
  // D = base + G - table; E = w_T D (Leibniz in ln T); F = table + w_X E (Leibniz in X, channels 1 and 4-6).
  constexpr unsigned binomial[4][4]={{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}};
  auto order=[](std::size_t k){return k==0?0u:(k<4?1u:2u);};
  std::array<HelmholtzJet,10> e{};
  for(std::size_t k=0;k<channels;++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j+order(k)<=3;++j)
    for(unsigned a=0;a<=i;++a)e[k][i][j]+=binomial[i][a]*wt[a]*(base[k][i-a][j]-old[k][i-a][j]);
  // Hessian channel of the (X, b) pair: 4 (X,X), 5 (X,Y3), 6 (X,Z).
  std::array<HelmholtzJet,10> out=old;
  for(std::size_t k=0;k<channels;++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j+order(k)<=3;++j) {
    double v=wx[0]*e[k][i][j];
    if(k==1)v+=wx[1]*e[0][i][j];
    else if(k==4)v+=2*wx[1]*e[1][i][j]+wx[2]*e[0][i][j];
    else if(k==5||k==6)v+=wx[1]*e[k-3][i][j];
    out[k][i][j]+=v;
  }
  return out;
}
} // namespace ember

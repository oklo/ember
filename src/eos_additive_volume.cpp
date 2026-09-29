#include "ember/eos_additive_volume.hpp"
#include "ember/constants.hpp"
#include "ember/detail/taylor3.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace ember {
namespace {
double val(double x){return x;}
template<std::size_t N>double val(const detail::Taylor3<N>& x){return x.value();}
using std::exp;using std::log;using detail::exp;using detail::log;
std::size_t cell(const std::vector<double>& nodes,double x){
  if(!std::isfinite(x)||x<nodes.front()||x>nodes.back())throw std::domain_error("additive-volume potential outside source support");
  return std::min<std::size_t>(std::upper_bound(nodes.begin(),nodes.end(),x)-nodes.begin()-1,nodes.size()-2);
}
}
AdditiveVolumePotential::AdditiveVolumePotential(const std::filesystem::path& path){
  std::ifstream in(path);std::string label;int version{};in>>label>>version;
  if(label!="EMBER_ADDITIVE_VOLUME_POTENTIAL"||(version!=1&&version!=2))throw std::runtime_error("invalid additive-volume potential header");
  for(auto entry:{std::pair{&h_,"H"},std::pair{&he_,"He"}}){
    auto& p=*entry.first;std::size_t nt{},nr{};in>>label>>nt>>nr;
    if(label!=entry.second||nt<2||nr<2||nt>1000||nr>1000)throw std::runtime_error("invalid pure-potential dimensions");
    p.t.resize(nt);p.r.resize(nr);p.cells.resize((nt-1)*(nr-1));
    for(auto* axis:{&p.t,&p.r})for(std::size_t i=0;i<axis->size();++i){
      in>>(*axis)[i];if(!in||!std::isfinite((*axis)[i])||(i&&(*axis)[i]<=(*axis)[i-1]))throw std::runtime_error("invalid pure-potential axis");
    }
    p.density_bounds.assign(nt-1,{p.r.front(),p.r.back()});
    if(version==2){
      in>>label;if(label!="density_bounds")throw std::runtime_error("missing pure-potential density bounds");
      for(auto& b:p.density_bounds){
        in>>b[0]>>b[1];
        if(!in||!std::isfinite(b[0])||!std::isfinite(b[1])||b[0]<p.r.front()||b[1]>p.r.back()||b[0]>=b[1])
          throw std::runtime_error("invalid pure-potential density bounds");
      }
    }
    for(auto& c:p.cells)for(auto& x:c){in>>x;if(!in||!std::isfinite(x))throw std::runtime_error("invalid pure-potential coefficient");}
  }
  if(in>>label)throw std::runtime_error("extra data in pure potential");
}
std::array<double,2> AdditiveVolumePotential::bounds(const Pure& p,double t){
  return p.density_bounds[cell(p.t,t)];
}
template<class J>J AdditiveVolumePotential::phi(const Pure& p,const J& t,const J& r,unsigned dr){
  const auto b=bounds(p,val(t));
  if(val(r)<b[0]||val(r)>b[1])throw std::domain_error("pure potential outside supported temperature-density region");
  const auto it=cell(p.t,val(t)),ir=cell(p.r,val(r));const auto& c=p.cells[it*(p.r.size()-1)+ir];
  const J dt=t-p.t[it],d=r-p.r[ir];J f(0);
  for(int i=3;i>=0;--i){
    J a(0);for(int j=3;j>=int(dr);--j){double fac=1;for(unsigned k=0;k<dr;++k)fac*=j-k;a=a*d+c[4*i+j]*fac;}
    f=f*dt+a;
  }
  return f;
}
double AdditiveVolumePotential::density(const Pure& p,double t,double pressure){
  auto state=[&](double r){const auto pr=phi(p,t,r,1),prr=phi(p,t,r,2);
    if(!(pr>0&&pr+prr>0))throw std::domain_error("nonpositive pure-potential compressibility");
    return std::array{t+r+std::log(pr),1+prr/pr};};
  const auto b=bounds(p,t);double lo=b[0],hi=b[1];const auto low=state(lo),high=state(hi);
  if(pressure<low[0]-2e-13||pressure>high[0]+2e-13)throw std::domain_error("pure pressure outside source support");
  if(pressure<=low[0])return lo;if(pressure>=high[0])return hi;
  double r=lo+(hi-lo)*(pressure-low[0])/(high[0]-low[0]);
  for(int i=0;i<60;++i){auto s=state(r);double f=s[0]-pressure;if(std::abs(f)<3e-14)return r;
    if(f<0)lo=r;else hi=r;double next=r-f/s[1];r=next>lo&&next<hi?next:.5*(lo+hi);
    if(hi-lo<5e-14)return r;
  }
  throw std::runtime_error("pure density inversion did not converge");
}
template<class J>J AdditiveVolumePotential::density_jet(const Pure& p,const J& t,const J& pressure){
  J r(density(p,val(t),val(pressure)));
  // The scalar root is already converged. Newton lifts its local Taylor
  // coefficients; three passes suffice through degree three.
  for(int i=0;i<3;++i){const auto pr=phi(p,t,r,1);r=r-(t+r+log(pr)-pressure)/(1+phi(p,t,r,2)/pr);}
  return r;
}
template<class J>J AdditiveVolumePotential::residual(const J& t,const J& r,const J& X,const J& Y3)const{
  const J helium=1-X+Y3/3,n=X+helium/4;
  const double td=val(t),rd=val(r),x=val(X),y=val(helium);
  if(!std::isfinite(rd+x+y)||x<0||val(Y3)<0||x+val(Y3)>1)throw std::domain_error("invalid additive-volume composition");
  auto pressure=[&](const Pure& a,double rr){return td+rr+std::log(phi(a,td,rr,1));};
  const auto bh=bounds(h_,td),bhe=bounds(he_,td);
  double lo=std::max(pressure(h_,bh[0]),pressure(he_,bhe[0]))+1e-12;
  double hi=std::min(pressure(h_,bh[1]),pressure(he_,bhe[1]))-1e-12;
  auto state=[&](double p){double rh=density(h_,td,p),rhe=density(he_,td,p);
    double vh=x*std::exp(-rh),vhe=y*std::exp(-rhe),v=vh+vhe;
    double d=vh/(1+phi(h_,td,rh,2)/phi(h_,td,rh,1))+vhe/(1+phi(he_,td,rhe,2)/phi(he_,td,rhe,1));
    return std::array{-std::log(v)-rd,d/v};};
  if(!(lo<hi)||state(lo)[0]>0||state(hi)[0]<0)throw std::domain_error("mixture density outside common pure support");
  double p=.5*(lo+hi);bool converged=false;
  for(int i=0;i<60;++i){const auto s=state(p);if(std::abs(s[0])<4e-14){converged=true;break;}
    if(s[0]<0)lo=p;else hi=p;const double next=p-s[0]/s[1];p=next>lo&&next<hi?next:.5*(lo+hi);
  }
  if(!converged)throw std::runtime_error("mixture pressure inversion did not converge");
  J lp(p),rh,rhe;
  for(int i=0;i<3;++i){
    rh=density_jet(h_,t,lp);rhe=density_jet(he_,t,lp);const auto vh=X*exp(-rh),vhe=helium*exp(-rhe),v=vh+vhe;
    const auto d=vh/(1+phi(h_,t,rh,2)/phi(h_,t,rh,1))+vhe/(1+phi(he_,t,rhe,2)/phi(he_,t,rhe,1));
    lp=lp-(-log(v)-r)/(d/v);
  }
  rh=density_jet(h_,t,lp);rhe=density_jet(he_,t,lp);
  return X*phi(h_,t,rh)+helium*phi(he_,t,rhe)-constants::R_gas*n*log(n);
}
std::array<HelmholtzJet,6> AdditiveVolumePotential::residual_jets(double T,double rho,double X,double Y3,bool composition)const{
  if(!(std::isfinite(T)&&T>0&&std::isfinite(rho)&&rho>0))throw std::domain_error("invalid additive-volume thermal state");
  std::array<HelmholtzJet,6> result{};
  auto evaluate=[&]<std::size_t N>(){
    using J=detail::Taylor3<N>;const auto t=J::variable(std::log(T),0),r=J::variable(std::log(rho),1);
    J x(X),y(Y3);if constexpr(N>2){x=J::variable(X,2);y=J::variable(Y3,3);}
    const auto f=residual(t,r,x,y);const std::array<std::array<unsigned,2>,6> powers{{{0,0},{1,0},{0,1},{2,0},{1,1},{0,2}}};
    for(std::size_t k=0;k<(N>2?6:1);++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j+powers[k][0]+powers[k][1]<=3;++j){
      typename J::Powers p{};p[0]=i;p[1]=j;if constexpr(N>2){p[2]=powers[k][0];p[3]=powers[k][1];}result[k][i][j]=f.derivative(p);
    }
  };
  if(composition)evaluate.template operator()<4>();else evaluate.template operator()<2>();
  return result;
}
Eos::DensityRange AdditiveVolumePotential::density_range(double T,double X,double Y3)const{
  if(!(T>0&&std::isfinite(T))||X<0||Y3<0||X+Y3>1)throw std::domain_error("invalid additive-volume density-range query");
  const double t=std::log(T),helium=1-X+Y3/3;
  auto pressure=[&](const Pure& a,double r){return t+r+std::log(phi(a,t,r,1));};
  const auto bh=bounds(h_,t),bhe=bounds(he_,t);
  const double lo=std::max(pressure(h_,bh[0]),pressure(he_,bhe[0]))+1e-12;
  const double hi=std::min(pressure(h_,bh[1]),pressure(he_,bhe[1]))-1e-12;
  if(!(lo<hi))throw std::domain_error("no common pure-potential pressure interval");
  auto rho=[&](double p){return 1/(X*std::exp(-density(h_,t,p))+helium*std::exp(-density(he_,t,p)));};
  return {rho(lo),rho(hi)};
}
double AdditiveVolumePotential::minimum_temperature()const{return std::exp(std::max(h_.t.front(),he_.t.front()));}
} // namespace ember

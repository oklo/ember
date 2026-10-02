#include "ember/phase_equilibrium.hpp"
#include "ember/constants.hpp"
#include "ember/detail/taylor3.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace ember;
using J=detail::Taylor3<5>;
namespace {
constexpr double a0=.1, gap=.6, stiffness=4, heat=1.5;
const double t0=std::log(1e5),r0=std::log(1e4);
J potential(const PhaseCoordinates&q,int phase){
  const J t=J::variable(q[0],0),r=J::variable(q[1],1);
  const J x=J::variable(q[2],2),y=J::variable(q[3],3),z=J::variable(q[4],4);
  const auto dt=t-t0;
  const auto base=exp(.6*(r-r0)-dt)+r-heat*t+x*(log(x)-1)+(y/3)*(log(y/3)-1);
  const auto a=a0+.05*dt+.02*dt*dt,offset=.2*dt+.07*dt*dt;
  if(phase<0){const auto mu=offset/gap;return base+mu*(z-a)-mu*mu/(2*stiffness);}
  const auto center=a+(phase?gap:0.);
  return base+.5*stiffness*(z-center)*(z-center)+(phase?offset:J(0));
}
PhasePotential branch(const PhaseCoordinates&q,int phase){
  const auto j=potential(q,phase);PhasePotential p;p.value=j.value();
  for(unsigned i=0;i<5;++i){J::Powers a{};a[i]++;p.gradient[i]=j.derivative(a);
    for(unsigned k=0;k<5;++k){auto b=a;b[k]++;p.hessian[i][k]=j.derivative(b);}}
  return p;
}
double relative(double a,double b){return std::abs(a-b)/std::max(1.,std::abs(b));}
}
int main(){
  unsigned failures=0;double derivative_error=0,response_error=0;
  const PhaseEvaluator first=[](const PhaseCoordinates&q){return branch(q,0);};
  const PhaseEvaluator second=[](const PhaseCoordinates&q){return branch(q,1);};
  constexpr std::array<std::array<unsigned,3>,10> powers{{
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},
    {1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}}};
  for(double dt:{-.1,0.,.1})for(double z:{.2,.4,.6}){
    const PhaseCoordinates q{t0+dt,r0+.3*dt,1e-4,1e-3,z};
    const PhaseSplit seed{q[1]+.01,q[1]-.01,{q[2],q[3],a0},{q[2],q[3],a0+gap},(z-a0)/gap};
    const auto result=phase_equilibrium_material(q,seed,first,second);
    const auto exact=potential(q,-1);HelmholtzJet material{};
    for(unsigned k=0;k<10;++k){const auto&c=powers[k];const unsigned order=c[0]+c[1]+c[2];
      for(unsigned t=0;t+order<=3;++t)for(unsigned r=0;t+r+order<=3;++r){
        const double expected=exact.derivative({t,r,c[0],c[1],c[2]});
        derivative_error=std::max(derivative_error,relative(result.jets[k][t][r]/constants::R_gas,expected));
        if(k==0)material[t][r]=constants::R_gas*expected;
      }
    }
    const auto obtained=helmholtz_response(std::exp(q[0]),std::exp(q[1]),result.jets[0]);
    const auto expected=helmholtz_response(std::exp(q[0]),std::exp(q[1]),material);
    const auto&a=obtained.state;const auto&b=expected.state;
    const std::array actual{a.P,a.E,a.S,a.cv,a.cp,a.grad_ad,a.chiT,a.chiRho,
      obtained.dcp_dlnT,obtained.dcp_dlnRho,obtained.dgrad_ad_dlnT,obtained.dgrad_ad_dlnRho,
      obtained.ddelta_dlnT,obtained.ddelta_dlnRho};
    const std::array control{b.P,b.E,b.S,b.cv,b.cp,b.grad_ad,b.chiT,b.chiRho,
      expected.dcp_dlnT,expected.dcp_dlnRho,expected.dgrad_ad_dlnT,expected.dgrad_ad_dlnRho,
      expected.ddelta_dlnT,expected.ddelta_dlnRho};
    for(std::size_t i=0;i<actual.size();++i)response_error=std::max(response_error,relative(actual[i],control[i]));
    if(!(a.cv>0&&a.chiRho>0))++failures;
    try{phase_equilibrium_material(q,seed,first,second,0);++failures;}catch(const std::invalid_argument&){}
    const PhaseEvaluator bounded=[&](const PhaseCoordinates&v){
      if(std::abs(v[0]-q[0])>5e-5)throw std::domain_error("unsupported test phase");
      return first(v);
    };
    try{phase_equilibrium_material(q,seed,bounded,second);++failures;}catch(const std::domain_error&){}
  }
  if(derivative_error>1e-5||response_error>1e-6)++failures;
  std::cout<<"phase material: derivative error "<<derivative_error<<", EOS response error "<<response_error
           <<", failures "<<failures<<'\n';
  return failures?1:0;
}

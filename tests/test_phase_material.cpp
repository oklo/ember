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
PhasePotential branch(const PhaseCoordinates&q,int phase,bool analytic=false){
  const auto j=potential(q,phase);PhasePotential p;p.value=j.value();
  for(unsigned i=0;i<5;++i){J::Powers a{};a[i]++;p.gradient[i]=j.derivative(a);
    for(unsigned k=0;k<5;++k){auto b=a;b[k]++;p.hessian[i][k]=j.derivative(b);}}
  if(analytic){p.has_third=true;
    for(unsigned i=0;i<5;++i)for(unsigned k=0;k<5;++k)for(unsigned n=0;n<5;++n){
      J::Powers a{};++a[i];++a[k];++a[n];p.third[i][k][n]=j.derivative(a);}
  }
  return p;
}
double relative(double a,double b){return std::abs(a-b)/std::max(1.,std::abs(b));}
}
int main(){
  unsigned failures=0;double derivative_error=0,response_error=0;
  constexpr std::array<std::array<unsigned,3>,10> powers{{
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},
    {1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}}};
  for(bool analytic:{false,true})for(double dt:{-.1,0.,.1})for(double z:{.2,.4,.6}){
    const PhaseEvaluator first=[analytic](const PhaseCoordinates&q){return branch(q,0,analytic);};
    const PhaseEvaluator second=[analytic](const PhaseCoordinates&q){return branch(q,1,analytic);};
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
    if(!analytic){try{phase_equilibrium_material(q,seed,bounded,second);++failures;}catch(const std::domain_error&){} }
    else{const auto local=phase_equilibrium_material(q,seed,bounded,second);
      if(!local.equilibrium.potential.has_third)++failures;}
  }
  // A vanishing mixture must join the homogeneous material without an energy
  // jump. Check both sides with an exact common-tangent construction, including
  // the large third derivatives of trace species and finite-difference mode.
  for(bool analytic:{false,true})for(int phase:{0,1}){
    const PhaseCoordinates q{t0+.05,r0,1e-8,1e-10,.4};
    const double dt=q[0]-t0,a=a0+.05*dt+.02*dt*dt;
    const double mu=(.2*dt+.07*dt*dt)/gap;
    auto boundary=q;boundary[4]=a+mu/stiffness+(phase?gap:0.);
    const PhaseEvaluator evaluator=[=](const PhaseCoordinates&x){return branch(x,phase,analytic);};
    // A single phase uses the bulk state, irrespective of an obsolete seed.
    const std::array<PhaseState,1> one{{{r0+.1,{.001,.002,.3},.2}}};
    const std::array<PhaseEvaluator,1> single{evaluator};
    const auto m=phase_equilibrium_material(boundary,one,single);
    const auto exact=potential(boundary,phase);
    for(unsigned k=0;k<10;++k){const auto&c=powers[k];const unsigned order=c[0]+c[1]+c[2];
      for(unsigned t=0;t+order<=3;++t)for(unsigned r=0;t+r+order<=3;++r)
        derivative_error=std::max(derivative_error,relative(m.jets[k][t][r]/constants::R_gas,
            exact.derivative({t,r,c[0],c[1],c[2]})));
    }
    if(m.equilibrium.phases[0].mass_fraction!=1 || m.analytic_third!=analytic)++failures;
    for(double fraction:{1e-3,1e-5,1e-7}){
      auto inside=boundary;inside[4]+=(phase?-gap:gap)*fraction;
      const PhaseEvaluator first=[=](const PhaseCoordinates&x){return branch(x,0,analytic);};
      const PhaseEvaluator second=[=](const PhaseCoordinates&x){return branch(x,1,analytic);};
      const PhaseSplit seed{r0,r0,{q[2],q[3],a+mu/stiffness},
        {q[2],q[3],a+mu/stiffness+gap},phase?1-fraction:fraction};
      const std::array<PhaseState,2> pair{{{seed.log_density_first,seed.composition_first,1-seed.second_mass_fraction},
        {seed.log_density_second,seed.composition_second,seed.second_mass_fraction}}};
      const std::array<PhaseEvaluator,2> branches{first,second};
      const auto mix=phase_equilibrium_material(inside,pair,branches);
      // F and its first derivatives approach the common tangent. Cv can jump
      // at the boundary and must not be forced to match the homogeneous value.
      for(unsigned k=0;k<5;++k)
        if(std::abs(mix.equilibrium.potential.gradient[k]-m.equilibrium.potential.gradient[k])
            >2*fraction)++failures;
      if(std::abs(mix.equilibrium.potential.value-m.equilibrium.potential.value)>fraction)++failures;
    }
  }
  if(derivative_error>1e-5||response_error>1e-6)++failures;
  std::cout<<"phase material: derivative error "<<derivative_error<<", EOS response error "<<response_error
           <<", failures "<<failures<<'\n';
  return failures?1:0;
}

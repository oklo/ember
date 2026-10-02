#include "ember/phase_equilibrium.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace ember;
namespace {
constexpr double a=.1,b=.7,kappa=4,d=.2,c=1.5;
const double t0=std::log(1e5);
PhasePotential branch(const PhaseCoordinates&q,bool second){
  const double t=q[0],r=q[1],x=q[2],y=q[3],z=q[4],center=second?b:a;
  PhasePotential p;
  p.value=r-c*t+x*(std::log(x)-1)+(y/3)*(std::log(y/3)-1)+kappa*.5*(z-center)*(z-center);
  p.gradient={-c,1,std::log(x),std::log(y/3)/3,kappa*(z-center)};
  p.hessian[2][2]=1/x;p.hessian[3][3]=1/(3*y);p.hessian[4][4]=kappa;
  if(second){p.value+=d*(t-t0);p.gradient[0]+=d;}
  return p;
}
double error(double x,double y){return std::abs(x-y)/std::max(1.,std::abs(y));}
}
int main(){
  double worst=0;unsigned failures=0;
  const PhaseEvaluator first=[](const PhaseCoordinates&q){return branch(q,false);};
  const PhaseEvaluator second=[](const PhaseCoordinates&q){return branch(q,true);};
  for(double dt:{-.2,0.,.2})for(double z:{.2,.4,.6}){
    const PhaseCoordinates q{t0+dt,std::log(1e4),1e-7,1e-9,z};
    const PhaseSplit seed{q[1]+.01,q[1]-.01,{q[2],q[3],a},{q[2],q[3],b},(z-a)/(b-a)};
    const auto result=equilibrate_phases(q,seed,first,second);
    const auto&p=result.potential;const double mu=d*dt/(b-a),slope=d/(b-a);
    // Exact common tangent of two equal-curvature composition wells.
    PhasePotential exact;
    exact.value=q[1]-c*q[0]+q[2]*(std::log(q[2])-1)+(q[3]/3)*(std::log(q[3]/3)-1)
      +mu*(z-a)-mu*mu/(2*kappa);
    exact.gradient={-c+slope*(z-a)-mu*slope/kappa,1,std::log(q[2]),std::log(q[3]/3)/3,mu};
    exact.hessian[0][0]=-slope*slope/kappa;exact.hessian[2][2]=1/q[2];exact.hessian[3][3]=1/(3*q[3]);
    exact.hessian[0][4]=exact.hessian[4][0]=slope;
    worst=std::max(worst,error(p.value,exact.value));
    for(std::size_t i=0;i<5;++i){
      worst=std::max(worst,error(p.gradient[i],exact.gradient[i]));
      for(std::size_t j=0;j<5;++j)worst=std::max(worst,error(p.hessian[i][j],exact.hessian[i][j]));
    }
    const auto&s=result.split;const double w=s.second_mass_fraction;
    if(std::abs((1-w)*std::exp(q[1]-s.log_density_first)+w*std::exp(q[1]-s.log_density_second)-1)>1e-8)++failures;
    for(std::size_t i=0;i<3;++i)
      if(std::abs(((1-w)*s.composition_first[i]+w*s.composition_second[i])/q[i+2]-1)>1e-8)++failures;
    if(error(s.composition_first[2],a+mu/kappa)>1e-8||error(s.composition_second[2],b+mu/kappa)>1e-8)++failures;
    auto invalid=q;invalid[2]=0;
    try{equilibrate_phases(invalid,seed,first,second);++failures;}catch(const std::domain_error&){}
    const PhaseEvaluator bad=[](const PhaseCoordinates&x){auto p=branch(x,false);
      p.hessian[2][2]=std::numeric_limits<double>::quiet_NaN();return p;};
    try{equilibrate_phases(q,seed,bad,second);++failures;}catch(const std::domain_error&){}
  }
  if(worst>1e-8)++failures;
  std::cout<<"phase equilibrium: analytic potential/response error "<<worst<<", failures "<<failures<<'\n';
  return failures?1:0;
}

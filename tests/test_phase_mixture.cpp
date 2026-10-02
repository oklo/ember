#include "ember/phase_equilibrium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace ember;
namespace {
constexpr double kappa=4,heat=1.5;
constexpr std::array<std::array<double,2>,3> centers{{{.1,.1},{.4,.1},{.1,.6}}};
constexpr std::array<double,2> moveT{.02,.03},moveR{.01,-.02};
const double t0=std::log(1e5),r0=std::log(1e4);
double dot(const std::array<double,2>&a,const std::array<double,2>&b){return a[0]*b[0]+a[1]*b[1];}
double error(double x,double y){return std::abs(x-y)/std::max(1.,std::abs(y));}
PhasePotential branch(const PhaseCoordinates&q,unsigned phase,const std::array<double,2>&slope){
  const double dt=q[0]-t0,dr=q[1]-r0,y=q[3];std::array<double,2> d{};double shift=0;
  for(unsigned j=0;j<2;++j){d[j]=q[2+2*j]-centers[phase][j]-moveT[j]*dt-moveR[j]*dr;
    shift+=(centers[phase][j]-centers[0][j])*slope[j];}
  PhasePotential p;
  p.value=q[1]-heat*q[0]+y/3*(std::log(y/3)-1)+.5*kappa*dot(d,d)+shift*dt;
  p.gradient={-heat-kappa*dot(d,moveT)+shift,1-kappa*dot(d,moveR),kappa*d[0],std::log(y/3)/3,kappa*d[1]};
  p.hessian[0][0]=kappa*dot(moveT,moveT);p.hessian[1][1]=kappa*dot(moveR,moveR);
  p.hessian[0][1]=p.hessian[1][0]=kappa*dot(moveT,moveR);p.hessian[3][3]=1/(3*y);
  for(unsigned j=0;j<2;++j){const auto i=2+2*j;p.hessian[i][i]=kappa;
    p.hessian[0][i]=p.hessian[i][0]=-kappa*moveT[j];p.hessian[1][i]=p.hessian[i][1]=-kappa*moveR[j];}
  p.value+=.1*dt*dt*dt+.03*dt*dt*dr+.02*dt*dr*dr+.01*dr*dr*dr;
  p.gradient[0]+=.3*dt*dt+.06*dt*dr+.02*dr*dr;
  p.gradient[1]+=.03*dt*dt+.04*dt*dr+.03*dr*dr;
  p.hessian[0][0]+=.6*dt+.06*dr;p.hessian[1][1]+=.04*dt+.06*dr;
  p.hessian[0][1]+=.06*dt+.04*dr;p.hessian[1][0]+=.06*dt+.04*dr;
  p.third[0][0][0]=.6;p.third[1][1][1]=.06;
  p.third[0][0][1]=p.third[0][1][0]=p.third[1][0][0]=.06;
  p.third[0][1][1]=p.third[1][0][1]=p.third[1][1][0]=.04;
  p.third[3][3][3]=-1/(3*y*y);
  return p;
}
}
int main(){
  double worst=0,material_error=0,third_error=0,stiff_third_error=0;unsigned failures=0,cases=0;
  for(bool analytic:{false,true})for(double scale:{1.,1e4})for(double dt:{-.15/scale,0.,.15/scale})for(double dr:{-.1,0.,.1}){
    ++cases;const std::array<double,2>slope{.3*scale,-.2*scale};
    const PhaseCoordinates q{t0+dt,r0+dr,.21,1e-9,.24};
    const double actualdt=q[0]-t0,actualdr=q[1]-r0;std::array<double,2>mu{},shift{},d{};
    for(unsigned j=0;j<2;++j){mu[j]=slope[j]*actualdt;shift[j]=moveT[j]*actualdt+moveR[j]*actualdr;
      d[j]=q[2+2*j]-centers[0][j]-shift[j];}
    std::array<double,3>w{0,(d[0]-mu[0]/kappa)/.3,(d[1]-mu[1]/kappa)/.5};w[0]=1-w[1]-w[2];
    std::array<PhaseState,3>guess;std::array<PhaseEvaluator,3>evaluators;
    for(unsigned i=0;i<3;++i){
      guess[i]={q[1]+.003*(int(i)-1),{(centers[i][0]+shift[0]+mu[0]/kappa)*1.001,q[3]*1.001,
        (centers[i][1]+shift[1]+mu[1]/kappa)*1.001},w[i]};
      evaluators[i]=[i,slope,analytic](const PhaseCoordinates&x){auto p=branch(x,i,slope);p.has_third=analytic;return p;};
    }
    const auto result=equilibrate_phases(q,guess,evaluators);const auto&p=result.potential;
    PhasePotential exact;
    exact.value=q[1]-heat*q[0]+q[3]/3*(std::log(q[3]/3)-1)+dot(mu,d)-dot(mu,mu)/(2*kappa);
    exact.gradient={-heat+dot(slope,d)-dot(mu,moveT)-dot(mu,slope)/kappa,1-dot(mu,moveR),mu[0],std::log(q[3]/3)/3,mu[1]};
    exact.hessian[0][0]=-2*dot(slope,moveT)-dot(slope,slope)/kappa;
    exact.hessian[0][1]=exact.hessian[1][0]=-dot(slope,moveR);
    exact.hessian[0][2]=exact.hessian[2][0]=slope[0];exact.hessian[0][4]=exact.hessian[4][0]=slope[1];
    exact.hessian[3][3]=1/(3*q[3]);
    const auto common=branch(q,0,slope);
    exact.value+=.1*actualdt*actualdt*actualdt+.03*actualdt*actualdt*actualdr+.02*actualdt*actualdr*actualdr+.01*actualdr*actualdr*actualdr;
    exact.gradient[0]+=.3*actualdt*actualdt+.06*actualdt*actualdr+.02*actualdr*actualdr;
    exact.gradient[1]+=.03*actualdt*actualdt+.04*actualdt*actualdr+.03*actualdr*actualdr;
    exact.hessian[0][0]+=.6*actualdt+.06*actualdr;exact.hessian[1][1]+=.04*actualdt+.06*actualdr;
    exact.hessian[0][1]+=.06*actualdt+.04*actualdr;exact.hessian[1][0]+=.06*actualdt+.04*actualdr;
    worst=std::max(worst,error(p.value,exact.value));
    for(unsigned i=0;i<5;++i){worst=std::max(worst,error(p.gradient[i],exact.gradient[i]));
      for(unsigned j=0;j<5;++j)worst=std::max(worst,error(p.hessian[i][j],exact.hessian[i][j]));}
    double volume=0;std::array<double,3>inventory{};
    for(unsigned i=0;i<3;++i){const auto&s=result.phases[i];
      worst=std::max(worst,error(s.mass_fraction,w[i]));volume+=s.mass_fraction*std::exp(q[1]-s.log_density);
      for(unsigned j=0;j<3;++j)inventory[j]+=s.mass_fraction*s.composition[j];
    }
    worst=std::max(worst,std::abs(volume-1));
    for(unsigned j=0;j<3;++j)worst=std::max(worst,std::abs(inventory[j]/q[j+2]-1));
    for(unsigned i=1;i<3;++i){const double width=i==1?.3:.5;
      const PhaseCoordinates dw{(-moveT[i-1]-slope[i-1]/kappa)/width,-moveR[i-1]/width,i==1?1/width:0,0,i==2?1/width:0};
      for(unsigned j=0;j<5;++j)worst=std::max(worst,error(result.fraction_response[i][j],dw[j]));}
    const auto material=phase_equilibrium_material(q,guess,evaluators);
    material_error=std::max(material_error,error(material.jets[0][2][0]/constants::R_gas,exact.hessian[0][0]));
    // With stiff thermal partitioning, a fixed 1e-4 stencil exits coexistence.
    if(scale>1&&!analytic&&!(material.log_steps[0]<1e-4&&material.log_steps[0]>0))++failures;
    if(material.analytic_third!=analytic)++failures;
    if(analytic){
      for(unsigned i=0;i<5;++i)for(unsigned j=0;j<5;++j)for(unsigned k=0;k<5;++k)
        {
          // A trace-species derivative is measured for a relative abundance
          // change; a zero cross derivative divided by a tiny X is ill-scaled.
          const double weight=(i<2?1:q[i])*(j<2?1:q[j])*(k<2?1:q[k]);
          const double e=error(material.equilibrium.potential.third[i][j][k]*weight,common.third[i][j][k]*weight);
          if(scale==1)third_error=std::max(third_error,e);else stiff_third_error=std::max(stiff_third_error,e);
        }
      if(material.log_steps[0]!=0||material.log_steps[1]!=0)++failures;
    }
    std::rotate(guess.begin(),guess.begin()+1,guess.end());std::rotate(evaluators.begin(),evaluators.begin()+1,evaluators.end());
    const auto reordered=equilibrate_phases(q,guess,evaluators);
    worst=std::max(worst,error(reordered.potential.value,p.value));
    for(unsigned i=0;i<5;++i)worst=std::max(worst,error(reordered.potential.gradient[i],p.gradient[i]));
    auto invalid=q;invalid[3]=0;
    try{equilibrate_phases(invalid,guess,evaluators);++failures;}catch(const std::domain_error&){}
  }
  // In the stiff case, order-unity third derivatives cancel terms of order
  // 1e12. Test them to the intended 0.1% physical accuracy, alongside the
  // well-conditioned case, where roundoff is much smaller.
  if(worst>2e-8||material_error>2e-8||third_error>1e-8||stiff_third_error>1e-3)++failures;
  std::cout<<"three-phase cases "<<cases<<", analytic error "<<worst<<", material error "<<material_error<<", third derivative error "<<third_error<<", stiff third derivative error "<<stiff_third_error<<", failures "<<failures<<'\n';
  return failures?1:0;
}

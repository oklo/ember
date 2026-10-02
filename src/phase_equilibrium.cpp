#include "ember/phase_equilibrium.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
using V=std::array<double,9>;
using M=std::array<V,9>;
using Responses=std::array<PhaseCoordinates,9>;
struct System {
  V residual{};
  M jacobian{};
  Responses external{};
  PhasePotential a,b;
  std::array<double,3> ca{},cb{};
};
double norm(const V&v){double value=0;for(double x:v)value=std::max(value,std::abs(x));return value;}
bool valid(const V&u){
  for(double x:u)if(!std::isfinite(x))return false;
  if(!(u[8]>0&&u[8]<1))return false;
  double a=0,b=0;
  for(std::size_t k=0;k<3;++k){
    if(u[2+k]<-600||u[5+k]<-600||u[2+k]>=0||u[5+k]>=0)return false;
    a+=std::exp(u[2+k]);b+=std::exp(u[5+k]);
  }
  return a<1&&b<1;
}
template<std::size_t N> std::array<std::array<double,N>,9> linear(
    M a,std::array<std::array<double,N>,9> b){
  for(std::size_t i=0;i<9;++i){
    double scale=0;for(double x:a[i])scale=std::max(scale,std::abs(x));
    if(!(scale>0&&std::isfinite(scale)))throw std::domain_error("phase equilibrium: singular constraints");
    for(auto&x:a[i])x/=scale;for(auto&x:b[i])x/=scale;
  }
  for(std::size_t k=0;k<9;++k){
    std::size_t pivot=k;for(std::size_t i=k+1;i<9;++i)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
    if(std::abs(a[pivot][k])<64*std::numeric_limits<double>::epsilon())
      throw std::domain_error("phase equilibrium: unresolved or coalescing phases");
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    for(std::size_t i=k+1;i<9;++i){
      const double f=a[i][k]/a[k][k];
      for(std::size_t j=k+1;j<9;++j)a[i][j]-=f*a[k][j];
      for(std::size_t j=0;j<N;++j)b[i][j]-=f*b[k][j];
    }
  }
  for(std::size_t ii=9;ii-->0;)for(std::size_t k=0;k<N;++k){
    for(std::size_t j=ii+1;j<9;++j)b[ii][k]-=a[ii][j]*b[j][k];
    b[ii][k]/=a[ii][ii];
  }
  return b;
}
System assemble(const PhaseCoordinates&q,const V&u,const PhaseEvaluator&first,const PhaseEvaluator&second){
  System s;auto qa=q,qb=q;qa[1]=u[0];qb[1]=u[1];
  for(std::size_t i=0;i<3;++i){s.ca[i]=qa[i+2]=std::exp(u[i+2]);s.cb[i]=qb[i+2]=std::exp(u[i+5]);}
  s.a=first(qa);s.b=second(qb);
  for(const auto*p:{&s.a,&s.b}){
    bool finite=std::isfinite(p->value);
    for(double x:p->gradient)finite=finite&&std::isfinite(x);
    for(const auto&row:p->hessian)for(double x:row)finite=finite&&std::isfinite(x);
    if(!finite)throw std::domain_error("phase equilibrium: nonfinite material derivative");
  }
  const auto&ga=s.a.gradient;const auto&gb=s.b.gradient;
  const auto&ha=s.a.hessian;const auto&hb=s.b.hessian;
  if(!(ga[1]>0&&gb[1]>0))throw std::domain_error("phase equilibrium: nonpositive phase pressure");
  const double w=u[8],wa=1-w,va=std::exp(q[1]-u[0]),vb=std::exp(q[1]-u[1]);
  auto&f=s.residual;auto&j=s.jacobian;auto&b=s.external;
  f[0]=u[0]-u[1]+std::log(ga[1]/gb[1]);
  f[4]=s.a.value+ga[1]-s.b.value-gb[1];
  f[5]=wa*va+w*vb-1;
  j[0][0]=1+ha[1][1]/ga[1];j[0][1]=-1-hb[1][1]/gb[1];
  j[4][0]=ga[1]+ha[1][1];j[4][1]=-gb[1]-hb[1][1];
  j[5][0]=-wa*va;j[5][1]=-w*vb;j[5][8]=vb-va;
  b[0][0]=ha[1][0]/ga[1]-hb[1][0]/gb[1];
  b[4][0]=ga[0]+ha[1][0]-gb[0]-hb[1][0];b[5][1]=wa*va+w*vb;
  for(std::size_t i=0;i<3;++i){
    const double ca=s.ca[i],cb=s.cb[i],total=wa*ca+w*cb;const auto k=i+2;
    f[i+1]=ga[k]-gb[k];f[4]-=ca*ga[k]-cb*gb[k];f[i+6]=total/q[k]-1;
    j[0][i+2]=ha[1][k]*ca/ga[1];j[0][i+5]=-hb[1][k]*cb/gb[1];
    j[i+1][0]=ha[k][1];j[i+1][1]=-hb[k][1];
    j[4][0]-=ca*ha[k][1];j[4][1]+=cb*hb[k][1];
    j[4][i+2]=ha[1][k]*ca;j[4][i+5]=-hb[1][k]*cb;
    j[i+6][i+2]=wa*ca/q[k];j[i+6][i+5]=w*cb/q[k];j[i+6][8]=(cb-ca)/q[k];
    b[i+1][0]=ha[k][0]-hb[k][0];b[4][0]-=ca*ha[k][0]-cb*hb[k][0];
    b[i+6][k]=-total/(q[k]*q[k]);
    for(std::size_t n=0;n<3;++n){
      j[i+1][n+2]=ha[k][n+2]*s.ca[n];j[i+1][n+5]=-hb[k][n+2]*s.cb[n];
      j[4][i+2]-=s.ca[n]*ha[n+2][k]*ca;j[4][i+5]+=s.cb[n]*hb[n+2][k]*cb;
    }
  }
  for(double x:f)if(!std::isfinite(x))throw std::domain_error("phase equilibrium: nonfinite phase response");
  return s;
}
} // namespace

PhaseEquilibrium equilibrate_phases(const PhaseCoordinates&q,const PhaseSplit&guess,
    const PhaseEvaluator&first,const PhaseEvaluator&second){
  for(double x:q)if(!std::isfinite(x))throw std::invalid_argument("phase equilibrium: invalid state");
  if(!(q[2]>0&&q[3]>0&&q[4]>0&&q[2]+q[3]+q[4]<1))
    throw std::domain_error("phase equilibrium: positive composition required");
  V u{};u[0]=guess.log_density_first;u[1]=guess.log_density_second;u[8]=guess.second_mass_fraction;
  for(std::size_t i=0;i<3;++i){u[i+2]=std::log(guess.composition_first[i]);u[i+5]=std::log(guess.composition_second[i]);}
  if(!valid(u))throw std::invalid_argument("phase equilibrium: invalid initial split");
  System s;std::size_t iteration=0;
  for(;iteration<40;++iteration){
    s=assemble(q,u,first,second);const double error=norm(s.residual);
    if(error<2e-9)break;
    std::array<std::array<double,1>,9> rhs{};for(std::size_t i=0;i<9;++i)rhs[i][0]=-s.residual[i];
    const auto delta=linear(s.jacobian,rhs);double scale=1;bool accepted=false;
    for(unsigned trial=0;trial<22;++trial,scale*=.5){
      V next=u;for(std::size_t i=0;i<9;++i)next[i]+=scale*delta[i][0];
      if(!valid(next))continue;
      try {if(norm(assemble(q,next,first,second).residual)<error){u=next;accepted=true;break;}}
      catch(const std::domain_error&){ } // shorten an unsupported trial; accepted states must pass
    }
    if(!accepted)throw std::runtime_error("phase equilibrium: Newton step did not decrease residual");
  }
  if(iteration==40)throw std::runtime_error("phase equilibrium: iteration limit");
  Responses rhs=s.external;for(auto&row:rhs)for(auto&x:row)x=-x;
  const auto response=linear(s.jacobian,rhs);
  const auto&ga=s.a.gradient;const auto&gb=s.b.gradient;const auto&ha=s.a.hessian;const auto&hb=s.b.hessian;
  const double w=u[8],wa=1-w,v=std::exp(u[0]-q[1]);
  PhaseEquilibrium result;auto&p=result.potential;
  result.residual=norm(s.residual);result.iterations=iteration;
  result.split={u[0],u[1],s.ca,s.cb,w};p.value=wa*s.a.value+w*s.b.value;
  p.gradient[0]=wa*ga[0]+w*gb[0];p.gradient[1]=v*ga[1];
  std::array<V,5> internal{};
  internal[0][0]=wa*ha[0][1];internal[0][1]=w*hb[0][1];internal[0][8]=gb[0]-ga[0];
  internal[1][0]=v*(ga[1]+ha[1][1]);
  p.hessian[0][0]=wa*ha[0][0]+w*hb[0][0];p.hessian[1][0]=v*ha[1][0];p.hessian[1][1]=-p.gradient[1];
  for(std::size_t i=0;i<3;++i){
    p.gradient[i+2]=ga[i+2];p.hessian[i+2][0]=ha[i+2][0];internal[i+2][0]=ha[i+2][1];
    internal[0][i+2]=wa*ha[0][i+2]*s.ca[i];internal[0][i+5]=w*hb[0][i+2]*s.cb[i];
    internal[1][i+2]=v*ha[1][i+2]*s.ca[i];
    for(std::size_t n=0;n<3;++n)internal[i+2][n+2]=ha[i+2][n+2]*s.ca[n];
  }
  for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)for(std::size_t k=0;k<9;++k)
    p.hessian[i][j]+=internal[i][k]*response[k][j];
  // The two derivative routes are identical analytically. Average their
  // roundoff, especially when a trace-species coordinate is very small.
  for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<i;++j)
    p.hessian[i][j]=p.hessian[j][i]=.5*(p.hessian[i][j]+p.hessian[j][i]);
  return result;
}

PhaseMaterial phase_equilibrium_material(const PhaseCoordinates&q,const PhaseSplit&guess,
    const PhaseEvaluator&first,const PhaseEvaluator&second,double step){
  if(!(std::isfinite(step)&&step>0))throw std::invalid_argument("phase material: invalid difference step");
  PhaseMaterial result;result.equilibrium=equilibrate_phases(q,guess,first,second);
  const auto&p=result.equilibrium.potential;
  std::array<std::array<PhaseCoordinates,5>,2> third{};
  for(std::size_t k=0;k<2;++k){
    auto plus=q,minus=q;plus[k]+=step;minus[k]-=step;
    if(plus[k]==q[k]||minus[k]==q[k])throw std::invalid_argument("phase material: unresolved difference step");
    const auto a=equilibrate_phases(plus,result.equilibrium.split,first,second);
    const auto b=equilibrate_phases(minus,result.equilibrium.split,first,second);
    for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)
      third[k][i][j]=(a.potential.hessian[i][j]-b.potential.hessian[i][j])/(2*step);
  }
  constexpr std::array<std::array<unsigned,3>,10> powers{{
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},
    {1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}}};
  for(std::size_t k=0;k<10;++k){
    const auto&c=powers[k];const unsigned order=c[0]+c[1]+c[2];
    for(unsigned t=0;t+order<=3;++t)for(unsigned r=0;t+r+order<=3;++r){
      std::array<unsigned,3> indices{};unsigned n=0;
      for(unsigned j=0;j<t;++j)indices[n++]=0;
      for(unsigned j=0;j<r;++j)indices[n++]=1;
      for(unsigned j=0;j<3;++j)for(unsigned l=0;l<c[j];++l)indices[n++]=j+2;
      double value=p.value;
      if(n==1)value=p.gradient[indices[0]];
      if(n==2)value=p.hessian[indices[0]][indices[1]];
      if(n==3){
        // Mixed derivatives have two independently evaluated routes. Average
        // only those available from a thermal/density Hessian difference.
        value=0;unsigned routes=0;
        for(unsigned j=0;j<3;++j)if(indices[j]<2){
          value+=third[indices[j]][indices[(j+1)%3]][indices[(j+2)%3]];++routes;
        }
        value/=routes;
      }
      value*=constants::R_gas;
      if(!std::isfinite(value))throw std::domain_error("phase material: nonfinite response");
      result.jets[k][t][r]=value;
    }
  }
  return result;
}
} // namespace ember

#include "ember/phase_equilibrium.hpp"
#include "ember/constants.hpp"
#include "ember/detail/taylor3.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace ember {
namespace {
template<std::size_t N> using Vector=std::array<double,N>;
template<std::size_t N> using Matrix=std::array<Vector<N>,N>;
template<std::size_t N> double norm(const Vector<N>&v){
  double value=0;for(double x:v)value=std::max(value,std::abs(x));return value;
}
template<std::size_t S,std::size_t K> std::array<Vector<K>,S> linear(
    Matrix<S> a,std::array<Vector<K>,S> b){
  for(std::size_t i=0;i<S;++i){
    double scale=norm(a[i]);
    if(!(scale>0&&std::isfinite(scale)))throw std::domain_error("phase equilibrium: singular constraints");
    for(auto&x:a[i])x/=scale;for(auto&x:b[i])x/=scale;
  }
  for(std::size_t k=0;k<S;++k){
    std::size_t pivot=k;for(std::size_t i=k+1;i<S;++i)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
    if(std::abs(a[pivot][k])<64*std::numeric_limits<double>::epsilon())
      throw std::domain_error("phase equilibrium: unresolved or coalescing phases");
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    for(std::size_t i=k+1;i<S;++i){
      const double f=a[i][k]/a[k][k];
      for(std::size_t j=k+1;j<S;++j)a[i][j]-=f*a[k][j];
      for(std::size_t j=0;j<K;++j)b[i][j]-=f*b[k][j];
    }
  }
  for(std::size_t ii=S;ii-->0;)for(std::size_t k=0;k<K;++k){
    for(std::size_t j=ii+1;j<S;++j)b[ii][k]-=a[ii][j]*b[j][k];
    b[ii][k]/=a[ii][ii];
  }
  return b;
}
template<std::size_t N> struct System {
  static constexpr std::size_t size=5*N-1;
  Vector<size> residual{};
  Matrix<size> jacobian{};
  std::array<PhaseCoordinates,size> external{};
  std::array<PhasePotential,N> potential;
  std::array<Vector<3>,N> composition{};
  Vector<N> weight{},volume{};
};
template<std::size_t N> bool valid(const Vector<5*N-1>&u){
  for(double x:u)if(!std::isfinite(x))return false;
  double total=0;
  for(std::size_t k=4*N;k<u.size();++k){if(u[k]<=0)return false;total+=u[k];}
  if(total>=1)return false;
  for(std::size_t k=0;k<N;++k){
    double sum=0;
    for(std::size_t j=1;j<4;++j){
      const double v=u[4*k+j];if(v<-600||v>=0)return false;sum+=std::exp(v);
    }
    if(sum>=1)return false;
  }
  return true;
}
void check_potential(const PhasePotential&p){
  bool finite=std::isfinite(p.value);
  for(double x:p.gradient)finite=finite&&std::isfinite(x);
  for(const auto&row:p.hessian)for(double x:row)finite=finite&&std::isfinite(x);
  if(p.has_third)for(const auto&plane:p.third)for(const auto&row:plane)
    for(double x:row)finite=finite&&std::isfinite(x);
  if(!finite)throw std::domain_error("phase equilibrium: nonfinite material derivative");
  if(!(p.gradient[1]>0))throw std::domain_error("phase equilibrium: nonpositive phase pressure");
}
template<std::size_t N> System<N> assemble(const PhaseCoordinates&q,const Vector<5*N-1>&u,
                                         std::span<const PhaseEvaluator> evaluators){
  System<N> s;s.weight[0]=1;
  for(std::size_t k=0;k<N;++k){
    auto point=q;point[1]=u[4*k];
    for(std::size_t j=0;j<3;++j)s.composition[k][j]=point[j+2]=std::exp(u[4*k+j+1]);
    s.potential[k]=evaluators[k](point);check_potential(s.potential[k]);
    s.volume[k]=std::exp(q[1]-point[1]);
    if(k){s.weight[k]=u[4*N+k-1];s.weight[0]-=s.weight[k];}
  }
  auto&f=s.residual;auto&j=s.jacobian;auto&b=s.external;
  const auto&c0=s.composition[0];const auto&p0=s.potential[0];const auto&g0=p0.gradient;const auto&h0=p0.hessian;
  for(std::size_t k=1;k<N;++k){
    const auto&c=s.composition[k];const auto&p=s.potential[k];const auto&g=p.gradient;const auto&h=p.hessian;
    const auto row=5*(k-1),col=4*k;
    f[row]=u[0]-u[col]+std::log(g0[1]/g[1]);
    j[row][0]=1+h0[1][1]/g0[1];j[row][col]=-1-h[1][1]/g[1];
    b[row][0]=h0[1][0]/g0[1]-h[1][0]/g[1];
    f[row+4]=p0.value+g0[1]-p.value-g[1];
    j[row+4][0]=g0[1]+h0[1][1];j[row+4][col]=-g[1]-h[1][1];
    b[row+4][0]=g0[0]+h0[1][0]-g[0]-h[1][0];
    for(std::size_t a=0;a<3;++a){
      const auto i=a+2;
      j[row][a+1]=h0[1][i]*c0[a]/g0[1];j[row][col+a+1]=-h[1][i]*c[a]/g[1];
      f[row+a+1]=g0[i]-g[i];j[row+a+1][0]=h0[i][1];j[row+a+1][col]=-h[i][1];
      b[row+a+1][0]=h0[i][0]-h[i][0];
      f[row+4]-=c0[a]*g0[i]-c[a]*g[i];
      j[row+4][0]-=c0[a]*h0[i][1];j[row+4][col]+=c[a]*h[i][1];
      b[row+4][0]-=c0[a]*h0[i][0]-c[a]*h[i][0];
      j[row+4][a+1]=h0[1][i]*c0[a];j[row+4][col+a+1]=-h[1][i]*c[a];
      for(std::size_t z=0;z<3;++z){
        j[row+a+1][z+1]=h0[i][z+2]*c0[z];j[row+a+1][col+z+1]=-h[i][z+2]*c[z];
        j[row+4][a+1]-=c0[z]*h0[z+2][i]*c0[a];j[row+4][col+a+1]+=c[z]*h[z+2][i]*c[a];
      }
    }
  }
  constexpr std::size_t vol=5*(N-1);
  f[vol]=-1;for(std::size_t a=0;a<3;++a)f[vol+1+a]=-1;
  for(std::size_t k=0;k<N;++k){
    const auto&c=s.composition[k];const double w=s.weight[k],v=s.volume[k];
    f[vol]+=w*v;b[vol][1]+=w*v;j[vol][4*k]=-w*v;
    if(k)j[vol][4*N+k-1]=v-s.volume[0];
    for(std::size_t a=0;a<3;++a){
      f[vol+1+a]+=w*c[a]/q[a+2];b[vol+1+a][a+2]-=w*c[a]/(q[a+2]*q[a+2]);
      j[vol+1+a][4*k+a+1]=w*c[a]/q[a+2];
      if(k)j[vol+1+a][4*N+k-1]=(c[a]-c0[a])/q[a+2];
    }
  }
  for(double x:f)if(!std::isfinite(x))throw std::domain_error("phase equilibrium: nonfinite phase response");
  return s;
}
template<std::size_t N> PhaseMixture solve(const PhaseCoordinates&q,std::span<const PhaseState> guess,
                                          std::span<const PhaseEvaluator> evaluators){
  constexpr std::size_t size=5*N-1;Vector<size> u{};
  double sum=0;
  for(std::size_t k=0;k<N;++k){
    if(!(guess[k].mass_fraction>0))throw std::invalid_argument("phase equilibrium: positive phase fractions required");
    u[4*k]=guess[k].log_density;sum+=guess[k].mass_fraction;
    for(std::size_t j=0;j<3;++j)u[4*k+j+1]=std::log(guess[k].composition[j]);
    if(k)u[4*N+k-1]=guess[k].mass_fraction;
  }
  if(!valid<N>(u)||!std::isfinite(sum)||std::abs(sum-1)>1e-8)
    throw std::invalid_argument("phase equilibrium: invalid initial split");
  System<N> s;std::size_t iteration=0;
  for(;iteration<50;++iteration){
    s=assemble<N>(q,u,evaluators);const double error=norm(s.residual);
    if(error<2e-9)break;
    std::array<Vector<1>,size> rhs{};for(std::size_t i=0;i<size;++i)rhs[i][0]=-s.residual[i];
    const auto delta=linear(s.jacobian,rhs);double scale=1;bool accepted=false;
    for(unsigned trial=0;trial<24;++trial,scale*=.5){
      auto next=u;for(std::size_t i=0;i<size;++i)next[i]+=scale*delta[i][0];
      if(!valid<N>(next))continue;
      try{if(norm(assemble<N>(q,next,evaluators).residual)<error){u=next;accepted=true;break;}}
      catch(const std::domain_error&){ } // shorten an unsupported trial; accepted states must pass
    }
    if(!accepted)throw std::runtime_error("phase equilibrium: Newton step did not decrease residual");
  }
  if(iteration==50)throw std::runtime_error("phase equilibrium: iteration limit");
  auto rhs=s.external;for(auto&row:rhs)for(auto&x:row)x=-x;
  const auto response=linear(s.jacobian,rhs);
  PhaseMixture result;result.residual=norm(s.residual);result.iterations=iteration;
  result.phases.resize(N);result.fraction_response.resize(N);result.state_response.resize(N);
  const auto&g0=s.potential[0].gradient;const auto&h0=s.potential[0].hessian;const auto&c0=s.composition[0];
  const double v=std::exp(u[0]-q[1]);auto&p=result.potential;
  std::array<Vector<size>,5> internal{};
  for(std::size_t k=0;k<N;++k){
    const auto&phase=s.potential[k];const auto&g=phase.gradient;const auto&h=phase.hessian;const double w=s.weight[k];
    result.phases[k]={u[4*k],s.composition[k],w};
    for(std::size_t a=0;a<4;++a)result.state_response[k][a]=response[4*k+a];
    p.value+=w*phase.value;p.gradient[0]+=w*g[0];p.hessian[0][0]+=w*h[0][0];
    internal[0][4*k]=w*h[0][1];
    for(std::size_t a=0;a<3;++a)internal[0][4*k+a+1]=w*h[0][a+2]*s.composition[k][a];
    if(k){
      internal[0][4*N+k-1]=g[0]-g0[0];
      result.fraction_response[k]=response[4*N+k-1];
      for(std::size_t a=0;a<5;++a)result.fraction_response[0][a]-=response[4*N+k-1][a];
    }
  }
  p.gradient[1]=v*g0[1];p.hessian[1][0]=v*h0[1][0];p.hessian[1][1]=-p.gradient[1];
  internal[1][0]=v*(g0[1]+h0[1][1]);
  for(std::size_t a=0;a<3;++a){
    p.gradient[a+2]=g0[a+2];p.hessian[a+2][0]=h0[a+2][0];internal[a+2][0]=h0[a+2][1];
    internal[1][a+1]=v*h0[1][a+2]*c0[a];
    for(std::size_t z=0;z<3;++z)internal[a+2][z+1]=h0[a+2][z+2]*c0[z];
  }
  for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)for(std::size_t k=0;k<size;++k)
    p.hessian[i][j]+=internal[i][k]*response[k][j];
  for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<i;++j)
    p.hessian[i][j]=p.hessian[j][i]=.5*(p.hessian[i][j]+p.hessian[j][i]);
  check_potential(p);
  return result;
}
PhaseEquilibrium pair_result(const PhaseMixture&r){
  const auto&a=r.phases[0];const auto&b=r.phases[1];
  return {r.potential,{a.log_density,b.log_density,a.composition,b.composition,b.mass_fraction},r.residual,r.iterations};
}
std::array<PhaseState,2> pair_guess(const PhaseSplit&s){
  return {{{s.log_density_first,s.composition_first,1-s.second_mass_fraction},
           {s.log_density_second,s.composition_second,s.second_mass_fraction}}};
}
std::array<HelmholtzJet,10> material_jets(const PhasePotential&p,
    const std::array<std::array<PhaseCoordinates,5>,2>&third){
  std::array<HelmholtzJet,10> jets{};
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
        value=0;unsigned routes=0;
        for(unsigned j=0;j<3;++j)if(indices[j]<2){
          value+=third[indices[j]][indices[(j+1)%3]][indices[(j+2)%3]];++routes;
        }
        value/=routes;
      }
      value*=constants::R_gas;
      if(!std::isfinite(value))throw std::domain_error("phase material: nonfinite response");
      jets[k][t][r]=value;
    }
  }
  return jets;
}
using Third=std::array<std::array<PhaseCoordinates,5>,5>;
using Jet=detail::Taylor3<5>;
struct ExpandedPhase {Jet potential;std::array<Jet,5> gradient;};
ExpandedPhase expand(const PhasePotential&p,const std::array<Jet,5>&delta){
  ExpandedPhase result;result.potential=Jet(p.value);
  for(std::size_t i=0;i<5;++i){
    result.gradient[i]=Jet(p.gradient[i]);result.potential=result.potential+p.gradient[i]*delta[i];
    for(std::size_t j=0;j<5;++j){
      result.potential=result.potential+.5*p.hessian[i][j]*delta[i]*delta[j];
      result.gradient[i]=result.gradient[i]+p.hessian[i][j]*delta[j];
      for(std::size_t k=0;k<5;++k)
        result.gradient[i]=result.gradient[i]+.5*p.third[i][j][k]*delta[j]*delta[k];
    }
  }
  return result;
}
template<std::size_t N> std::optional<Third> analytic_material(const PhaseCoordinates&q,
    const PhaseMixture&mixture,std::span<const PhaseEvaluator>evaluators){
  constexpr std::size_t size=5*N-1;Vector<size> u{};std::array<Jet,size> variable;
  const PhaseCoordinates scale{1,1,q[2],q[3],q[4]};std::array<Jet,5> bulk;
  for(std::size_t j=0;j<5;++j)bulk[j]=Jet(q[j])+Jet::variable(0,j)*scale[j];
  auto variable_with_response=[&](double value,const PhaseCoordinates&response){
    Jet v(value);for(std::size_t j=0;j<5;++j)v=v+Jet::variable(0,j)*(response[j]*scale[j]);return v;
  };
  for(std::size_t k=0;k<N;++k){
    const auto&phase=mixture.phases[k];u[4*k]=phase.log_density;
    for(std::size_t a=0;a<3;++a)u[4*k+a+1]=std::log(phase.composition[a]);
    for(std::size_t a=0;a<4;++a)variable[4*k+a]=variable_with_response(u[4*k+a],mixture.state_response[k][a]);
    if(k){u[4*N+k-1]=phase.mass_fraction;variable[4*N+k-1]=variable_with_response(phase.mass_fraction,mixture.fraction_response[k]);}
  }
  const auto s=assemble<N>(q,u,evaluators);
  for(const auto&p:s.potential)if(!p.has_third)return std::nullopt;
  std::array<Jet,N> weight;std::array<std::array<Jet,3>,N> composition;std::array<ExpandedPhase,N> phases;
  auto evaluate=[&](){
    weight[0]=Jet(1);
    for(std::size_t k=0;k<N;++k){
      if(k){weight[k]=variable[4*N+k-1];weight[0]=weight[0]-weight[k];}
      std::array<Jet,5> delta{bulk[0]-q[0],variable[4*k]-u[4*k],Jet(0),Jet(0),Jet(0)};
      for(std::size_t a=0;a<3;++a){composition[k][a]=detail::exp(variable[4*k+a+1]);delta[a+2]=composition[k][a]-s.composition[k][a];}
      phases[k]=expand(s.potential[k],delta);
    }
  };
  evaluate();std::array<Jet,size> residual;
  for(std::size_t k=1;k<N;++k){
    const std::size_t row=5*(k-1);const auto&a=phases[0];const auto&b=phases[k];
    residual[row]=variable[0]-variable[4*k]+detail::log(a.gradient[1]/b.gradient[1]);
    residual[row+4]=a.potential+a.gradient[1]-b.potential-b.gradient[1];
    for(std::size_t z=0;z<3;++z){
      residual[row+1+z]=a.gradient[z+2]-b.gradient[z+2];
      residual[row+4]=residual[row+4]-composition[0][z]*a.gradient[z+2]+composition[k][z]*b.gradient[z+2];
    }
  }
  constexpr std::size_t vol=5*(N-1);residual[vol]=Jet(-1);
  for(std::size_t z=0;z<3;++z)residual[vol+z+1]=Jet(-1);
  for(std::size_t k=0;k<N;++k){
    residual[vol]=residual[vol]+weight[k]*detail::exp(bulk[1]-variable[4*k]);
    for(std::size_t z=0;z<3;++z)residual[vol+z+1]=residual[vol+z+1]+weight[k]*composition[k][z]/bulk[z+2];
  }
  // Known first responses leave the second-order constraint residual. Solve
  // the same Newton matrix for all 15 quadratic coefficients simultaneously.
  std::array<unsigned,15> coefficients{};unsigned count=0;
  for(unsigned i=0;i<Jet::count;++i)if(Jet::layout.degree[i]==2)coefficients[count++]=i;
  std::array<Vector<15>,size> rhs{};
  for(std::size_t i=0;i<size;++i)for(std::size_t j=0;j<15;++j)rhs[i][j]=-residual[i].c[coefficients[j]];
  const auto correction=linear(s.jacobian,rhs);
  for(std::size_t i=0;i<size;++i)for(std::size_t j=0;j<15;++j)variable[i].c[coefficients[j]]+=correction[i][j];
  evaluate();std::array<Jet,5> gradient;
  for(std::size_t k=0;k<N;++k)gradient[0]=gradient[0]+weight[k]*phases[k].gradient[0];
  gradient[1]=detail::exp(variable[0]-bulk[1])*phases[0].gradient[1];
  for(std::size_t z=0;z<3;++z)gradient[z+2]=phases[0].gradient[z+2];
  Third result{};
  for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)for(std::size_t k=0;k<5;++k){
    Jet::Powers powers{};++powers[j];++powers[k];
    result[i][j][k]=gradient[i].derivative(powers)/scale[j]/scale[k];
  }
  for(std::size_t i=0;i<5;++i)for(std::size_t j=i;j<5;++j)for(std::size_t k=j;k<5;++k){
    const double value=(result[i][j][k]+result[j][i][k]+result[k][i][j])/3;
    if(!std::isfinite(value))throw std::domain_error("phase material: nonfinite analytic response");
    result[i][j][k]=result[i][k][j]=result[j][i][k]=result[j][k][i]=result[k][i][j]=result[k][j][i]=value;
  }
  return result;
}

PhaseMixtureMaterial material(const PhaseCoordinates&q,std::span<const PhaseState>guess,
    std::span<const PhaseEvaluator>evaluators,double step,bool shorten){
  if(!(std::isfinite(step)&&step>0))throw std::invalid_argument("phase material: invalid difference step");
  PhaseMixtureMaterial result;result.equilibrium=equilibrate_phases(q,guess,evaluators);
  std::array<std::array<PhaseCoordinates,5>,2> third{};
  const auto exact=guess.size()==1
      ? (result.equilibrium.potential.has_third
          ? std::optional<Third>(result.equilibrium.potential.third) : std::nullopt)
      : guess.size()==2?analytic_material<2>(q,result.equilibrium,evaluators):
                        analytic_material<3>(q,result.equilibrium,evaluators);
  if(exact){
    result.equilibrium.potential.third=*exact;result.equilibrium.potential.has_third=true;
    third[0]=(*exact)[0];third[1]=(*exact)[1];result.analytic_third=true;
    result.jets=material_jets(result.equilibrium.potential,third);return result;
  }
  for(std::size_t k=0;k<2;++k){
    double h=step;
    if(shorten)for(std::size_t j=0;j<result.equilibrium.phases.size();++j){
      const double dw=std::abs(result.equilibrium.fraction_response[j][k]);
      if(dw>0)h=std::min(h,.02*result.equilibrium.phases[j].mass_fraction/dw);
      for(const auto&response:result.equilibrium.state_response[j])
        if(std::abs(response[k])>0)h=std::min(h,.02/std::abs(response[k]));
    }
    auto plus=q,minus=q;plus[k]+=h;minus[k]-=h;
    if(plus[k]==q[k]||minus[k]==q[k])throw std::domain_error("phase material: unresolved difference step");
    auto predict=[&](double delta){
      auto states=result.equilibrium.phases;
      for(std::size_t j=0;j<states.size();++j){
        states[j].log_density+=result.equilibrium.state_response[j][0][k]*delta;
        for(std::size_t c=0;c<3;++c)
          states[j].composition[c]*=std::exp(result.equilibrium.state_response[j][c+1][k]*delta);
        states[j].mass_fraction+=result.equilibrium.fraction_response[j][k]*delta;
      }
      return states;
    };
    const auto a=equilibrate_phases(plus,predict(h),evaluators);
    const auto b=equilibrate_phases(minus,predict(-h),evaluators);
    result.log_steps[k]=h;
    for(std::size_t i=0;i<5;++i)for(std::size_t j=0;j<5;++j)
      third[k][i][j]=(a.potential.hessian[i][j]-b.potential.hessian[i][j])/(2*h);
  }
  result.jets=material_jets(result.equilibrium.potential,third);return result;
}
} // namespace

PhaseMixture equilibrate_phases(const PhaseCoordinates&q,std::span<const PhaseState>guess,
                                std::span<const PhaseEvaluator>evaluators){
  for(double x:q)if(!std::isfinite(x))throw std::invalid_argument("phase equilibrium: invalid state");
  if(!(q[2]>0&&q[3]>0&&q[4]>0&&q[2]+q[3]+q[4]<1))
    throw std::domain_error("phase equilibrium: positive composition required");
  if(guess.size()!=evaluators.size())throw std::invalid_argument("phase equilibrium: evaluator count");
  if(guess.size()==1){
    PhaseMixture result;result.potential=evaluators[0](q);check_potential(result.potential);
    result.phases.push_back({q[1],{q[2],q[3],q[4]},1});
    result.fraction_response.resize(1);result.state_response.resize(1);
    result.state_response[0][0][1]=1;
    for(std::size_t k=0;k<3;++k)result.state_response[0][k+1][k+2]=1/q[k+2];
    return result;
  }
  if(guess.size()==2)return solve<2>(q,guess,evaluators);
  if(guess.size()==3)return solve<3>(q,guess,evaluators);
  throw std::invalid_argument("phase equilibrium: one to three phases required");
}
PhaseEquilibrium equilibrate_phases(const PhaseCoordinates&q,const PhaseSplit&guess,
    const PhaseEvaluator&first,const PhaseEvaluator&second){
  const std::array<PhaseEvaluator,2> evaluators{first,second};
  return pair_result(equilibrate_phases(q,pair_guess(guess),evaluators));
}
PhaseMaterial phase_equilibrium_material(const PhaseCoordinates&q,const PhaseSplit&guess,
    const PhaseEvaluator&first,const PhaseEvaluator&second,double step){
  const std::array<PhaseEvaluator,2> evaluators{first,second};
  const auto r=material(q,pair_guess(guess),evaluators,step,false);
  return {pair_result(r.equilibrium),r.jets};
}
PhaseMixtureMaterial phase_equilibrium_material(const PhaseCoordinates&q,std::span<const PhaseState>guess,
    std::span<const PhaseEvaluator>evaluators,double step){return material(q,guess,evaluators,step,true);}
} // namespace ember

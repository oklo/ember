#include "ember/material_transport.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ember {
namespace {
using V = MaterialVector;
using M = MaterialMatrix;
void size_check(std::size_t n) {
  if(n<1 || n>3)throw std::invalid_argument("material transport: components must be 1, 2 or 3");
}
M identity(std::size_t n) {M a{};for(std::size_t i=0;i<n;++i)a[i][i]=1;return a;}
M transpose(const M& a,std::size_t n) {
  M b{};for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)b[i][j]=a[j][i];return b;
}
M product(const M& a,const M& b,std::size_t n) {
  M c{};for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
    for(std::size_t k=0;k<n;++k)c[i][j]+=a[i][k]*b[k][j];return c;
}
V product(const M& a,const V& b,std::size_t n) {
  V c{};for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)c[i]+=a[i][j]*b[j];return c;
}
M symmetric(M a,std::size_t n) {
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<i;++j)
    a[i][j]=a[j][i]=.5*a[i][j]+.5*a[j][i];return a;
}
void finite(const V& a,std::size_t n) {
  for(std::size_t i=0;i<n;++i)if(!std::isfinite(a[i]))
    throw std::domain_error("material transport: non-finite vector");
}
void finite(const M& a,std::size_t n) {for(std::size_t i=0;i<n;++i)finite(a[i],n);}
void symmetric_check(const M& a,std::size_t n) {
  finite(a,n);
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<i;++j)
    if(a[i][j]!=a[j][i])throw std::domain_error("material transport: asymmetric matrix");
}
bool zero(const M& a,std::size_t n) {
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)if(a[i][j]!=0)return false;
  return true;
}
M cholesky(const M& a,std::size_t n) {
  finite(a,n);M l{};
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j) {
    double v=a[i][j];for(std::size_t k=0;k<j;++k)v-=l[i][k]*l[j][k];
    if(i==j) {
      if(!(v>0) || !std::isfinite(v))throw std::domain_error("material transport: matrix is not resolved as positive");
      l[i][j]=std::sqrt(v);
    } else l[i][j]=v/l[j][j];
  }
  return l;
}
M lower_inverse(const M& l,std::size_t n) {
  M a{};for(std::size_t j=0;j<n;++j)for(std::size_t i=j;i<n;++i) {
    double v=i==j?1:0;for(std::size_t k=j;k<i;++k)v-=l[i][k]*a[k][j];a[i][j]=v/l[i][i];
  }return a;
}

// Symmetric Jacobi rotations, scaled before diagonalization. This small kernel
// is independently compared with dense solves and noncommuting stiff blocks.
M fraction_of_positive(const M& input,std::size_t n) {
  double scale=0;for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
    scale=std::max(scale,std::abs(input[i][j]));
  if(scale==0)return {};
  if(!std::isfinite(scale))throw std::domain_error("material transport: non-finite whitened conductance");
  M a=input,q=identity(n);
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)a[i][j]/=scale;
  bool converged=n==1;
  for(unsigned iteration=0;iteration<64 && n>1;++iteration) {
    std::size_t p=0,r=1;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j)
      if(std::abs(a[i][j])>std::abs(a[p][r])){p=i;r=j;}
    if(std::abs(a[p][r])<=4*std::numeric_limits<double>::epsilon()) {converged=true;break;}
    const double tau=(a[r][r]-a[p][p])/(2*a[p][r]);
    const double t=std::copysign(1.0,tau)/(std::abs(tau)+std::hypot(1.0,tau));
    const double c=1/std::sqrt(1+t*t),s=t*c,off=a[p][r];
    a[p][p]-=t*off;a[r][r]+=t*off;a[p][r]=a[r][p]=0;
    for(std::size_t i=0;i<n;++i) {
      if(i!=p && i!=r) {
        const double ip=a[i][p],ir=a[i][r];
        a[i][p]=a[p][i]=c*ip-s*ir;a[i][r]=a[r][i]=s*ip+c*ir;
      }
      const double ip=q[i][p],ir=q[i][r];q[i][p]=c*ip-s*ir;q[i][r]=s*ip+c*ir;
    }
  }
  if(!converged)throw std::domain_error("material transport: symmetric eigensolve did not converge");
  M f{};
  for(std::size_t k=0;k<n;++k) {
    const double value=a[k][k]*scale;
    if(!(value>0) || !std::isfinite(value))
      throw std::domain_error("material transport: conductance eigenvalue is not resolved as positive");
    const double v=value/(1+value);
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)f[i][j]+=q[i][k]*v*q[j][k];
  }
  return symmetric(f,n);
}
} // namespace

MaterialMatrix material_positive_inverse(const M& a,std::size_t n) {
  size_check(n);finite(a,n);
  V scale{};M normalized{};
  for(std::size_t i=0;i<n;++i) {
    if(!(a[i][i]>0))throw std::domain_error("material transport: nonpositive diagonal");
    scale[i]=std::sqrt(a[i][i]);
  }
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)
    normalized[i][j]=(.5*a[i][j]+.5*a[j][i])/scale[i]/scale[j];
  const auto inv_l=lower_inverse(cholesky(normalized,n),n);
  auto answer=product(transpose(inv_l,n),inv_l,n);
  for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)answer[i][j]/=scale[i]*scale[j];
  finite(answer,n);return symmetric(answer,n);
}

std::vector<V> solve_material_chain(std::span<const M> capacities,
    std::span<const M> conductances,std::span<const V> rhs,std::size_t n) {
  size_check(n);const auto count=rhs.size();
  if(count==0 || capacities.size()!=count || conductances.size()!=count-1)
    throw std::invalid_argument("material transport: mismatched block chain sizes");
  std::vector<M> effective(capacities.begin(),capacities.end()),transfers(count-1),inverses(count-1);
  std::vector<V> right(rhs.begin(),rhs.end()),answer(count);
  for(const auto& a:capacities){symmetric_check(a,n);cholesky(a,n);}
  for(const auto& k:conductances){symmetric_check(k,n);if(!zero(k,n))cholesky(k,n);}
  for(const auto& b:rhs)finite(b,n);
  for(std::size_t i=0;i+1<count;++i) {
    const auto e=symmetric(effective[i],n),l=cholesky(e,n),inv_l=lower_inverse(l,n);
    const auto& k=conductances[i];
    const auto c=symmetric(product(product(inv_l,k,n),transpose(inv_l,n),n),n);
    const auto f=fraction_of_positive(c,n);
    const auto lf=product(l,f,n),parallel=product(lf,transpose(l,n),n);
    const auto transfer=product(lf,inv_l,n);
    const auto update=product(transfer,right[i],n);
    M sum=e;
    for(std::size_t j=0;j<n;++j) {
      right[i+1][j]+=update[j];
      for(std::size_t v=0;v<n;++v){effective[i+1][j][v]+=parallel[j][v];sum[j][v]+=k[j][v];}
    }
    transfers[i]=transpose(transfer,n);inverses[i]=material_positive_inverse(sum,n);
  }
  answer.back()=product(material_positive_inverse(effective.back(),n),right.back(),n);
  for(std::size_t i=count-1;i-->0;) {
    answer[i]=product(inverses[i],right[i],n);
    const auto next=product(transfers[i],answer[i+1],n);
    for(std::size_t j=0;j<n;++j)answer[i][j]+=next[j];
  }
  for(const auto& a:answer)finite(a,n);
  return answer;
}

MaterialTransportResult transport_material(std::span<const V> old_primitives,
    std::span<const V> initial_guess,std::span<const double> mass,double dt,
    const MaterialThermodynamics& thermodynamics,const MaterialConductance& conductance,
    const MaterialTransportOptions& options) {
  const auto n=options.components,count=mass.size();size_check(n);
  if(count==0 || old_primitives.size()!=count || initial_guess.size()!=count || !thermodynamics || !conductance
      || !(dt>0) || !std::isfinite(dt) || !(options.tolerance>0) || !std::isfinite(options.tolerance)
      || !(options.composition_sum_limit>0) || options.composition_sum_limit>1
      || !options.max_iterations || !options.max_backtracks)
    throw std::invalid_argument("material transport: invalid sizes, callbacks, time or options");
  V scale{};scale.fill(1);std::vector<MaterialPoint> old(count);double old_entropy=0;
  auto check_primitive=[&](const V& q,bool interior) {
    finite(q,n);double sum=0;
    for(std::size_t j=0;j+1<n;++j) {
      if(q[j]<0 || (interior && q[j]==0))throw std::domain_error("material transport: composition leaves interior");
      sum+=q[j];
    }
    if(sum>options.composition_sum_limit || (interior && sum==options.composition_sum_limit))
      throw std::domain_error("material transport: reference species leaves interior");
  };
  for(std::size_t i=0;i<count;++i) {
    if(!(mass[i]>0) || !std::isfinite(mass[i]))throw std::invalid_argument("material transport: invalid cell mass");
    check_primitive(old_primitives[i],false);
    old[i]=thermodynamics(i,old_primitives[i],false);finite(old[i].conserved,n);
    if(!std::isfinite(old[i].entropy))throw std::domain_error("material transport: non-finite old entropy");
    old_entropy+=mass[i]*old[i].entropy;
    for(std::size_t j=0;j<n;++j)scale[j]=std::max(scale[j],std::abs(old[i].conserved[j]));
  }
  struct Evaluation {
    std::vector<MaterialPoint> points;
    std::vector<M> k;
    std::vector<V> flux,residual;
    double norm{};
  };
  auto evaluate=[&](std::span<const V> q) {
    Evaluation a;a.points.reserve(count);a.flux.resize(count-1);a.residual.resize(count);
    for(std::size_t i=0;i<count;++i) {
      check_primitive(q[i],true);auto p=thermodynamics(i,q[i],true);
      finite(p.conserved,n);finite(p.potential,n);finite(p.primitive_from_conserved,n);
      symmetric_check(p.capacity,n);cholesky(p.capacity,n);
      if(!std::isfinite(p.entropy))throw std::domain_error("material transport: non-finite entropy");
      for(std::size_t j=0;j<n;++j)a.residual[i][j]=mass[i]*(p.conserved[j]-old[i].conserved[j]);
      a.points.push_back(p);
    }
    a.k=conductance(q);
    if(a.k.size()!=count-1)throw std::invalid_argument("material transport: wrong number of faces");
    for(std::size_t i=0;i+1<count;++i) {
      symmetric_check(a.k[i],n);if(!zero(a.k[i],n))cholesky(a.k[i],n);
      V difference{};for(std::size_t j=0;j<n;++j)difference[j]=a.points[i].potential[j]-a.points[i+1].potential[j];
      a.flux[i]=product(a.k[i],difference,n);
      for(std::size_t j=0;j<n;++j) {
        a.residual[i][j]+=dt*a.flux[i][j];a.residual[i+1][j]-=dt*a.flux[i][j];
      }
    }
    for(std::size_t i=0;i<count;++i) {
      finite(a.residual[i],n);
      for(std::size_t j=0;j<n;++j)a.norm=std::max(a.norm,std::abs(a.residual[i][j])/mass[i]/scale[j]);
    }
    if(!std::isfinite(a.norm))throw std::domain_error("material transport: non-finite residual norm");
    return a;
  };
  MaterialTransportResult result;result.primitives.assign(initial_guess.begin(),initial_guess.end());
  auto state=evaluate(result.primitives);
  for(std::size_t iteration=0;iteration<options.max_iterations;++iteration) {
    result.residual_history.push_back(state.norm);
    if(state.norm<=options.tolerance) {
      result.iterations=iteration;result.residual=state.norm;result.face_flux=state.flux;
      double entropy=0;
      for(std::size_t i=0;i<count;++i) {
        result.conserved.push_back(state.points[i].conserved);entropy+=mass[i]*state.points[i].entropy;
        for(std::size_t j=0;j<n;++j)result.integrated_conservation_error[j]+=mass[i]*(state.points[i].conserved[j]-old[i].conserved[j]);
      }
      result.entropy_change=entropy-old_entropy;
      for(std::size_t i=0;i+1<count;++i)for(std::size_t j=0;j<n;++j)
        result.backward_euler_entropy_bound-=dt*state.flux[i][j]*(state.points[i+1].potential[j]-state.points[i].potential[j]);
      return result;
    }
    std::vector<M> capacity(count),k=state.k;std::vector<V> rhs(count),delta(count);
    for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<n;++j) {
      rhs[i][j]=-state.residual[i][j];
      for(std::size_t v=0;v<n;++v)capacity[i][j][v]=mass[i]*state.points[i].capacity[j][v];
    }
    for(auto& a:k)for(std::size_t j=0;j<n;++j)for(std::size_t v=0;v<n;++v)a[j][v]*=dt;
    const auto dw=solve_material_chain(capacity,k,rhs,n);
    for(std::size_t i=0;i<count;++i)delta[i]=product(state.points[i].primitive_from_conserved,product(state.points[i].capacity,dw[i],n),n);
    bool accepted=false;double damping=1;
    for(std::size_t backtrack=0;backtrack<options.max_backtracks;++backtrack,damping*=.5) {
      auto candidate=result.primitives;
      for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<n;++j)candidate[i][j]+=damping*delta[i][j];
      try {
        auto trial=evaluate(candidate);
        if(trial.norm<state.norm*(1-1e-4*damping) || trial.norm<=options.tolerance) {
          result.primitives=std::move(candidate);state=std::move(trial);accepted=true;break;
        }
      } catch(const std::domain_error&) { /* reject the trial without changing the accepted state */ }
    }
    if(!accepted)throw std::runtime_error("material transport: line search failed at residual "+std::to_string(state.norm));
  }
  throw std::runtime_error("material transport: iteration limit at residual "+std::to_string(state.norm));
}
} // namespace ember

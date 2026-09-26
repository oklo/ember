#include "ember/cn_burning.hpp"
#include "ember/cn_transport.hpp"
#include "cn_source.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
using Vec=detail::CNVector;
using Mat=detail::CNMatrix;
constexpr Vec scale{1,1,12,13};
Vec solve(Mat a,Vec b) {
  for(std::size_t k=0;k<4;++k) {
    std::size_t pivot=k;
    for(std::size_t j=k+1;j<4;++j)if(std::abs(a[j][k])>std::abs(a[pivot][k]))pivot=j;
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    if(!std::isfinite(a[k][k]) || a[k][k]==0)
      throw std::runtime_error("CN transport: singular chemical Jacobian");
    for(std::size_t i=k+1;i<4;++i) {
      const double f=a[i][k]/a[k][k];
      for(std::size_t j=k;j<4;++j)a[i][j]-=f*a[k][j];
      b[i]-=f*b[k];
    }
  }
  Vec x{};
  for(std::size_t i=4;i-->0;) {
    double rhs=b[i];for(std::size_t j=i+1;j<4;++j)rhs-=a[i][j]*x[j];
    x[i]=rhs/a[i][i];
    if(!std::isfinite(x[i]))throw std::runtime_error("CN transport: nonfinite solution");
  }
  return x;
}
double matrix_scale(const Mat& a,double g) {
  for(const auto& row:a)for(double x:row)g=std::max(g,std::abs(x));
  return g;
}
Mat shifted(const Mat& a,double g,double s) {
  Mat out=a;for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)
    out[i][j]=a[i][j]/s+(i==j?g/s:0);
  return out;
}
Vec multiply(const Mat& a,const Vec& b) {
  Vec out{};for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)out[i]+=a[i][j]*b[j];return out;
}
// W = g (E+gI)^-1. Form the next Schur term as W E: this retains
// the finite mass/reaction matrix even when g exceeds it by many decades.
Mat transfer(const Mat& E,double g) {
  const double s=matrix_scale(E,g);const auto a=shifted(E,g,s);Mat w{};
  for(std::size_t col=0;col<4;++col) {
    Vec rhs{};rhs[col]=g/s;const auto x=solve(a,rhs);
    for(std::size_t row=0;row<4;++row)w[row][col]=x[row];
  }
  return w;
}
}

std::vector<Composition> burn_cn_and_transport(const Model& thermal,const Model& previous,
    const PPCNNetwork& nuclear,const MixingRegions& regions,const std::vector<double>& D,
    double dt,double tolerance) {
  const auto weights=nodal_mass_weights(thermal);
  if(thermal.m!=previous.m || thermal.comp.size()!=thermal.size() || previous.comp.size()!=thermal.size()
      || D.size()+1!=thermal.size() || !std::isfinite(dt+tolerance) || dt<=0 || tolerance<=0)
    throw std::invalid_argument("CN transport: inconsistent mesh, step or diffusivities");
  for(double d:D)if(!std::isfinite(d) || d<0)throw std::invalid_argument("CN transport: invalid diffusivity");
  if(std::all_of(D.begin(),D.end(),[](double d){return d==0;}))
    return burn_and_mix(thermal,previous,nuclear,regions,dt,tolerance);
  for(const auto& c:previous.comp) {
    if(c[Species::H2]!=0)
      throw std::invalid_argument("CN transport: initial D requires the physical seven-species solve");
    if(!c.cn_molality)throw std::invalid_argument("CN transport: missing catalyst inventory");
    (void)cn_physical_ledger(c,*c.cn_molality);
    for(std::size_t j=3;j<NSPEC;++j)if(c.X[j]!=previous.comp.front().X[j])
      throw std::invalid_argument("CN transport: uniform fixed material metal lookup required");
  }
  const auto& first=*previous.comp.front().cn_molality;
  const double first_total=first[0]+first[1]+first[2];
  const bool variable_total=std::any_of(previous.comp.begin(),previous.comp.end(),[&](const Composition& c) {
    const auto& y=*c.cn_molality;
    return std::abs(y[0]+y[1]+y[2]-first_total)>64*std::numeric_limits<double>::epsilon()*first_total;
  });
  if(variable_total) {
    std::vector<double> mixing(D.size());
    for(std::size_t i=0;i<D.size();++i) {
      const double r=.5*(thermal.r(i)+thermal.r(i+1)),rho=.5*(thermal.rho(i)+thermal.rho(i+1));
      const double area=4*M_PI*r*r*rho;mixing[i]=area*area*D[i]/(thermal.m[i+1]-thermal.m[i]);
    }
    CNSpeciesFlux flux=[](std::size_t,const Composition&,const Composition&,bool){return CNSpeciesFaceResponse{};};
    SpeciesTransportOptions options;options.abundance_tolerance=tolerance;
    return burn_cn_and_diffuse(thermal,previous,nuclear,regions,flux,dt,options,mixing).composition;
  }
  std::vector<Vec> old;std::vector<Composition> lookup;std::vector<double> mass,total,g;
  std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())throw std::invalid_argument("CN transport: invalid partition");
    next=end;Vec u{};double w=0,n=0;
    for(std::size_t i=begin;i<end;++i) {
      const auto& c=previous.comp[i];const auto& y=*c.cn_molality;
      w+=weights[i];n+=weights[i]*(y[0]+y[1]+y[2]);
      u[0]+=weights[i]*c.X[0];u[1]+=weights[i]*c.X[1];
      u[2]+=weights[i]*12*y[0];u[3]+=weights[i]*13*y[1];
    }
    for(double& v:u)v/=w;old.push_back(u);mass.push_back(w/thermal.M);
    total.push_back(n/w);lookup.push_back(previous.comp[begin]);
    if(end<thermal.size()) {
      const double r=.5*(thermal.r(end-1)+thermal.r(end)),rho=.5*(thermal.rho(end-1)+thermal.rho(end));
      const double a=4*M_PI*r*r*rho;
      const double coupling=dt/thermal.M*a*a*D[end-1]/(thermal.m[end]-thermal.m[end-1]);
      if(!std::isfinite(coupling))throw std::domain_error("CN transport: unrepresentable diffusion coupling");
      g.push_back(coupling);
    }
  }
  if(next!=thermal.size())throw std::invalid_argument("CN transport: incomplete partition");
  // The total catalyst number is spatially uniform for this fixed-Z model.
  // All three catalysts share D; their summed flux therefore vanishes.
  const auto initial=initial_gs98_cn(previous.comp.front());
  const auto& first_cn=*previous.comp.front().cn_molality;
  const double uniform_total=first_cn[0]+first_cn[1]+first_cn[2];
  for(const auto& c:previous.comp) {
    const auto& y=*c.cn_molality;
    if(std::abs(y[0]+y[1]+y[2]-uniform_total)>64*std::numeric_limits<double>::epsilon()*uniform_total)
      throw std::invalid_argument("CN transport: catalyst-number gradients require a larger network");
  }
  std::fill(total.begin(),total.end(),uniform_total);
  auto composition=[&](std::size_t i,const Vec& u) {
    auto c=lookup[i];c.X[0]=u[0];c.X[1]=u[1];c.X[2]=1-c.Z()-u[0]-u[1];
    c.cn_molality=CNAbundances{u[2]/12,u[3]/13,total[i]-u[2]/12-u[3]/13};return c;
  };
  auto positive=[&](std::size_t i,const Vec& u) {
    for(double v:u)if(!std::isfinite(v) || v<0)return false;
    const auto c=composition(i,u);const auto& y=*c.cn_molality;
    if(c.X[2]<0 || y[2]<0)return false;
    const double extra=12*(y[0]-initial[0])+13*y[1]+14*(y[2]-initial[2]);
    return c.X[2]-extra>=0;
  };
  auto current=old;const auto count=regions.size();
  for(int iteration=0;iteration<80;++iteration) {
    std::vector<Mat> E(count);std::vector<Vec> rhs(count),answer(count);
    for(std::size_t i=0;i<count;++i) {
      for(std::size_t k=0;k<4;++k){E[i][k][k]=mass[i];rhs[i][k]=mass[i]*old[i][k];}
      const auto c=composition(i,current[i]);
      for(std::size_t j=regions[i].first;j<regions[i].second;++j) {
        const auto rates=detail::cn_reduced_source(thermal.T(j),thermal.rho(j),c,*c.cn_molality,nuclear.pp(),nuclear.cn());
        const double factor=dt*weights[j]/thermal.M;
        for(std::size_t row=0;row<4;++row) {
          rhs[i][row]+=factor*scale[row]*rates.source[row];
          for(std::size_t col=0;col<4;++col) {
            const double J=scale[row]/scale[col]*rates.jacobian[row][col];
            E[i][row][col]-=factor*J;rhs[i][row]-=factor*J*current[i][col];
          }
        }
      }
      if(i && g[i-1]>0) {
        const auto W=transfer(E[i-1],g[i-1]);const auto b=multiply(W,rhs[i-1]);
        for(std::size_t row=0;row<4;++row) {
          rhs[i][row]+=b[row];
          for(std::size_t col=0;col<4;++col)for(std::size_t k=0;k<4;++k)
            E[i][row][col]+=W[row][k]*E[i-1][k][col];
        }
      }
    }
    for(std::size_t i=count;i-->0;) {
      const double coupling=i+1<count?g[i]:0,s=matrix_scale(E[i],coupling);
      Vec b=rhs[i];for(std::size_t k=0;k<4;++k) {
        b[k]/=s;if(i+1<count)b[k]+=coupling/s*answer[i+1][k];
      }
      answer[i]=solve(shifted(E[i],coupling,s),b);
    }
    double damping=1;bool accepted=false;
    std::vector<Vec> candidate(count);
    for(int trial=0;trial<50;++trial,damping*=.5) {
      accepted=true;
      for(std::size_t i=0;i<count;++i) {
        for(std::size_t k=0;k<4;++k)candidate[i][k]=current[i][k]+damping*(answer[i][k]-current[i][k]);
        if(!positive(i,candidate[i])){accepted=false;break;}
      }
      if(accepted)break;
    }
    if(!accepted)throw std::runtime_error("CN transport: no positive Newton step");
    double change=0;
    for(std::size_t i=0;i<count;++i)for(std::size_t k=0;k<4;++k)
      change=std::max(change,std::abs(answer[i][k]-current[i][k]));
    current=std::move(candidate);
    if(change<tolerance && damping==1) {
      Vec balance{};
      auto result=previous.comp;
      for(std::size_t i=0;i<count;++i) {
        const auto c=composition(i,current[i]);(void)cn_physical_ledger(c,*c.cn_molality);
        for(std::size_t k=0;k<4;++k)balance[k]+=mass[i]*(current[i][k]-old[i][k]);
        for(std::size_t j=regions[i].first;j<regions[i].second;++j) {
          result[j]=c;
          const auto rates=detail::cn_reduced_source(thermal.T(j),thermal.rho(j),c,*c.cn_molality,nuclear.pp(),nuclear.cn());
          for(std::size_t k=0;k<4;++k)balance[k]-=dt*weights[j]/thermal.M*scale[k]*rates.source[k];
        }
      }
      for(double b:balance)if(std::abs(b)>10*tolerance)
        throw std::runtime_error("CN transport: integrated reaction/transport imbalance");
      return result;
    }
  }
  throw std::runtime_error("CN transport: implicit chemical solve did not converge");
}
} // namespace ember

#include "ember/deuterium_burning.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace ember {
namespace {
constexpr std::array<std::size_t,3> active{0,1,static_cast<std::size_t>(Species::H2)};
using Vec=std::array<double,3>;
using Mat=std::array<Vec,3>;
double norm(const Vec& f){return std::max({std::abs(f[0]),std::abs(f[1]),std::abs(f[2])});}
Vec solve(Mat a,Vec b) {
  for(std::size_t k=0;k<3;++k) {
    std::size_t pivot=k;
    for(std::size_t i=k+1;i<3;++i)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    if(a[k][k]==0 || !std::isfinite(a[k][k]))throw std::runtime_error("D burning: singular Jacobian");
    for(std::size_t i=k+1;i<3;++i) {
      const double q=a[i][k]/a[k][k];
      for(std::size_t j=k;j<3;++j)a[i][j]-=q*a[k][j];
      b[i]-=q*b[k];
    }
  }
  Vec x{};
  for(std::size_t i=3;i-->0;) {
    double v=b[i];for(std::size_t j=i+1;j<3;++j)v-=a[i][j]*x[j];x[i]=v/a[i][i];
  }
  return x;
}
Mat inverse(const Mat& a) {
  Mat result{};
  for(std::size_t j=0;j<3;++j) {
    Vec unit{};unit[j]=1;const auto column=solve(a,unit);
    for(std::size_t i=0;i<3;++i)result[i][j]=column[i];
  }
  return result;
}
Vec multiply(const Mat& a,const Vec& b) {
  Vec result{};
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)result[i]+=a[i][j]*b[j];
  return result;
}
Mat multiply(const Mat& a,const Mat& b) {
  Mat result{};
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
    for(std::size_t k=0;k<3;++k)result[i][j]+=a[i][k]*b[k][j];
  return result;
}
double scale(const Mat& a,double g) {
  for(const auto& row:a)for(double v:row)g=std::max(g,std::abs(v));
  return g;
}
// g (E+gI)^-1. The Schur update W E preserves the reaction/mass matrix
// when g greatly exceeds it; subtracting two diffusion matrices would not.
Mat transfer(const Mat& E,double g) {
  const double s=scale(E,g);Mat a=E;
  for(std::size_t i=0;i<3;++i) {
    for(double& v:a[i])v/=s;
    a[i][i]+=g/s;
  }
  auto result=inverse(a);
  for(auto& row:result)for(double& v:row)v*=g/s;
  return result;
}
Composition region_mean(const Model& model,const std::vector<double>& weights,
                        std::size_t begin,std::size_t end,double& mass) {
  const auto& reference=model.comp[begin];auto result=reference;
  std::array<double,NSPEC> difference{};mass=0;
  for(std::size_t i=begin;i<end;++i) {
    if(model.comp[i].metal_inventory!=reference.metal_inventory)
      throw std::invalid_argument("D burning: inconsistent metal convention");
    mass+=weights[i];
    for(std::size_t j=0;j<NSPEC;++j)
      difference[j]+=weights[i]*(model.comp[i].X[j]-reference.X[j]);
  }
  // Center the average on an actual composition. Uniform abundances then
  // survive mixing exactly, including when the nuclear changes are tiny.
  for(std::size_t j=0;j<NSPEC;++j)result.X[j]+=difference[j]/mass;
  return result;
}
}
std::vector<Composition> burn_deuterium_and_mix(const Model& thermal,const Model& previous,
    const PPDeuterium& nuclear,const MixingRegions& regions,double dt,double tolerance) {
  if(thermal.m!=previous.m || thermal.size()!=previous.size() || previous.comp.size()!=thermal.size()
      || thermal.comp.size()!=thermal.size() || !std::isfinite(dt+tolerance) || dt<=0 || tolerance<=0)
    throw std::invalid_argument("D burning: inconsistent mesh or step");
  const auto weights=nodal_mass_weights(thermal);
  for(const auto& c:previous.comp) {
    if(c.basis!=AbundanceBasis::baryon_mass || std::abs(c.sum()-1)>1e-10 || c.cn_molality)
      throw std::invalid_argument("D burning: normalized baryon inventory without CN ledger required");
    for(double x:c.X)if(!std::isfinite(x) || x<0)throw std::invalid_argument("D burning: invalid isotope abundance");
  }
  auto result=previous.comp;std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())throw std::invalid_argument("D burning: invalid mixing partition");
    next=end;
    double mass=0;const auto old=region_mean(previous,weights,begin,end,mass);
    Vec prior{old.X[active[0]],old.X[active[1]],old.X[active[2]]},u=prior;
    auto composition=[&](const Vec& v) {
      auto c=old;for(std::size_t k=0;k<3;++k)c.X[active[k]]=v[k];
      c.X[2]=old.X[2]-((v[0]-prior[0])+(v[1]-prior[1])+(v[2]-prior[2]));return c;
    };
    auto positive=[&](const Vec& v) {
      for(double x:v)if(!std::isfinite(x) || x<0)return false;
      return old.X[2]-((v[0]-prior[0])+(v[1]-prior[1])+(v[2]-prior[2]))>=0;
    };
    auto equation=[&](const Vec& v,Mat* j) {
      Vec f{};for(std::size_t k=0;k<3;++k)f[k]=v[k]-prior[k];
      if(j){*j={};for(std::size_t k=0;k<3;++k)(*j)[k][k]=1;}
      const auto c=composition(v);
      for(std::size_t i=begin;i<end;++i) {
        const auto response=nuclear.composition_response(thermal.T(i),thermal.rho(i),c);
        const double factor=dt*weights[i]/mass;
        for(std::size_t row=0;row<3;++row) {
          f[row]-=factor*response.state.dXdt[active[row]];
          if(j)for(std::size_t col=0;col<3;++col)
            (*j)[row][col]-=factor*(response.d_dXdt_dX[active[row]][active[col]]
                                  -response.d_dXdt_dX[active[row]][2]);
        }
      }
      return f;
    };
    // Resolve initial-D destruction relative to its remaining reservoir,
    // even when its absolute abundance is below the other species' tolerance.
    // This changes Newton convergence, not the reaction equations or rates.
    const double dscale=std::min(tolerance,std::max(1e-12*prior[2],
        16*std::numeric_limits<double>::denorm_min()));
    const auto residual_norm=[&](const Vec& f){return std::max({
        std::abs(f[0])/tolerance,std::abs(f[1])/tolerance,std::abs(f[2])/dscale});};
    bool converged=false;
    for(int iteration=0;iteration<60;++iteration) {
      Mat j;const auto f=equation(u,&j);const double error=residual_norm(f);
      if(!std::isfinite(error))throw std::runtime_error("D burning: nonfinite residual");
      if(error<1){converged=true;break;}
      const auto correction=solve(j,{-f[0],-f[1],-f[2]});bool accepted=false;
      for(double damping=1;damping>1e-14;damping*=.5) {
        Vec trial=u;for(std::size_t k=0;k<3;++k)trial[k]+=damping*correction[k];
        if(positive(trial) && residual_norm(equation(trial,nullptr))<error){u=trial;accepted=true;break;}
      }
      if(!accepted)throw std::runtime_error("D burning: abundance line search failed");
    }
    if(!converged)throw std::runtime_error("D burning: abundance iteration limit");
    for(std::size_t i=begin;i<end;++i)result[i]=composition(u);
  }
  if(next!=thermal.size())throw std::invalid_argument("D burning: incomplete mixing partition");
  return result;
}

std::vector<Composition> burn_deuterium_and_transport(const Model& thermal,const Model& previous,
    const PPDeuterium& nuclear,const MixingRegions& regions,const std::vector<double>& D,
    double dt,double tolerance) {
  if(thermal.m!=previous.m || thermal.size()!=previous.size() || previous.comp.size()!=thermal.size()
      || thermal.comp.size()!=thermal.size() || D.size()+1!=thermal.size()
      || !std::isfinite(dt+tolerance) || dt<=0 || tolerance<=0)
    throw std::invalid_argument("D transport: inconsistent mesh, diffusivities or step");
  for(double v:D)if(!std::isfinite(v) || v<0)
    throw std::invalid_argument("D transport: invalid diffusivity");
  if(std::all_of(D.begin(),D.end(),[](double v){return v==0;}))
    return burn_deuterium_and_mix(thermal,previous,nuclear,regions,dt,tolerance);
  const auto weights=nodal_mass_weights(thermal);
  for(const auto& c:previous.comp) {
    if(c.basis!=AbundanceBasis::baryon_mass || std::abs(c.sum()-1)>1e-10 || c.cn_molality
        || c.cn_mass_convention!=CNMassConvention::fixed_metal_proxy
        || c.metal_inventory!=previous.comp.front().metal_inventory)
      throw std::invalid_argument("D transport: normalized baryon inventory with common fixed metals required");
    for(std::size_t j=0;j<NSPEC;++j) {
      if(!std::isfinite(c.X[j]) || c.X[j]<0)
        throw std::invalid_argument("D transport: invalid isotope abundance");
      if(is_metal_species(j) && c.X[j]!=previous.comp.front().X[j])
        throw std::invalid_argument("D transport: metal gradients need a larger transport system");
    }
  }
  std::vector<Composition> old;std::vector<double> mass,g;std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())
      throw std::invalid_argument("D transport: invalid mixing partition");
    next=end;double w=0;const auto c=region_mean(previous,weights,begin,end,w);
    old.push_back(c);mass.push_back(w/thermal.M);
    if(end<thermal.size()) {
      const double r=.5*(thermal.r(end-1)+thermal.r(end)),rho=.5*(thermal.rho(end-1)+thermal.rho(end));
      const double area_mass=4*std::acos(-1.)*r*r*rho;
      const double coupling=dt/thermal.M*area_mass*area_mass*D[end-1]/(thermal.m[end]-thermal.m[end-1]);
      if(!std::isfinite(coupling))throw std::domain_error("D transport: unrepresentable mixing coefficient");
      g.push_back(coupling);
    }
  }
  if(next!=thermal.size())throw std::invalid_argument("D transport: incomplete mixing partition");
  auto current=old;const auto count=regions.size();
  for(int iteration=0;iteration<60;++iteration) {
    std::vector<Mat> E(count);std::vector<Vec> rhs(count),answer(count);
    for(std::size_t i=0;i<count;++i) {
      for(std::size_t k=0;k<3;++k){E[i][k][k]=mass[i];rhs[i][k]=mass[i]*old[i].X[active[k]];}
      for(std::size_t j=regions[i].first;j<regions[i].second;++j) {
        const auto response=nuclear.composition_response(thermal.T(j),thermal.rho(j),current[i]);
        const double f=dt*weights[j]/thermal.M;
        for(std::size_t row=0;row<3;++row) {
          rhs[i][row]+=f*response.state.dXdt[active[row]];
          for(std::size_t col=0;col<3;++col) {
            const double J=response.d_dXdt_dX[active[row]][active[col]]-response.d_dXdt_dX[active[row]][2];
            E[i][row][col]-=f*J;rhs[i][row]-=f*J*current[i].X[active[col]];
          }
        }
      }
      if(i && g[i-1]>0) {
        const auto W=transfer(E[i-1],g[i-1]),addition=multiply(W,E[i-1]);
        const auto b=multiply(W,rhs[i-1]);
        for(std::size_t row=0;row<3;++row) {
          rhs[i][row]+=b[row];
          for(std::size_t col=0;col<3;++col)E[i][row][col]+=addition[row][col];
        }
      }
    }
    for(std::size_t i=count;i-->0;) {
      const double coupling=i+1<count?g[i]:0,s=scale(E[i],coupling);Mat a=E[i];Vec b=rhs[i];
      for(std::size_t k=0;k<3;++k) {
        for(double& v:a[k])v/=s;
        a[k][k]+=coupling/s;b[k]/=s;
        if(i+1<count)b[k]+=coupling/s*answer[i+1][k];
      }
      answer[i]=solve(a,b);
    }
    double damping=1;
    for(std::size_t i=0;i<count;++i) {
      double helium_change=0;
      for(std::size_t k=0;k<3;++k) {
        const auto species=active[k];const double change=answer[i][k]-current[i].X[species];
        helium_change-=change;
        if(change<-current[i].X[species])damping=std::min(damping,-.9*current[i].X[species]/change);
      }
      if(helium_change<-current[i].X[2])damping=std::min(damping,-.9*current[i].X[2]/helium_change);
    }
    if(!std::isfinite(damping) || damping<=0)throw std::runtime_error("D transport: no positive Newton step");
    double change=0;
    for(std::size_t i=0;i<count;++i) {
      for(std::size_t k=0;k<3;++k) {
        const double delta=damping*(answer[i][k]-current[i].X[active[k]]);
        current[i].X[active[k]]+=delta;change=std::max(change,std::abs(delta));
      }
      double consumed=0;
      for(auto species:active)consumed+=current[i].X[species]-old[i].X[species];
      current[i].X[2]=old[i].X[2]-consumed;
      for(double v:current[i].X)if(!std::isfinite(v) || v<0)
        throw std::runtime_error("D transport: invalid chemical iterate");
    }
    if(change<tolerance) {
      // Sum the source balance independently of the block elimination.
      Vec balance{};
      for(std::size_t i=0;i<count;++i) {
        for(std::size_t k=0;k<3;++k)balance[k]+=mass[i]*(current[i].X[active[k]]-old[i].X[active[k]]);
        for(std::size_t j=regions[i].first;j<regions[i].second;++j) {
          const auto rates=nuclear.eval(thermal.T(j),thermal.rho(j),current[i]);
          for(std::size_t k=0;k<3;++k)balance[k]-=dt*weights[j]/thermal.M*rates.dXdt[active[k]];
        }
      }
      if(norm(balance)>10*tolerance)throw std::runtime_error("D transport: integrated reaction imbalance");
      auto result=previous.comp;
      for(std::size_t i=0;i<count;++i)for(std::size_t j=regions[i].first;j<regions[i].second;++j)result[j]=current[i];
      return result;
    }
  }
  throw std::runtime_error("D transport: abundance iteration limit");
}
}

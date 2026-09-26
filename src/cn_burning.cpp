#include "ember/cn_burning.hpp"
#include "ember/constants.hpp"
#include "cn_source.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
using Vec=std::array<double,4>;
using Mat=std::array<Vec,4>;
Vec solve(Mat a,Vec b) {
  for(std::size_t k=0;k<4;++k) {
    std::size_t pivot=k;
    for(std::size_t j=k+1;j<4;++j)if(std::abs(a[j][k])>std::abs(a[pivot][k]))pivot=j;
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    if(!std::isfinite(a[k][k]) || a[k][k]==0)throw std::runtime_error("CN burning: singular Jacobian");
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
    if(!std::isfinite(x[i]))throw std::runtime_error("CN burning: nonfinite Newton correction");
  }
  return x;
}
// Express catalyst residuals as mass fractions when comparing with H/He.
double norm(const Vec& f) {return std::max({std::abs(f[0]),std::abs(f[1]),12*std::abs(f[2]),13*std::abs(f[3])});}
}
CNBurnResult burn_cn_and_mix(const Model& thermal,const Model& previous,
    const std::vector<CNAbundances>& old_cn,const PPChains& pp,const CNNetwork& cn,
    const MixingRegions& regions,double dt,double tolerance) {
  const auto weights=nodal_mass_weights(thermal);
  if(thermal.m!=previous.m || thermal.comp.size()!=thermal.size() || previous.comp.size()!=thermal.size()
      || old_cn.size()!=thermal.size() || !std::isfinite(dt) || dt<=0
      || !std::isfinite(tolerance) || tolerance<=0)
    throw std::invalid_argument("CN burning: inconsistent mesh, composition or timestep");
  CNBurnResult result;result.lookup=previous.comp;result.catalysts=old_cn;
  for(std::size_t i=0;i<old_cn.size();++i) {
    if(previous.comp[i][Species::H2]!=0)
      throw std::invalid_argument("CN burning: initial D requires the physical seven-species solve");
    if(previous.comp[i].cn_mass_convention!=CNMassConvention::fixed_metal_proxy)
      throw std::invalid_argument("CN burning: explicit metal mass requires the six-species solve");
    (void)cn_physical_ledger(previous.comp[i],old_cn[i]);
    for(std::size_t j=3;j<NSPEC;++j)
      if(previous.comp[i].X[j]!=previous.comp.front().X[j])
        throw std::invalid_argument("CN burning: fixed uniform metal lookup required");
  }
  std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())throw std::invalid_argument("CN burning: invalid mixing partition");
    next=end;Vec old{};double mass=0,total=0;
    for(std::size_t i=begin;i<end;++i) {
      mass+=weights[i];old[0]+=weights[i]*previous.comp[i].X[0];old[1]+=weights[i]*previous.comp[i].X[1];
      old[2]+=weights[i]*old_cn[i][0];old[3]+=weights[i]*old_cn[i][1];
      total+=weights[i]*(old_cn[i][0]+old_cn[i][1]+old_cn[i][2]);
    }
    for(double& v:old)v/=mass;total/=mass;
    auto lookup=previous.comp[begin];const double available=1-lookup.Z();
    auto composition=[&](const Vec& u) {
      auto c=lookup;c.X[0]=u[0];c.X[1]=u[1];c.X[2]=available-u[0]-u[1];return c;
    };
    auto catalysts=[&](const Vec& u)->CNAbundances{return {u[2],u[3],total-u[2]-u[3]};};
    auto positive=[&](const Vec& u) {
      for(double v:u)if(!std::isfinite(v) || v<0)return false;
      if(u[0]+u[1]>available || u[2]+u[3]>total)return false;
      const auto y=catalysts(u);const auto initial=initial_gs98_cn(lookup);
      const double extra=12*(y[0]-initial[0])+13*y[1]+14*(y[2]-initial[2]);
      return available-u[0]-u[1]-extra>=0;
    };
    auto equation=[&](const Vec& u,Mat* jacobian) {
      Vec f{};for(std::size_t row=0;row<4;++row)f[row]=u[row]-old[row];
      if(jacobian){*jacobian={};for(std::size_t j=0;j<4;++j)(*jacobian)[j][j]=1;}
      const auto c=composition(u);const auto y=catalysts(u);
      for(std::size_t i=begin;i<end;++i) {
        const auto rates=detail::cn_reduced_source(thermal.T(i),thermal.rho(i),c,y,pp,cn);
        const double factor=dt*weights[i]/mass;
        for(std::size_t row=0;row<4;++row) {
          f[row]-=factor*rates.source[row];
          if(jacobian)for(std::size_t col=0;col<4;++col)
            (*jacobian)[row][col]-=factor*rates.jacobian[row][col];
        }
      }
      return f;
    };
    Vec u=old;bool converged=false;
    for(std::size_t iteration=0;iteration<60;++iteration) {
      Mat j;const auto f=equation(u,&j);const double residual=norm(f);
      result.maximum_iterations=std::max(result.maximum_iterations,iteration+1);
      if(!std::isfinite(residual))throw std::runtime_error("CN burning: nonfinite residual");
      if(residual<tolerance){converged=true;result.maximum_equation_residual=std::max(result.maximum_equation_residual,residual);break;}
      Vec minus{};for(std::size_t k=0;k<4;++k)minus[k]=-f[k];
      const auto correction=solve(j,minus);double damping=1;bool accepted=false;
      for(int trial=0;trial<50;++trial,damping*=.5) {
        Vec candidate=u;for(std::size_t k=0;k<4;++k)candidate[k]+=damping*correction[k];
        if(!positive(candidate))continue;
        if(norm(equation(candidate,nullptr))<residual){u=candidate;accepted=true;break;}
      }
      if(!accepted)throw std::runtime_error("CN burning: abundance line search failed");
    }
    if(!converged)throw std::runtime_error("CN burning: abundance iteration limit");
    for(std::size_t i=begin;i<end;++i) {
      result.lookup[i]=composition(u);result.catalysts[i]=catalysts(u);
      result.lookup[i].cn_molality=result.catalysts[i];
    }
  }
  if(next!=thermal.size())throw std::invalid_argument("CN burning: incomplete mixing partition");
  double mass_release=0;
  for(std::size_t i=0;i<thermal.size();++i) {
    const auto p=pp.eval(thermal.T(i),thermal.rho(i),result.lookup[i]);
    const auto n=cn.response(thermal.T(i),thermal.rho(i),result.lookup[i],result.catalysts[i]).physical.state;
    result.nuclear_luminosity+=weights[i]*(p.eps+n.eps);
    result.neutrino_luminosity+=weights[i]*(p.eps_neutrino+n.eps_neutrino);
    double delta=0;
    for(std::size_t j=0;j<3;++j)
      delta+=(nuclides[j].A/mass_numbers[j]-1)*(result.lookup[i].X[j]-previous.comp[i].X[j])
        *constants::c*constants::c;
    delta+=cn_physical_ledger(result.lookup[i],result.catalysts[i]).rest_energy_difference
      -cn_physical_ledger(previous.comp[i],old_cn[i]).rest_energy_difference;
    mass_release-=weights[i]*delta/dt;
  }
  const double released=result.nuclear_luminosity+result.neutrino_luminosity;
  result.nuclear_mass_balance=released>0?mass_release/released-1:0;
  return result;
}
} // namespace ember

#include "ember/eos_smooth_mixture.hpp"
#include "ember/constants.hpp"
#include "composition_spline.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace ember {
namespace {
using detail::SplineFactor;
double xlogx(double n) {return n>0?n*std::log(n):0.;}
double full_mixing(const Composition& c) {
  const double n1=c.X[0],n3=c.X[1]/3,n4=c.X[2]/4;
  return constants::R_gas*(xlogx(n1)+xlogx(n3)+xlogx(n4)
      -1.5*n3*std::log(nuclides[1].A/nuclides[2].A));
}
double source_mixing(const Composition& c) {
  return constants::R_gas*(xlogx(c.X[0])+xlogx(c.X[1]/3+c.X[2]/4));
}

// H Hermite basis in order: left value, left slope, right value, right slope.
std::array<double,4> hydrogen_weights(double u,double h,unsigned d) {
  if(d==0)return {1-3*u*u+2*u*u*u,h*(u-2*u*u+u*u*u),
                  3*u*u-2*u*u*u,h*(-u*u+u*u*u)};
  if(d==1)return {(-6*u+6*u*u)/h,1-4*u+3*u*u,(6*u-6*u*u)/h,-2*u+3*u*u};
  if(d==2)return {(-6+12*u)/(h*h),(-4+6*u)/h,(6-12*u)/(h*h),(-2+6*u)/h};
  throw std::invalid_argument("unsupported composition derivative");
}
} // namespace

SmoothMetalHelmholtzEos::SmoothMetalHelmholtzEos(const std::filesystem::path& path,
    HelmholtzTableEos::Mixture mixture):residual_(path,mixture) {
  const auto nx=residual_.x_.size(),ny=residual_.y_.size();
  if(nx<4 || (ny!=3 && ny!=4))
    throw std::invalid_argument("smooth EOS requires >=4 H planes and 3 or 4 He3 planes");
  for(auto& p:residual_.tables_) {
    const double mix=source_mixing(p->composition_);
    for(auto& node:p->nodes_)if(node.valid)node.d[0]-=mix;
    auto slope=std::make_unique<HelmholtzTableEos>(*p);
    for(auto& node:slope->nodes_)node={};
    slopes_.push_back(std::move(slope));
  }
  std::map<std::pair<std::size_t,std::size_t>,SplineFactor> factors;
  std::vector<double> values(nx),derivatives(nx),seconds(nx);
  for(std::size_t iy=0;iy<ny;++iy)
    for(std::size_t node=0;node<residual_.tables_[iy]->nodes_.size();++node) {
      std::size_t begin=0;
      while(begin<nx) {
        while(begin<nx && !residual_.table(begin,iy).nodes_[node].valid)++begin;
        auto end=begin;
        while(end<nx && residual_.table(end,iy).nodes_[node].valid)++end;
        const auto count=end-begin;
        if(count>=4) {
          const auto key=std::pair{begin,end};
          auto f=factors.find(key);
          if(f==factors.end())f=factors.emplace(key,SplineFactor(
              std::span<const double>(residual_.x_).subspan(begin,count))).first;
          for(std::size_t k=0;k<9;++k) {
            for(std::size_t i=0;i<count;++i)values[i]=residual_.table(begin+i,iy).nodes_[node].d[k];
            f->second.slopes(std::span<const double>(values).first(count),
                std::span<double>(derivatives).first(count),std::span<double>(seconds).first(count));
            for(std::size_t i=0;i<count;++i) {
              if(!std::isfinite(derivatives[i]))throw std::runtime_error("nonfinite composition slope");
              auto& output=slopes_[(begin+i)*ny+iy]->nodes_[node];
              output.d[k]=derivatives[i];output.valid=true;
            }
          }
        }
        begin=end<nx?end+1:end;
      }
    }
  // Slope masks can be narrower than source masks near an isolated fragment.
  // Density inversion must see exactly the same contiguous support as eval.
  for(auto& p:slopes_)p->initialize_support();
}

SmoothMetalHelmholtzEos::WeightedTables
SmoothMetalHelmholtzEos::weights(const Composition& c,unsigned dx,unsigned dy) const {
  const auto q=residual_.coordinates(c);
  const double h=residual_.x_[q.x+1]-residual_.x_[q.x];
  const auto wx=hydrogen_weights(q.u,h,dx);
  const auto& y=residual_.y_;
  const double value=std::clamp(c.X[1],y.front(),y.back());
  if(dy>2)throw std::invalid_argument("unsupported composition derivative");
  WeightedTables result{};result.count=4*y.size();
  for(std::size_t j=0;j<y.size();++j) {
    double wy{};
    if(y.size()==3) {
      // Preserve the original three-plane arithmetic.
      const auto a=(j+1)%3,b=(j+2)%3;
      const double denominator=(y[j]-y[a])*(y[j]-y[b]);
      if(dy==0)wy=(value-y[a])*(value-y[b])/denominator;
      else if(dy==1)wy=(2*value-y[a]-y[b])/denominator;
      else wy=2/denominator;
    } else {
      // Four-node not-a-knot interpolation is a single cubic. Differentiate
      // its Lagrange products directly; this remains finite at each node.
      double p=1,first=0,second=0,denominator=1;
      for(std::size_t k=0;k<y.size();++k)if(k!=j) {
        const double factor=value-y[k];
        second=second*factor+2*first;first=first*factor+p;p*=factor;
        denominator*=y[j]-y[k];
      }
      wy=(dy==0?p:dy==1?first:second)/denominator;
    }
    for(std::size_t i=0;i<2;++i) {
      result.tables[j*4+i*2]={&residual_.table(q.x+i,j),wx[2*i]*wy};
      result.tables[j*4+i*2+1]={slopes_[(q.x+i)*y.size()+j].get(),wx[2*i+1]*wy};
    }
  }
  return result;
}

HelmholtzJet SmoothMetalHelmholtzEos::residual_jet(double T,double rho,const Composition& c,
    unsigned dx,unsigned dy) const {
  const auto w=weights(c,dx,dy);
  return HelmholtzTableEos::mixed_material_jet(T,rho,w.span());
}

std::optional<Eos::DensityRange> SmoothMetalHelmholtzEos::density_range(double T,const Composition& c) const {
  DensityRange result{0,std::numeric_limits<double>::infinity()};
  const auto tables=weights(c,0,0);
  for(const auto& w:tables.span()) {
    const auto range=w.table->material_density_range(T);
    result.min=std::max(result.min,range.min);result.max=std::min(result.max,range.max);
  }
  if(result.min>=result.max)throw std::domain_error("smooth EOS: empty source overlap");
  return result;
}

std::optional<Eos::DensityRange> SmoothMetalHelmholtzEos::density_range_near(double T,const Composition& c,double rho) const {
  DensityRange result{0,std::numeric_limits<double>::infinity()};
  const auto tables=weights(c,0,0);
  for(const auto& w:tables.span()) {
    const auto range=w.table->material_density_range_near(T,rho);
    result.min=std::max(result.min,range.min);result.max=std::min(result.max,range.max);
  }
  if(result.min>=result.max)throw std::domain_error("smooth EOS: empty source overlap");
  return result;
}

EosResponse SmoothMetalHelmholtzEos::eval_with_derivatives(double T,double rho,const Composition& c) const {
  auto f=residual_jet(T,rho,c,0,0);f[0][0]+=full_mixing(c);
  return helmholtz_response(T,rho,f);
}

EosCompositionResponse SmoothMetalHelmholtzEos::composition_response(double T,double rho,const Composition& c) const {
  EosCompositionResponse result{};
  for(unsigned k=0;k<2;++k) {
    const auto f=residual_jet(T,rho,c,k==0?1:0,k==1?1:0);
    result.dP[k]=rho*T*f[0][1];result.dE[k]=-T*f[1][0];
  }
  return result;
}

CompositionPotentialResponse SmoothMetalHelmholtzEos::composition_potential(double T,double rho,const Composition& c) const {
  const auto f=residual_jet(T,rho,c,0,0);
  if(c.X[0]<=0 || c.X[1]<=0 || c.X[2]<=0)
    throw std::domain_error("chemical derivatives require positive active H1, He3 and He4 populations");
  CompositionPotentialResponse result{};result.phi=f[0][0]+full_mixing(c);
  const std::array<double,3> n{c.X[0],c.X[1]/3,c.X[2]/4};
  constexpr double b[3][2]={{1,0},{0,1./3},{-.25,-.25}};
  for(unsigned k=0;k<2;++k) {
    const auto g=residual_jet(T,rho,c,k==0?1:0,k==1?1:0);
    result.gradient[k]=g[0][0];result.dgradient_dlnT[k]=g[1][0];result.dgradient_dlnRho[k]=g[0][1];
    for(unsigned i=0;i<3;++i)result.gradient[k]+=constants::R_gas*b[i][k]*(std::log(n[i])+1);
    for(unsigned l=k;l<2;++l) {
      const auto second=residual_jet(T,rho,c,(k==0?1:0)+(l==0?1:0),(k==1?1:0)+(l==1?1:0));
      double value=second[0][0];
      for(unsigned i=0;i<3;++i)value+=constants::R_gas*b[i][k]*b[i][l]/n[i];
      result.hessian[k][l]=result.hessian[l][k]=value;
    }
  }
  result.gradient[1]-=.5*constants::R_gas*std::log(nuclides[1].A/nuclides[2].A);
  return result;
}

CompositionPotentialResponse SmoothMetalHelmholtzEos::regular_composition_potential(
    double T,double rho,const Composition& c) const {
  const auto f=residual_jet(T,rho,c,0,0);
  if(c.X[2]<=0)throw std::domain_error("regular chemical derivatives require positive reference He4");
  CompositionPotentialResponse result;result.phi=f[0][0]+full_mixing(c);
  const double reference=-.25*constants::R_gas*(std::log(c.X[2]/4)+1);
  const double reference_curvature=.25*constants::R_gas/c.X[2];
  constexpr double mass[2]={1,3};
  for(unsigned k=0;k<2;++k) {
    const auto g=residual_jet(T,rho,c,k==0?1:0,k==1?1:0);
    result.gradient[k]=g[0][0]+constants::R_gas/mass[k]*(1-std::log(mass[k]))+reference;
    result.dgradient_dlnT[k]=g[1][0];result.dgradient_dlnRho[k]=g[0][1];
    for(unsigned l=k;l<2;++l) {
      const auto second=residual_jet(T,rho,c,(k==0?1:0)+(l==0?1:0),(k==1?1:0)+(l==1?1:0));
      result.hessian[k][l]=result.hessian[l][k]=second[0][0]+reference_curvature;
    }
  }
  result.gradient[1]-=.5*constants::R_gas*std::log(nuclides[1].A/nuclides[2].A);
  return result;
}

CompositionPotentialResponse SmoothMetalHelmholtzEos::active_composition_potential(
    double T,double rho,const Composition& c,std::array<bool,2> active) const {
  if(active[0] && active[1])return composition_potential(T,rho,c);
  const auto f=residual_jet(T,rho,c,0,0);
  if((active[0] && c.X[0]<=0) || (active[1] && c.X[1]<=0) ||
      ((active[0] || active[1]) && c.X[2]<=0))
    throw std::domain_error("chemical derivatives require positive active H1, He3 and He4 populations");
  CompositionPotentialResponse result{};result.phi=f[0][0]+full_mixing(c);
  const double absent=std::numeric_limits<double>::quiet_NaN();
  result.gradient.fill(absent);result.dgradient_dlnT.fill(absent);
  result.dgradient_dlnRho.fill(absent);
  for(auto& row:result.hessian)row.fill(absent);
  const std::array<double,3> n{c.X[0],c.X[1]/3,c.X[2]/4};
  constexpr double b[3][2]={{1,0},{0,1./3},{-.25,-.25}};
  for(unsigned k=0;k<2;++k)if(active[k]) {
    const auto g=residual_jet(T,rho,c,k==0?1:0,k==1?1:0);
    result.gradient[k]=g[0][0];result.dgradient_dlnT[k]=g[1][0];result.dgradient_dlnRho[k]=g[0][1];
    // Absent, unselected ions contribute neither log(0) nor 0/0. Retain the
    // original arithmetic for all positive populations, including zero b.
    for(unsigned i=0;i<3;++i)if(n[i]>0)
      result.gradient[k]+=constants::R_gas*b[i][k]*(std::log(n[i])+1);
    for(unsigned l=k;l<2;++l)if(active[l]) {
      const auto second=residual_jet(T,rho,c,(k==0?1:0)+(l==0?1:0),(k==1?1:0)+(l==1?1:0));
      double value=second[0][0];
      for(unsigned i=0;i<3;++i)if(n[i]>0)value+=constants::R_gas*b[i][k]*b[i][l]/n[i];
      result.hessian[k][l]=result.hessian[l][k]=value;
    }
  }
  if(active[1])result.gradient[1]-=.5*constants::R_gas*std::log(nuclides[1].A/nuclides[2].A);
  return result;
}

CompositionHeatResponse SmoothMetalHelmholtzEos::composition_heat(double T,double rho,
    const Composition& c,std::array<bool,2> active,bool derivatives) const {
  const auto f=residual_jet(T,rho,c,0,0);
  // Validate the same caloric/stability conditions used by the structural EOS.
  // The ideal mixing addition is constant in T,rho and cannot alter these.
  (void)helmholtz_response(T,rho,f);
  if((active[0] && c.X[0]<=0) || (active[1] && c.X[1]<=0) ||
      ((active[0] || active[1]) && c.X[2]<=0))
    throw std::domain_error("enthalpy derivatives require positive active species and He4");
  const double denominator=f[0][1]+f[0][2];
  if(!(denominator>0) || !std::isfinite(denominator))
    throw std::domain_error("composition heat: nonpositive material compressibility");
  CompositionHeatResponse result;
  const double delta=(f[0][1]+f[1][1])/denominator;
  const double radiation_delta=4*constants::a_rad*T*T*T/(3*rho*denominator);
  if(!std::isfinite(delta))throw std::domain_error("composition heat: invalid thermal expansion");
  if(!std::isfinite(radiation_delta))throw std::domain_error("composition heat: invalid radiation expansion");
  result.material_delta=delta;
  const double absent=std::numeric_limits<double>::quiet_NaN();
  result.exchange_enthalpy.fill(absent);result.delta_partials.fill(absent);
  result.radiation_enthalpy.fill(absent);
  for(auto& row:result.enthalpy_partials)row.fill(absent);
  for(auto& row:result.radiation_enthalpy_partials)row.fill(absent);
  std::array<HelmholtzJet,2> gradients{};
  for(unsigned k=0;k<2;++k)if(active[k]) {
    const auto g=residual_jet(T,rho,c,k==0?1:0,k==1?1:0);
    gradients[k]=g;
    result.exchange_enthalpy[k]=T*(delta*g[0][1]-g[1][0]);
    result.radiation_enthalpy[k]=T*radiation_delta*g[0][1];
    if(!std::isfinite(result.exchange_enthalpy[k]))throw std::domain_error("composition heat: nonfinite enthalpy");
    if(!std::isfinite(result.radiation_enthalpy[k]))throw std::domain_error("composition heat: nonfinite radiation enthalpy");
  }
  if(!derivatives)return result;
  auto& d=result.delta_partials;
  std::array<double,4> dr{radiation_delta*(3-(f[1][1]+f[1][2])/denominator),
                        radiation_delta*(-1-(f[0][2]+f[0][3])/denominator),absent,absent};
  d[0]=(f[1][1]+f[2][1]-delta*(f[1][1]+f[1][2]))/denominator;
  d[1]=(f[0][2]+f[1][2]-delta*(f[0][2]+f[0][3]))/denominator;
  if(!std::isfinite(d[0]) || !std::isfinite(d[1]))throw std::domain_error("composition heat: nonfinite thermal derivative");
  for(unsigned k=0;k<2;++k)if(active[k]) {
    const auto& g=gradients[k];
    d[2+k]=(g[0][1]+g[1][1]-delta*(g[0][1]+g[0][2]))/denominator;
    dr[2+k]=-radiation_delta*(g[0][1]+g[0][2])/denominator;
    auto& h=result.enthalpy_partials[k];
    h[0]=result.exchange_enthalpy[k]+T*(d[0]*g[0][1]+delta*g[1][1]-g[2][0]);
    h[1]=T*(d[1]*g[0][1]+delta*g[0][2]-g[1][1]);
    auto& hr=result.radiation_enthalpy_partials[k];
    hr[0]=result.radiation_enthalpy[k]+T*(dr[0]*g[0][1]+radiation_delta*g[1][1]);
    hr[1]=T*(dr[1]*g[0][1]+radiation_delta*g[0][2]);
  }
  // Reuse each symmetric composition Hessian jet for both enthalpy rows.
  for(unsigned k=0;k<2;++k)if(active[k])for(unsigned l=k;l<2;++l)if(active[l]) {
    const auto second=residual_jet(T,rho,c,(k==0?1:0)+(l==0?1:0),(k==1?1:0)+(l==1?1:0));
    const double fixed_delta=T*(delta*second[0][1]-second[1][0]);
    result.enthalpy_partials[k][2+l]=fixed_delta+T*d[2+l]*gradients[k][0][1];
    result.enthalpy_partials[l][2+k]=fixed_delta+T*d[2+k]*gradients[l][0][1];
    const double radiation_fixed_delta=T*radiation_delta*second[0][1];
    result.radiation_enthalpy_partials[k][2+l]=radiation_fixed_delta+T*dr[2+l]*gradients[k][0][1];
    result.radiation_enthalpy_partials[l][2+k]=radiation_fixed_delta+T*dr[2+k]*gradients[l][0][1];
  }
  for(unsigned k=0;k<2;++k)if(active[k]) {
    if(!std::isfinite(result.exchange_enthalpy[k]))throw std::domain_error("composition heat: nonfinite enthalpy");
    for(unsigned j=0;j<4;++j)if(j<2 || active[j-2])
      if(!std::isfinite(result.enthalpy_partials[k][j]) || !std::isfinite(d[j])
          || !std::isfinite(result.radiation_enthalpy_partials[k][j]) || !std::isfinite(dr[j]))
        throw std::domain_error("composition heat: nonfinite derivative");
  }
  return result;
}
} // namespace ember

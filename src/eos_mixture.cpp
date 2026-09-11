#include "ember/eos_mixture.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ember {
MetalHelmholtzEos::MetalHelmholtzEos(const std::filesystem::path& path,HelmholtzTableEos::Mixture mixture) {
  if(mixture!=HelmholtzTableEos::Mixture::allow_documented_proxy)
    throw std::invalid_argument("MetalHelmholtzEos: explicit source approximation selection required");
  std::ifstream in(path);std::string key;int version{};
  in>>key>>version;
  if(!in || key!="EMBER_METAL_HELMHOLTZ" || version!=1)
    throw std::runtime_error("MetalHelmholtzEos: invalid manifest");
  auto axis=[&](const char* label,std::vector<double>& values) {
    std::size_t count{};in>>key>>count;
    if(!in || key!=label || count<2 || count>100)throw std::runtime_error("MetalHelmholtzEos: invalid axis");
    values.resize(count);
    for(std::size_t i=0;i<count;++i) {
      in>>values[i];
      if(!in || !std::isfinite(values[i]) || values[i]<0 || values[i]>.98 || (i && values[i]<=values[i-1]))
        throw std::runtime_error("MetalHelmholtzEos: invalid composition coordinate");
    }
  };
  axis("hydrogen",x_);axis("helium3",y_);
  for(std::size_t i=0;i<x_.size();++i)for(std::size_t j=0;j<y_.size();++j) {
    std::string file;in>>std::quoted(file);
    if(!in || file.empty())throw std::runtime_error("MetalHelmholtzEos: missing plane");
    auto p=std::make_unique<HelmholtzTableEos>(path.parent_path()/file,mixture);
    if(!tables_.empty() && !p->same_material_grid(*tables_.front()))
      throw std::runtime_error("MetalHelmholtzEos: material grids must agree across composition planes");
    const auto& c=p->composition();
    if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98
        || std::abs(c.X[0]-x_[i])>1e-12 || std::abs(c.X[1]-y_[j])>1e-12)
      throw std::runtime_error("MetalHelmholtzEos: source composition mismatch");
    if(tables_.empty())metals_=c;
    for(std::size_t k=3;k<NSPEC;++k)if(std::abs(c.X[k]-metals_.X[k])>1e-12)
      throw std::runtime_error("MetalHelmholtzEos: inconsistent metal inventory");
    tables_.push_back(std::move(p));
  }
  if(in>>key)throw std::runtime_error("MetalHelmholtzEos: trailing manifest data");
}

MetalHelmholtzEos::Coordinates MetalHelmholtzEos::coordinates(const Composition& c) const {
  if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98 || std::abs(c.sum()-1)>1e-10)
    throw std::domain_error("MetalHelmholtzEos: normalized baryonic composition required");
  for(std::size_t k=0;k<NSPEC;++k)
    if(!std::isfinite(c.X[k]) || c.X[k]<0 || (k>=3 && std::abs(c.X[k]-metals_.X[k])>1e-10))
      throw std::domain_error("MetalHelmholtzEos: invalid composition or changed metals");
  auto bracket=[](double v,const std::vector<double>& axis) {
    if(v<axis.front()-2e-14 || v>axis.back()+2e-14)
      throw std::domain_error("MetalHelmholtzEos: composition outside family");
    v=std::clamp(v,axis.front(),axis.back());
    const auto i=interp::locate(axis,v);
    return std::pair{i,(v-axis[i])/(axis[i+1]-axis[i])};
  };
  const auto [i,u]=bracket(c.X[0],x_);const auto [j,v]=bracket(c.X[1],y_);
  return {i,j,u,v};
}

std::size_t MetalHelmholtzEos::check_temperature_extension(const MetalHelmholtzEos& old) const {
  if(x_!=old.x_ || y_!=old.y_)
    throw std::runtime_error("EOS extension: composition axes changed");
  auto physical_source=[](const std::string& source) {
    const std::string label="; direct-source SHA256 ";
    const auto at=source.rfind(label);
    if(at==std::string::npos)return source;
    const auto hash=source.substr(at+label.size());
    if(hash.size()!=64 || hash.find_first_not_of("0123456789abcdef")!=std::string::npos)
      throw std::runtime_error("EOS extension: malformed raw source identity");
    return source.substr(0,at);
  };
  std::size_t added=0;
  for(std::size_t k=0;k<tables_.size();++k) {
    const auto& a=*tables_[k];const auto& b=*old.tables_[k];
    if(a.proxy_!=b.proxy_ || physical_source(a.source_)!=physical_source(b.source_)
        || a.composition_.X!=b.composition_.X || a.composition_.basis!=b.composition_.basis
        || a.composition_.metal_inventory!=b.composition_.metal_inventory)
      throw std::runtime_error("EOS extension: source physics or composition changed");
    if(a.q_!=b.q_ || a.t_.size()<=b.t_.size()
        || !std::equal(b.t_.begin(),b.t_.end(),a.t_.begin()))
      throw std::runtime_error("EOS extension: original material axes changed");
    for(std::size_t i=0;i<b.nodes_.size();++i)
      if(a.nodes_[i].valid!=b.nodes_[i].valid || a.nodes_[i].d!=b.nodes_[i].d)
        throw std::runtime_error("EOS extension: original potential values or masks changed");
    for(std::size_t i=b.nodes_.size();i<a.nodes_.size();++i)added+=a.nodes_[i].valid;
  }
  if(!added)throw std::runtime_error("EOS extension: no valid hotter states added");
  return added;
}

std::optional<Eos::DensityRange> MetalHelmholtzEos::density_range(double T,const Composition& c) const {
  const auto q=coordinates(c);DensityRange r{0,std::numeric_limits<double>::infinity()};
  for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j) {
    const auto p=table(q.x+i,q.y+j).material_density_range(T);
    r.min=std::max(r.min,p.min);r.max=std::min(r.max,p.max);
  }
  if(r.min>=r.max)throw std::domain_error("MetalHelmholtzEos: empty plane support overlap");
  return r;
}

EosResponse MetalHelmholtzEos::eval_with_derivatives(double T,double rho,const Composition& c) const {
  const auto q=coordinates(c);std::array<HelmholtzTableEos::WeightedTable,4> planes{};
  for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j) {
    const double w=(i?q.u:1-q.u)*(j?q.v:1-q.v);
    planes[2*i+j]={&table(q.x+i,q.y+j),w};
  }
  auto f=HelmholtzTableEos::mixed_material_jet(T,rho,planes);
  const double n3=c.X[1]/3,n4=c.X[2]/4,he=n3+n4;
  const auto xlog=[](double n,double sum){return n>0?n*std::log(n/sum):0.;};
  f[0][0]+=constants::R_gas*(xlog(n3,he)+xlog(n4,he)-1.5*n3*std::log(nuclides[1].A/nuclides[2].A));
  return helmholtz_response(T,rho,f);
}

EosCompositionResponse MetalHelmholtzEos::composition_response(double T,double rho,const Composition& c) const {
  const auto q=coordinates(c);EosCompositionResponse out{};
  for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j) {
    const auto p=table(q.x+i,q.y+j).material_jet(T,rho);
    const std::array<double,2> w{(i?1.:-1.)*(j?q.v:1-q.v)/(x_[q.x+1]-x_[q.x]),
                               (j?1.:-1.)*(i?q.u:1-q.u)/(y_[q.y+1]-y_[q.y])};
    for(int k=0;k<2;++k) {out.dP[k]+=rho*T*w[k]*p[0][1];out.dE[k]-=T*w[k]*p[1][0];}
  }
  return out;
}
} // namespace ember

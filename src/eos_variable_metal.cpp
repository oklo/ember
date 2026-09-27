#include "ember/eos_variable_metal.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include "composition_spline.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>

namespace ember {
namespace {
constexpr std::size_t second_channel(std::size_t i,std::size_t j) {
  if(i>j)std::swap(i,j);
  return 4+i*3-i*(i-1)/2+j-i;
}
// Second-order chain rule only for the three composition coordinates.
// Thermodynamic derivatives come from the source potential independently.
struct CJet {
  double value{};std::array<double,3> d{};std::array<std::array<double,3>,3> h{};
  CJet()=default;CJet(double v):value(v){}
  static CJet variable(double v,std::size_t k){CJet a(v);a.d[k]=1;return a;}
  friend CJet operator+(const CJet& a,const CJet& b) {
    CJet c(a.value+b.value);
    for(std::size_t i=0;i<3;++i){c.d[i]=a.d[i]+b.d[i];for(std::size_t j=0;j<3;++j)c.h[i][j]=a.h[i][j]+b.h[i][j];}
    return c;
  }
  friend CJet operator-(const CJet& a,const CJet& b) {
    CJet c(a.value-b.value);
    for(std::size_t i=0;i<3;++i){c.d[i]=a.d[i]-b.d[i];for(std::size_t j=0;j<3;++j)c.h[i][j]=a.h[i][j]-b.h[i][j];}
    return c;
  }
  friend CJet operator*(const CJet& a,const CJet& b) {
    CJet c(a.value*b.value);
    for(std::size_t i=0;i<3;++i){c.d[i]=a.d[i]*b.value+a.value*b.d[i];
      for(std::size_t j=0;j<3;++j)c.h[i][j]=a.h[i][j]*b.value+a.d[i]*b.d[j]+a.d[j]*b.d[i]+a.value*b.h[i][j];}
    return c;
  }
  friend CJet operator/(const CJet& a,const CJet& b) {
    if(b.value==0)throw std::domain_error("variable EOS: zero coordinate denominator");
    CJet inverse(1/b.value);
    for(std::size_t i=0;i<3;++i){inverse.d[i]=-b.d[i]/(b.value*b.value);
      for(std::size_t j=0;j<3;++j)inverse.h[i][j]=(2*b.d[i]*b.d[j]/b.value-b.h[i][j])/(b.value*b.value);}
    return a*inverse;
  }
};
std::vector<CJet> lagrange(const CJet& x,std::span<const double> nodes) {
  std::vector<CJet> result(nodes.size(),CJet(1));
  for(std::size_t i=0;i<nodes.size();++i)for(std::size_t j=0;j<nodes.size();++j)if(i!=j)
    result[i]=result[i]*(x-nodes[j])/(nodes[i]-nodes[j]);
  return result;
}
double xlogx(double x){return x>0?x*std::log(x):0.;}
bool potassium(const Gs98Metal& m){return m.charge==19 && m.mass_number==39;}
double metal_mixing(double Z) {
  double sum=0;for(const auto& m:gs98_metals)if(!potassium(m))sum+=xlogx(Z*m.fraction/m.mass_number);
  return sum;
}
double source_mixing(const Composition& c) {
  return constants::R_gas*(xlogx(c.X[0])+xlogx(c.X[1]/3+c.X[2]/4)+metal_mixing(c.Z()));
}
double full_mixing(const Composition& c) {
  const double n3=c.X[1]/3;
  return constants::R_gas*(xlogx(c.X[0])+xlogx(n3)+xlogx(c.X[2]/4)+metal_mixing(c.Z())
      -1.5*n3*std::log(nuclides[1].A/nuclides[2].A));
}
bool active_any(std::array<bool,3> active){return active[0] || active[1] || active[2];}
void check_active(const Composition& c,std::array<bool,3> active) {
  if((active[0] && c.X[0]<=0) || (active[1] && c.X[1]<=0) || (active[2] && c.Z()<=0)
      || (active_any(active) && c.X[2]<=0))
    throw std::domain_error("variable EOS: positive active species and reference He4 required");
}
}

VariableMetalHelmholtzEos::VariableMetalHelmholtzEos(const std::filesystem::path& path,
    HelmholtzTableEos::Mixture mixture) {
  if(mixture!=HelmholtzTableEos::Mixture::allow_documented_proxy)
    throw std::invalid_argument("variable EOS: source approximations require explicit selection");
  std::ifstream in(path,std::ios::binary);std::string key;int version{};in>>key>>version;
  const bool binary=key=="EMBER_VARIABLE_METAL_HELMHOLTZ_BINARY";
  if(!in || (!binary && key!="EMBER_VARIABLE_METAL_HELMHOLTZ") || (version!=1 && version!=2))
    throw std::runtime_error("variable EOS: invalid manifest");
  extend_metals_=version==2;
  auto axis=[&](const char* label,std::vector<double>& values,std::size_t minimum,std::size_t maximum) {
    std::size_t n{};in>>key>>n;
    if(!in || key!=label || n<minimum || n>maximum)throw std::runtime_error("variable EOS: invalid axis");
    values.resize(n);
    for(std::size_t i=0;i<n;++i){in>>values[i];if(!in || !std::isfinite(values[i]) || values[i]<0
        || values[i]>1 || (i && values[i]<=values[i-1]))throw std::runtime_error("variable EOS: invalid axis coordinate");}
  };
  axis("metals",z_,extend_metals_?4:3,extend_metals_?32:4);
  axis("hydrogen_share",u_,4,100);axis("helium3_share",v_,3,4);
  if(binary) {
    in>>key;
    if(key!="planes_binary_v1" || in.get()!='\n')throw std::runtime_error("variable EOS: invalid binary plane marker");
  }
  if(z_.back()>=1)throw std::runtime_error("variable EOS: material needs H or helium");
  // Preserve the original cubic on the first four planes. Each subsequent
  // interval matches value, first derivative and second derivative at both
  // ends. At node i, derivatives use precisely nodes i-3..i, so a later
  // source extension cannot change an earlier interval or its support.
  constexpr double basis[6][6]={
    {1,0,0,-10,15,-6},{0,1,0,-6,8,-3},{0,0,.5,-1.5,1.5,-.5},
    {0,0,0,10,-15,6},{0,0,0,-4,7,-3},{0,0,0,.5,-1,.5}};
  if(extend_metals_)for(std::size_t i=3;i+1<z_.size();++i) {
    const auto left=lagrange(CJet::variable(z_[i],0),std::span<const double>(z_).subspan(i-3,4));
    const auto right=lagrange(CJet::variable(z_[i+1],0),std::span<const double>(z_).subspan(i-2,4));
    const double h=z_[i+1]-z_[i];MetalExtension extension{};
    for(std::size_t k=0;k<5;++k) {
      const CJet a=k<4?left[k]:CJet{},b=k>0?right[k-1]:CJet{};
      const std::array<double,6> endpoint{a.value,h*a.d[0],h*h*a.h[0][0],
                                         b.value,h*b.d[0],h*h*b.h[0][0]};
      for(std::size_t n=0;n<6;++n)for(std::size_t j=0;j<6;++j)
        extension[k][n]+=endpoint[j]*basis[j][n];
    }
    metal_extensions_.push_back(extension);
  }
  bool have_pattern=false;
  for(std::size_t iz=0;iz<z_.size();++iz)for(std::size_t iu=0;iu<u_.size();++iu)
    for(std::size_t iv=0;iv<v_.size();++iv) {
      std::unique_ptr<HelmholtzTableEos> table;
      if(binary)table.reset(new HelmholtzTableEos(in,mixture));
      else {
        std::string filename;in>>std::quoted(filename);
        if(!in || filename.empty())throw std::runtime_error("variable EOS: missing source plane");
        table=std::make_unique<HelmholtzTableEos>(path.parent_path()/filename,mixture);
      }
      if(!tables_.empty() && !table->same_material_grid(*tables_.front()))
        throw std::runtime_error("variable EOS: source material grids differ");
      const auto& c=table->composition();
      const double X=(1-z_[iz])*u_[iu],Y=(1-z_[iz])*(1-u_[iu])*v_[iv];
      if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98
          || std::abs(c.X[0]-X)>1e-12 || std::abs(c.X[1]-Y)>1e-12 || std::abs(c.Z()-z_[iz])>1e-12)
        throw std::runtime_error("variable EOS: source composition mismatch");
      if(c.Z()>0) {
        if(!have_pattern){for(std::size_t k=0;k<5;++k)metal_pattern_[k]=c.X[3+k]/c.Z();have_pattern=true;}
        for(std::size_t k=0;k<5;++k)if(std::abs(c.X[3+k]/c.Z()-metal_pattern_[k])>1e-12)
          throw std::runtime_error("variable EOS: source metal patterns differ");
      }
      const double mix=source_mixing(c);
      for(auto& node:table->nodes_)if(node.valid)node.d[0]-=mix;
      auto slope=std::make_unique<HelmholtzTableEos>(*table);
      for(auto& node:slope->nodes_)node={};
      tables_.push_back(std::move(table));slopes_.push_back(std::move(slope));
    }
  if(in>>key)throw std::runtime_error("variable EOS: trailing manifest data");
  std::map<std::pair<std::size_t,std::size_t>,detail::SplineFactor> factors;
  const auto nu=u_.size();std::vector<double> values(nu),slopes(nu),seconds(nu);
  for(std::size_t iz=0;iz<z_.size();++iz)for(std::size_t iv=0;iv<v_.size();++iv)
    for(std::size_t node=0;node<tables_[0]->nodes_.size();++node) {
      std::size_t begin=0;
      while(begin<nu) {
        while(begin<nu && !tables_[index(iz,begin,iv)]->nodes_[node].valid)++begin;
        auto end=begin;while(end<nu && tables_[index(iz,end,iv)]->nodes_[node].valid)++end;
        const auto count=end-begin;
        if(count>=4) {
          const auto interval=std::pair{begin,end};auto factor=factors.find(interval);
          if(factor==factors.end())factor=factors.emplace(interval,detail::SplineFactor(std::span<const double>(u_).subspan(begin,count))).first;
          for(std::size_t k=0;k<9;++k) {
            for(std::size_t i=0;i<count;++i)values[i]=tables_[index(iz,begin+i,iv)]->nodes_[node].d[k];
            factor->second.slopes(std::span<const double>(values).first(count),std::span<double>(slopes).first(count),std::span<double>(seconds).first(count));
            for(std::size_t i=0;i<count;++i) {
              if(!std::isfinite(slopes[i]))throw std::runtime_error("variable EOS: nonfinite composition slope");
              auto& out=slopes_[index(iz,begin+i,iv)]->nodes_[node];out.d[k]=slopes[i];out.valid=true;
            }
          }
        }
        begin=end<nu?end+1:end;
      }
    }
  for(auto& p:slopes_)for(std::size_t it=0;it+1<p->t_.size();++it) {
    std::size_t hi=0;while(hi<p->q_.size() && p->nodes_[it*p->q_.size()+hi].valid
        && p->nodes_[(it+1)*p->q_.size()+hi].valid)++hi;
    p->supported_hi_[it]=hi;
  }
}

VariableMetalHelmholtzEos::Weights VariableMetalHelmholtzEos::weights(const Composition& c,std::size_t channels) const {
  if(c[Species::H2]!=0)throw std::domain_error("variable EOS: select an explicit deuterium material approximation");
  if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98 || std::abs(c.sum()-1)>1e-10)
    throw std::domain_error("variable EOS: normalized baryonic GS98 composition required");
  for(double x:c.X)if(!std::isfinite(x) || x<0)throw std::domain_error("variable EOS: invalid abundance");
  const double metal=c.Z();
  if(metal>=1)throw std::domain_error("variable EOS: empty H/He material");
  for(std::size_t k=0;k<5;++k)if(std::abs(c.X[3+k]-metal*metal_pattern_[k])>1e-12)
    throw std::domain_error("variable EOS: changed relative metal pattern");
  const auto H=CJet::variable(c.X[0],0),Y=CJet::variable(c.X[1],1),Z=CJet::variable(metal,2);
  auto u=H/(1-Z);CJet v;
  if(c.X[1]+c.X[2]>0) {
    auto helium=1-H-Z;helium.value=c.X[1]+c.X[2];v=Y/helium;
  }else if(channels!=1)throw std::domain_error("variable EOS: pure-H composition derivatives require a different reference species");
  auto bracket=[](CJet& x,const std::vector<double>& axis) {
    if(x.value<axis.front()-2e-14 || x.value>axis.back()+2e-14)
      throw std::domain_error("variable EOS: composition outside source family");
    x.value=std::clamp(x.value,axis.front(),axis.back());
  };
  auto z=Z;bracket(u,u_);bracket(v,v_);bracket(z,z_);
  const auto iu=interp::locate(u_,u.value);const double h=u_[iu+1]-u_[iu];
  const auto q=(u-u_[iu])/h,q2=q*q,q3=q2*q;
  const std::array<CJet,4> wu{1-3*q2+2*q3,h*(q-2*q2+q3),3*q2-2*q3,h*(q3-q2)};
  const auto wv=lagrange(v,v_);std::vector<CJet> wz;std::size_t first_z=0;
  if(!extend_metals_ || z.value<=z_[3]) {
    wz=lagrange(z,std::span<const double>(z_).first(extend_metals_?4:z_.size()));
  }else {
    const auto iz=interp::locate(z_,z.value);first_z=iz-3;
    const auto qz=(z-z_[iz])/(z_[iz+1]-z_[iz]);
    const auto& polynomial=metal_extensions_[iz-3];wz.resize(5);
    for(std::size_t k=0;k<5;++k) {
      wz[k]=polynomial[k][5];
      for(int n=4;n>=0;--n)wz[k]=wz[k]*qz+polynomial[k][n];
    }
  }
  Weights out;
  for(std::size_t iz=0;iz<wz.size();++iz)for(std::size_t iv=0;iv<v_.size();++iv)
    for(std::size_t side=0;side<2;++side)for(std::size_t slope=0;slope<2;++slope) {
      const auto w=wz[iz]*wv[iv]*wu[2*side+slope];
      // A plane whose weight is exactly zero in every REQUESTED channel adds exactly zero to every
      // accumulated sum, so its table is not consulted for any value here. Emitting it anyway still
      // subjects the query to its density support, because mixed_composition_jets checks supported_q
      // for every plane in the span. That turns a masked cell in a plane the mixture does not use into
      // a hard domain limit. Lagrange weights vanish exactly on their own node, so this is not a rare
      // case: a star whose metal mass sits on a metals-axis node carries three such planes at all
      // times. Channels above the requested count are never read, so they cannot keep a plane alive.
      bool contributes=w.value!=0;
      if(channels>1)for(std::size_t k=0;k<3;++k)contributes|=w.d[k]!=0;
      if(channels>4)for(std::size_t k=0;k<3;++k)for(std::size_t l=k;l<3;++l)contributes|=w.h[k][l]!=0;
      if(!contributes)continue;
      auto& item=out.tables[out.count++];
      item.table=(slope?slopes_:tables_)[index(first_z+iz,iu+side,iv)].get();item.weight[0]=w.value;
      if(channels>1)for(std::size_t k=0;k<3;++k)item.weight[1+k]=w.d[k];
      if(channels>4)for(std::size_t k=0;k<3;++k)for(std::size_t l=k;l<3;++l)item.weight[second_channel(k,l)]=w.h[k][l];
    }
  return out;
}
std::array<HelmholtzJet,10> VariableMetalHelmholtzEos::jets(double T,double rho,const Composition& c,std::size_t channels) const {
  const auto w=weights(c,channels);
  return HelmholtzTableEos::mixed_composition_jets(T,rho,std::span<const WeightedTable>(w.tables).first(w.count),channels);
}
std::optional<Eos::DensityRange> VariableMetalHelmholtzEos::density_range(double T,const Composition& c) const {
  const auto w=weights(c,1);DensityRange result{0,std::numeric_limits<double>::infinity()};
  for(std::size_t i=0;i<w.count;++i){const auto r=w.tables[i].table->material_density_range(T);result.min=std::max(result.min,r.min);result.max=std::min(result.max,r.max);}
  if(result.min>=result.max)throw std::domain_error("variable EOS: empty source density overlap");
  return result;
}
EosResponse VariableMetalHelmholtzEos::eval_with_derivatives(double T,double rho,const Composition& c) const {
  auto f=jets(T,rho,c,1)[0];f[0][0]+=full_mixing(c);return helmholtz_response(T,rho,f);
}
EosCompositionResponse VariableMetalHelmholtzEos::composition_response(double T,double rho,const Composition& c) const {
  const auto j=jets(T,rho,c,4);EosCompositionResponse out;
  for(std::size_t k=0;k<2;++k){out.dP[k]=rho*T*j[1+k][0][1];out.dE[k]=-T*j[1+k][1][0];}
  return out;
}
MetalCompositionPotentialResponse VariableMetalHelmholtzEos::composition_potential(double T,double rho,
    const Composition& c,std::array<bool,3> active) const {
  check_active(c,active);const auto j=jets(T,rho,c,active_any(active)?10:1);
  MetalCompositionPotentialResponse out;out.phi=j[0][0][0]+full_mixing(c);
  const double nan=std::numeric_limits<double>::quiet_NaN();
  out.gradient.fill(nan);out.dgradient_dlnT.fill(nan);out.dgradient_dlnRho.fill(nan);for(auto& h:out.hessian)h.fill(nan);
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    out.gradient[k]=j[1+k][0][0];out.dgradient_dlnT[k]=j[1+k][1][0];out.dgradient_dlnRho[k]=j[1+k][0][1];
    for(std::size_t l=0;l<3;++l)if(active[l])out.hessian[k][l]=j[second_channel(k,l)][0][0];
  }
  auto add=[&](double n,std::array<double,3> b) {
    if(n<=0)return;
    for(std::size_t k=0;k<3;++k)if(active[k]) {
      out.gradient[k]+=constants::R_gas*b[k]*(std::log(n)+1);
      for(std::size_t l=0;l<3;++l)if(active[l])out.hessian[k][l]+=constants::R_gas*b[k]*b[l]/n;
    }
  };
  add(c.X[0],{1,0,0});add(c.X[1]/3,{0,1./3,0});add(c.X[2]/4,{-.25,-.25,-.25});
  for(const auto& m:gs98_metals)if(!potassium(m))add(c.Z()*m.fraction/m.mass_number,{0,0,m.fraction/m.mass_number});
  if(active[1])out.gradient[1]-=.5*constants::R_gas*std::log(nuclides[1].A/nuclides[2].A);
  return out;
}
MetalCompositionHeatResponse VariableMetalHelmholtzEos::composition_heat(double T,double rho,
    const Composition& c,std::array<bool,3> active,bool derivatives) const {
  check_active(c,active);const auto j=jets(T,rho,c,active_any(active)?(derivatives?10:4):1);const auto& f=j[0];
  (void)helmholtz_response(T,rho,f);
  const double denominator=f[0][1]+f[0][2],delta=(f[0][1]+f[1][1])/denominator;
  const double radiation=4*constants::a_rad*T*T*T/(3*rho*denominator);
  if(!(denominator>0) || !std::isfinite(delta+radiation))throw std::domain_error("variable EOS: invalid material heat response");
  MetalCompositionHeatResponse out;out.material_delta=delta;
  const double nan=std::numeric_limits<double>::quiet_NaN();out.exchange_enthalpy.fill(nan);out.radiation_enthalpy.fill(nan);out.delta_partials.fill(nan);
  for(auto& h:out.enthalpy_partials)h.fill(nan);for(auto& h:out.radiation_enthalpy_partials)h.fill(nan);
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    const auto& g=j[1+k];out.exchange_enthalpy[k]=T*(delta*g[0][1]-g[1][0]);out.radiation_enthalpy[k]=T*radiation*g[0][1];
  }
  if(!derivatives)return out;
  auto& d=out.delta_partials;
  std::array<double,5> dr{radiation*(3-(f[1][1]+f[1][2])/denominator),radiation*(-1-(f[0][2]+f[0][3])/denominator),nan,nan,nan};
  d[0]=(f[1][1]+f[2][1]-delta*(f[1][1]+f[1][2]))/denominator;
  d[1]=(f[0][2]+f[1][2]-delta*(f[0][2]+f[0][3]))/denominator;
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    const auto& g=j[1+k];d[2+k]=(g[0][1]+g[1][1]-delta*(g[0][1]+g[0][2]))/denominator;
    dr[2+k]=-radiation*(g[0][1]+g[0][2])/denominator;
    auto& h=out.enthalpy_partials[k];h[0]=out.exchange_enthalpy[k]+T*(d[0]*g[0][1]+delta*g[1][1]-g[2][0]);
    h[1]=T*(d[1]*g[0][1]+delta*g[0][2]-g[1][1]);
    auto& hr=out.radiation_enthalpy_partials[k];hr[0]=out.radiation_enthalpy[k]+T*(dr[0]*g[0][1]+radiation*g[1][1]);
    hr[1]=T*(dr[1]*g[0][1]+radiation*g[0][2]);
  }
  for(std::size_t k=0;k<3;++k)if(active[k])for(std::size_t l=k;l<3;++l)if(active[l]) {
    const auto& second=j[second_channel(k,l)];const double fixed=T*(delta*second[0][1]-second[1][0]);
    out.enthalpy_partials[k][2+l]=fixed+T*d[2+l]*j[1+k][0][1];out.enthalpy_partials[l][2+k]=fixed+T*d[2+k]*j[1+l][0][1];
    const double radfixed=T*radiation*second[0][1];
    out.radiation_enthalpy_partials[k][2+l]=radfixed+T*dr[2+l]*j[1+k][0][1];out.radiation_enthalpy_partials[l][2+k]=radfixed+T*dr[2+k]*j[1+l][0][1];
  }
  for(std::size_t k=0;k<3;++k)if(active[k]) {
    if(!std::isfinite(out.exchange_enthalpy[k]+out.radiation_enthalpy[k]))throw std::domain_error("variable EOS: nonfinite enthalpy");
    for(std::size_t l=0;l<5;++l)if(l<2 || active[l-2])
      if(!std::isfinite(out.enthalpy_partials[k][l]+out.radiation_enthalpy_partials[k][l]+d[l]))
        throw std::domain_error("variable EOS: nonfinite enthalpy derivative");
  }
  return out;
}
} // namespace ember

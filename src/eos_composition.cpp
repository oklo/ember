#include "ember/eos_composition.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace ember {
CompositionHelmholtzEos::CompositionHelmholtzEos(const std::filesystem::path& path,
    HelmholtzTableEos::Mixture mixture) {
  if(mixture!=HelmholtzTableEos::Mixture::allow_documented_proxy)
    throw std::invalid_argument("CompositionHelmholtzEos: explicit metals/isotope approximation required");
  std::ifstream in(path);std::string magic;int version{};std::size_t count{};
  in>>magic>>version>>count;
  if(!in || magic!="EMBER_COMPOSITION_HELMHOLTZ" || version!=1 || count<2 || count>100)
    throw std::runtime_error("CompositionHelmholtzEos: invalid manifest");
  for(std::size_t i=0;i<count;++i) {
    double x{};std::string file;in>>x>>std::quoted(file);
    if(!in || !std::isfinite(x) || x<0 || x>.98 || (i && x<=x_.back()) || file.empty())
      throw std::runtime_error("CompositionHelmholtzEos: invalid composition axis");
    auto table=std::make_unique<HelmholtzTableEos>(path.parent_path()/file,mixture);
    const auto& c=table->composition();
    if(std::abs(c.h1()-x)>1e-12 || c[Species::He3]!=0)
      throw std::runtime_error("CompositionHelmholtzEos: source composition mismatch");
    if(!i) metals_=c;
    for(std::size_t j=3;j<NSPEC;++j) if(std::abs(c.X[j]-metals_.X[j])>1e-12)
      throw std::runtime_error("CompositionHelmholtzEos: metal pattern mismatch");
    x_.push_back(x);tables_.push_back(std::move(table));
  }
  if(in>>magic) throw std::runtime_error("CompositionHelmholtzEos: trailing manifest data");
}

CompositionHelmholtzEos::Coordinates CompositionHelmholtzEos::coordinates(const Composition& c) const {
  if(c.metal_inventory!=MetalInventory::carried_isotopes)
    throw std::domain_error("CompositionHelmholtzEos: GS98 requires the metal-bearing EOS family");
  if(std::abs(c.sum()-1)>1e-10) throw std::domain_error("CompositionHelmholtzEos: abundances must sum to one");
  for(std::size_t i=0;i<NSPEC;++i)
    if(!std::isfinite(c.X[i]) || c.X[i]<0 || (i>=3 && std::abs(c.X[i]-metals_.X[i])>1e-10))
      throw std::domain_error("CompositionHelmholtzEos: invalid abundances or changed metals");
  const double A1=nuclides[0].A,A3=nuclides[1].A,A4=nuclides[2].A;
  const double w1=c.abundance_weight(0),w3=c.abundance_weight(1),w4=c.abundance_weight(2);
  const double h=c.X[0]/w1,n3=c.X[1]/w3,n4=(c.X[2]+c.Z())/w4,he=n3+n4;
  Coordinates q;
  q.scale=A1*h+A4*he;q.x=A1*h/q.scale;
  q.dscale={A1/w1-A4/w4,A4/w3-A4/w4};
  q.dx={(A1/w1-q.x*q.dscale[0])/q.scale,-q.x*q.dscale[1]/q.scale};
  // The reference He4 table has one isotope. Distinguishable He3/He4
  // have ideal mixing entropy and different translational masses. Ignore
  // tiny isotope shifts of electronic levels and nuclear spin constants.
  const auto xlog=[](double n,double total) {return n>0?n*std::log(n/total):0.;};
  q.isotope_phi=constants::kB*constants::NA*(xlog(n3,he)+xlog(n4,he)-1.5*n3*std::log(A3/A4));
  const double tol=2e-14;
  if(q.x<x_.front()-tol || q.x>x_.back()+tol)
    throw std::domain_error("CompositionHelmholtzEos: equivalent hydrogen abundance outside table family");
  q.x=std::clamp(q.x,x_.front(),x_.back());
  q.interval=interp::locate(x_,q.x);
  q.fraction=(q.x-x_[q.interval])/(x_[q.interval+1]-x_[q.interval]);
  return q;
}

std::optional<Eos::DensityRange> CompositionHelmholtzEos::density_range(double T,const Composition& c) const {
  const auto q=coordinates(c);
  const auto a=tables_[q.interval]->material_density_range(T),b=tables_[q.interval+1]->material_density_range(T);
  const DensityRange result{std::max(a.min,b.min)/q.scale,std::min(a.max,b.max)/q.scale};
  if(result.min>=result.max) throw std::domain_error("CompositionHelmholtzEos: empty source overlap");
  return result;
}

EosResponse CompositionHelmholtzEos::eval_with_derivatives(double T,double rho,const Composition& c) const {
  const auto q=coordinates(c);
  const auto a=tables_[q.interval]->material_jet(T,rho*q.scale),b=tables_[q.interval+1]->material_jet(T,rho*q.scale);
  HelmholtzJet f{};
  for(int i=0;i<4;++i) for(int j=0;j<4-i;++j)
    f[i][j]=q.scale*((1-q.fraction)*a[i][j]+q.fraction*b[i][j]);
  f[0][0]+=q.isotope_phi;
  return helmholtz_response(T,rho,f);
}

EosCompositionResponse CompositionHelmholtzEos::composition_response(double T,double rho,const Composition& c) const {
  const auto q=coordinates(c);
  const auto a=tables_[q.interval]->material_jet(T,rho*q.scale),b=tables_[q.interval+1]->material_jet(T,rho*q.scale);
  const double dx=x_[q.interval+1]-x_[q.interval];
  EosCompositionResponse result;
  for(int k=0;k<2;++k) {
    auto partial=[&](int i,int j) {
      return q.dscale[k]*((1-q.fraction)*(a[i][j]+a[i][j+1])+q.fraction*(b[i][j]+b[i][j+1]))
           +q.scale*q.dx[k]*(b[i][j]-a[i][j])/dx;
    };
    result.dP[k]=rho*T*partial(0,1);result.dE[k]=-T*partial(1,0);
  }
  return result;
}
} // namespace ember

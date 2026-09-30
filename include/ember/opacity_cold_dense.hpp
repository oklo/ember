#pragma once
#include "ember/opacity_table.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {

// Optional computed dense-gas table for a nearly metal-free H envelope.
// Input abundances are elemental atomic-mass fractions; place ElementalOpacity
// outside this join. The source remains an approximation to nonideal chemistry.
class ColdDenseOpacity final : public Opacity {
public:
  ColdDenseOpacity(const Opacity& original,const std::filesystem::path& table,double scale=1)
      :original_(original),dense_(table,"cold dense gas",TabulatedOpacity::DensityAxis::logRho,
          AbundanceBasis::baryon_mass),scale_(scale) {
    if(original.includes_conduction() || dense_.metallicity()!=0 || !std::isfinite(scale)
        || scale<.1 || scale>10)
      throw std::invalid_argument("ColdDenseOpacity: radiative zero-metal source and scale .1–10 required");
  }
  OpacityState eval(double T,double rho,const Composition& c) const override {
    if(!(T>0 && rho>0) || !std::isfinite(T+rho))
      throw std::domain_error("ColdDenseOpacity: invalid temperature or density");
    const double Z=c.Z();
    // The computed gas source reaches 20 kK. Returning to the original
    // table at 10–12 kK leaves its log R <= 6 support in cool envelopes.
    // Use the warmer overlap; both sources must still cover every query
    // whenever their weights are nonzero.
    const auto lo=step((T-3000)/500),hi=step((20000-T)/4000);
    const auto x=step((c.X[0]-.98)/.005),z=step((1e-10-Z)/(1e-10-1e-12));
    const double factor=lo.value*hi.value*x.value*z.value;
    if(factor==0)return original_.eval(T,rho,c);
    const double ln10=std::log(10.),logR=std::log10(rho)-3*(std::log10(T)-6);
    const auto d=step((logR-5.6)/.3);
    const double w=factor*d.value;
    if(w==0)return original_.eval(T,rho,c);
    if(c.basis!=AbundanceBasis::atomic_mass || c.X[1]!=0 || c[Species::H2]!=0)
      throw std::domain_error("ColdDenseOpacity: elemental atomic-mass input required");
    // SYNSPEC uses integer H/He masses. Preserve the number densities and
    // extinction per length when converting the elemental input to that basis.
    // The negligible metals are represented by helium at fixed input H.
    constexpr double h=1/nuclides[0].A,he=4/nuclides[2].A,ds=h-he;
    const double s=c.X[0]*h+(1-c.X[0])*he;
    Composition source;source.basis=AbundanceBasis::baryon_mass;
    source.X[0]=c.X[0]*h/s;source.X[2]=1-source.X[0];
    auto b=dense_.eval(T,rho*s,source);
    b.dlnk_dX=b.dlnk_dX*(h-source.X[0]*ds)/s+(1+b.dlnk_dlnRho)*ds/s;
    b.kappa*=scale_*s;b.dlnk_dZ=0;
    if(w==1)return b;
    const auto a=original_.eval(T,rho,c);
    const double delta=std::log(b.kappa/a.kappa);
    const double dwT=factor*d.derivative*(-3/(.3*ln10))
        +d.value*x.value*z.value*T*(lo.derivative*hi.value/500-lo.value*hi.derivative/4000);
    const double dwR=factor*d.derivative/(.3*ln10);
    const double dwX=d.value*lo.value*hi.value*z.value*x.derivative/.005;
    const double dwZ=-d.value*lo.value*hi.value*x.value*z.derivative/(1e-10-1e-12);
    return {std::exp((1-w)*std::log(a.kappa)+w*std::log(b.kappa)),
        (1-w)*a.dlnk_dlnT+w*b.dlnk_dlnT+dwT*delta,
        (1-w)*a.dlnk_dlnRho+w*b.dlnk_dlnRho+dwR*delta,
        (1-w)*a.dlnk_dX+w*b.dlnk_dX+dwX*delta,
        (1-w)*a.dlnk_dZ+dwZ*delta};
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    if(T<=3000 || T>=20000 || c.X[0]<=.98 || c.Z()>=1e-10)
      return original_.density_range(T,c);
    // The union depends on both join weights and each source's stencil.
    // eval enforces both domains wherever their weights are nonzero.
    return std::nullopt;
  }
  const char* name() const override {return "radiative opacity with cold dense-gas source";}
private:
  struct Weight {double value,derivative;};
  static Weight step(double u) {
    if(u<=0)return {0,0};
    if(u>=1)return {1,0};
    return {u*u*u*(10+u*(-15+6*u)),30*u*u*(1-u)*(1-u)};
  }
  const Opacity& original_;
  TabulatedOpacity dense_;
  double scale_;
};
} // namespace ember

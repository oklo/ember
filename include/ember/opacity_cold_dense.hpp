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
    // Return over 17.5–20 kK only where the original log R <= 6 table
    // is supported. At higher density retain the computed source, which
    // still enforces its own temperature and density limits. Every source with
    // nonzero weight must support the query.
    // Below 3500 K the original (AESOPUS) table ends at log R = 6, while cool
    // envelopes reach log R ~ 6.4. Complete the density join to the computed
    // source by log R = 5.98 at every T >= 3000 K (e), so the original is read
    // only where it is supported; unchanged for log R <= 5.9.
    const auto lo=step((T-3000)/500),hi=step((20000-T)/2500);
    const auto x=step((c.X[0]-.98)/.005),z=step((1e-10-Z)/(1e-10-1e-12));
    const double ln10=std::log(10.),logR=std::log10(rho)-3*(std::log10(T)-6);
    const auto d=step((logR-5.6)/.3),e=step((logR-5.9)/.08);
    const double hot=hi.value+(1-hi.value)*e.value;
    const double outer=hot*x.value*z.value;
    if(outer==0)return original_.eval(T,rho,c);
    const double m=lo.value+(1-lo.value)*e.value;
    const double w=outer*d.value*m;
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
    // Derivatives of w = outer d m in ln T (log R moves as -3 log T), ln rho, X and Z.
    const double dd_lnrho=d.derivative/(.3*ln10),de_lnrho=e.derivative/(.08*ln10);
    const double dm_lnT=T*lo.derivative/500*(1-e.value)-3*(1-lo.value)*de_lnrho;
    const double dm_lnrho=(1-lo.value)*de_lnrho;
    const double dh_lnT=-T*hi.derivative/2500*(1-e.value)-3*(1-hi.value)*de_lnrho;
    const double dh_lnrho=(1-hi.value)*de_lnrho;
    const double dwT=outer*(-3*dd_lnrho*m+d.value*dm_lnT)
        +x.value*z.value*d.value*m*dh_lnT;
    const double dwR=outer*(dd_lnrho*m+d.value*dm_lnrho)
        +x.value*z.value*d.value*m*dh_lnrho;
    const double dwX=d.value*m*hot*z.value*x.derivative/.005;
    const double dwZ=-d.value*m*hot*x.value*z.derivative/(1e-10-1e-12);
    return {std::exp((1-w)*std::log(a.kappa)+w*std::log(b.kappa)),
        (1-w)*a.dlnk_dlnT+w*b.dlnk_dlnT+dwT*delta,
        (1-w)*a.dlnk_dlnRho+w*b.dlnk_dlnRho+dwR*delta,
        (1-w)*a.dlnk_dX+w*b.dlnk_dX+dwX*delta,
        (1-w)*a.dlnk_dZ+dwZ*delta};
  }
  std::optional<DensityRange> density_range(double T,const Composition& c) const override {
    if(c.X[0]<=.98 || c.Z()>=1e-10)
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

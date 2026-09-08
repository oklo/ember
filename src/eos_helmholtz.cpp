#include "ember/eos_helmholtz.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include "differential.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
constexpr double ln10 = 2.302585092994045684;
bool positive(double x) { return std::isfinite(x) && x > 0; }
// Quintic endpoint bases: value, first derivative, second derivative,
// at left then right. Derivative slots multiply h and h^2 respectively.
constexpr double basis[6][6] = {
  {1,0,0,-10,15,-6}, {0,1,0,-6,8,-3}, {0,0,.5,-1.5,1.5,-.5},
  {0,0,0,10,-15,6}, {0,0,0,-4,7,-3}, {0,0,0,.5,-1,.5}
};
using Bases = std::array<std::array<double,4>,6>;
Bases evaluate_basis(double u, double h) {
  Bases b{};
  for (int k=0;k<6;++k) for (int d=0;d<4;++d) {
    double value=0;
    for (int n=5;n>=d;--n) {
      double c=basis[k][n];
      for (int j=0;j<d;++j) c*=n-j;
      value=value*u+c;
    }
    b[k][d]=value*std::pow(h,k%3-d);
  }
  return b;
}
}

HelmholtzTableEos::HelmholtzTableEos(const std::filesystem::path& file, Mixture mixture) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("HelmholtzTableEos: cannot open "+file.string());
  auto label=[&](const char* want) {
    std::string key; in>>key;
    if (!in || key!=want) throw std::runtime_error(std::string("HelmholtzTableEos: expected ")+want);
  };
  label("EMBER_HELMHOLTZ"); int version=0; in>>version;
  if (!in || version!=1) throw std::runtime_error("HelmholtzTableEos: unsupported version");
  label("source"); in>>std::quoted(source_);
  if (!in || source_.empty()) throw std::runtime_error("HelmholtzTableEos: missing provenance");
  label("composition_proxy"); std::string proxy; in>>std::quoted(proxy);
  if (!in || proxy.empty()) throw std::runtime_error("HelmholtzTableEos: missing mixture declaration");
  if (proxy!="none" && mixture!=Mixture::allow_documented_proxy)
    throw std::invalid_argument("HelmholtzTableEos: explicit composition proxy selection required");
  label("composition");
  for (double& v:composition_.X) {
    in>>v;
    if (!in || !std::isfinite(v) || v<0 || v>1) throw std::runtime_error("HelmholtzTableEos: bad composition");
  }
  if (std::abs(composition_.sum()-1)>1e-10) throw std::runtime_error("HelmholtzTableEos: composition sum");
  auto axis=[&](const char* key,std::vector<double>& a) {
    label(key); std::size_t n=0; in>>n;
    if (!in || n<2 || n>2000) throw std::runtime_error("HelmholtzTableEos: bad axis size");
    a.resize(n);
    for (std::size_t i=0;i<n;++i) {
      in>>a[i]; a[i]*=ln10;
      if (!in || !std::isfinite(a[i]) || (i && a[i]<=a[i-1]))
        throw std::runtime_error("HelmholtzTableEos: bad axis");
    }
  };
  axis("log_t",t_); axis("log_q",q_);
  if (t_.size()*q_.size()>1000000) throw std::runtime_error("HelmholtzTableEos: too many nodes");
  label("data"); nodes_.resize(t_.size()*q_.size());
  for (auto& n:nodes_) {
    int valid=-1; in>>valid;
    if (!in || (valid!=0 && valid!=1)) throw std::runtime_error("HelmholtzTableEos: invalid mask");
    n.valid=valid==1;
    for (double& v:n.d) {
      in>>v;
      if (!in || !std::isfinite(v)) throw std::runtime_error("HelmholtzTableEos: invalid potential data");
    }
  }
  std::string extra;
  if (in>>extra) throw std::runtime_error("HelmholtzTableEos: trailing data");
}

void HelmholtzTableEos::check_composition(const Composition& c) const {
  for (std::size_t i=0;i<NSPEC;++i)
    if (!std::isfinite(c.X[i]) || c.X[i]<0 || std::abs(c.X[i]-composition_.X[i])>1e-10)
      throw std::domain_error("HelmholtzTableEos: fixed composition only; He3/evolution unsupported");
}

std::pair<std::size_t,std::size_t> HelmholtzTableEos::supported_q(std::size_t it) const {
  // Only the contiguous branch attached to the dilute boundary is exposed.
  // Do not jump across a masked phase/unstable region during PT inversion.
  std::size_t hi=0;
  while (hi<q_.size() && nodes_[it*q_.size()+hi].valid && nodes_[(it+1)*q_.size()+hi].valid) ++hi;
  if (hi<2) throw std::domain_error("HelmholtzTableEos: no supported density interval");
  return {0,hi-1};
}

std::optional<Eos::DensityRange> HelmholtzTableEos::density_range(double T, const Composition& c) const {
  check_composition(c);
  if (!positive(T) || std::log(T)<t_.front() || std::log(T)>t_.back())
    throw std::domain_error("HelmholtzTableEos: temperature outside table");
  const double t=std::log(T);
  const auto [lo,hi]=supported_q(interp::locate(t_,t));
  const double offset=1.5*(t-6*ln10);
  // Inset by roundoff so transformation back to q cannot leave the grid.
  return DensityRange{std::exp(q_[lo]+offset+1e-12),std::exp(q_[hi]+offset-1e-12)};
}

EosResponse HelmholtzTableEos::eval_with_derivatives(double T, double rho, const Composition& c) const {
  check_composition(c);
  if (!positive(T) || !positive(rho)) throw std::domain_error("HelmholtzTableEos: invalid state");
  const double t=std::log(T), q=std::log(rho)-1.5*(t-6*ln10);
  if (t<t_.front() || t>t_.back() || q<q_.front() || q>q_.back())
    throw std::domain_error("HelmholtzTableEos: state outside table");
  const std::size_t it=interp::locate(t_,t), iq=interp::locate(q_,q);
  const auto [lo,hi]=supported_q(it);
  if (iq<lo || iq>=hi) throw std::domain_error("HelmholtzTableEos: masked density region");
  const double ht=t_[it+1]-t_[it], hq=q_[iq+1]-q_[iq];
  const auto bt=evaluate_basis((t-t_[it])/ht,ht), bq=evaluate_basis((q-q_[iq])/hq,hq);
  double f[4][4]{};
  for (std::size_t si=0;si<2;++si) for (std::size_t sj=0;sj<2;++sj) {
    const auto& node=nodes_[(it+si)*q_.size()+iq+sj];
    for (std::size_t i=0;i<3;++i) for (std::size_t j=0;j<3;++j)
      for (std::size_t a=0;a<4;++a) for (std::size_t b=0;b<4-a;++b)
        f[a][b]+=node.d[3*i+j]*bt[3*si+i][a]*bq[3*sj+j][b];
  }
  // Transform x derivatives at fixed Q to t derivatives at fixed rho.
  // D_t = D_x - 1.5 D_q; D_r = D_q. Add phi_rad analytically.
  const double rad=-constants::a_rad*std::pow(T,3)/(3*rho);
  const double phi=f[0][0]+rad;
  const double ft=f[1][0]-1.5*f[0][1]+3*rad, fr=f[0][1]-rad;
  const double ftt=f[2][0]-3*f[1][1]+2.25*f[0][2]+9*rad;
  const double ftr=f[1][1]-1.5*f[0][2]-3*rad, frr=f[0][2]+rad;
  const double fttt=f[3][0]-4.5*f[2][1]+6.75*f[1][2]-3.375*f[0][3]+27*rad;
  const double fttr=f[2][1]-3*f[1][2]+2.25*f[0][3]-9*rad;
  const double ftrr=f[1][2]-1.5*f[0][3]+3*rad, frrr=f[0][3]-rad;
  using D=detail::Differential<2>;
  auto dual=[](double v,double dt,double dr) { D a(v); a.d={dt,dr}; return a; };
  const auto C=dual(fr,ftr,frr), cr=1+dual(frr,ftrr,frrr)/C;
  const auto ct=1+dual(ftr,fttr,ftrr)/C;
  const auto cv=dual(-ft-ftt,-ftt-fttt,-ftr-fttr);
  const auto delta=ct/cr, cp=cv+C*ct*delta, ad=C*delta/cp;
  EosResponse a{}; auto& e=a.state;
  e.P=rho*T*fr; e.E=-T*ft; e.S=-phi-ft;
  e.chiT=ct.value; e.chiRho=cr.value; e.cv=cv.value; e.cp=cp.value;
  e.delta=delta.value; e.grad_ad=ad.value; e.Gamma1=(cr*cp/cv).value;
  // Particle fractions are not derivatives of the equilibrium potential.
  e.mu=e.free_e=std::numeric_limits<double>::quiet_NaN();
  if (!positive(e.P) || !positive(e.cv) || !positive(e.cp) || !positive(e.chiRho)
      || !positive(e.delta) || !positive(e.grad_ad) || !std::isfinite(e.E) || !std::isfinite(e.S))
    throw std::domain_error("HelmholtzTableEos: unstable potential interpolant");
  a.dE_dlnRho=-T*ftr;
  a.dcp_dlnT=cp.d[0]; a.dcp_dlnRho=cp.d[1];
  a.ddelta_dlnT=delta.d[0]; a.ddelta_dlnRho=delta.d[1];
  a.dgrad_ad_dlnT=ad.d[0]; a.dgrad_ad_dlnRho=ad.d[1];
  return a;
}

EosState HelmholtzTableEos::eval(double T, double rho, const Composition& c) const {
  return eval_with_derivatives(T,rho,c).state;
}

double HelmholtzTableEos::rho_from_PT(double T,double P,const Composition& c,double guess) const {
  if (!positive(P) || !std::isfinite(guess) || guess<0)
    throw std::domain_error("HelmholtzTableEos: invalid pressure or density guess");
  const auto range=density_range(T,c).value();
  double lo=std::log(range.min), hi=std::log(range.max);
  const double target=std::log(P);
  if (std::log(eval(T,range.min,c).P)>target || std::log(eval(T,range.max,c).P)<target)
    throw std::domain_error("HelmholtzTableEos: pressure outside supported density interval");
  double r=positive(guess)?std::clamp(std::log(guess),lo,hi):(lo+hi)/2;
  for (int i=0;i<100;++i) {
    const auto e=eval(T,std::exp(r),c);
    const double f=std::abs(e.P-P)<.5*P ? std::log1p((e.P-P)/P) : std::log(e.P)-target;
    if (std::abs(f/e.chiRho)<2e-13) return std::exp(r);
    if (f<0) lo=r; else hi=r;
    const double next=r-f/e.chiRho;
    r=next>lo && next<hi?next:(lo+hi)/2;
  }
  throw std::runtime_error("HelmholtzTableEos: PT inversion did not converge");
}
} // namespace ember

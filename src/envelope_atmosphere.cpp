#include "ember/envelope_atmosphere.hpp"
#include "ember/constants.hpp"
#include "ember/convection.hpp"
#include "ember/eos.hpp"
#include "parallel_evaluate.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace ember {
namespace {
void require(bool ok, const char* why) { if (!ok) throw std::domain_error(std::string("envelope atmosphere: ") + why); }
// Catmull-Rom weights for a point at fractional position u in [0,1] between nodes 1 and 2 of a 4-node stencil.
std::array<double,4> cr(double u) {
  const double u2 = u * u, u3 = u2 * u;
  return {-.5 * u3 + u2 - .5 * u, 1.5 * u3 - 2.5 * u2 + 1, -1.5 * u3 + 2 * u2 + .5 * u, .5 * u3 - .5 * u2};
}
std::size_t cell(const std::vector<double>& x, double v, const char* what) {
  require(v >= x[1] && v <= x[x.size() - 2], what);   // a full 4-point stencil is required: no edge extrapolation
  const auto i = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), v) - x.begin()) - 1;
  return std::min(i, x.size() - 3);
}
}  // namespace

EnvelopeSource::EnvelopeSource(const std::string& path) {
  std::ifstream f(path); std::string tag; f >> tag;
  require(tag == "EMBER_ENVELOPE_SOURCE_V1", "unsupported source format");
  std::size_t nx, ny, nt, nr; f >> nx >> ny >> nt >> nr;
  require(f.good() && nx >= 2 && ny >= 2 && nt >= 4 && nr >= 4, "invalid source dimensions");
  auto read = [&](std::vector<double>& a, std::size_t n) { a.resize(n); for (auto& v : a) { f >> v; require(f.good() && std::isfinite(v), "invalid source value"); } };
  read(X_, nx); read(Y3_, ny); read(lt_, nt); read(lr_, nr);
  planes_.resize(nx * ny);
  for (auto& p : planes_) for (auto& q : p) read(q, nt * nr);
  for (const auto* a : {&X_, &Y3_, &lt_, &lr_}) for (std::size_t i = 1; i < a->size(); ++i) require((*a)[i] > (*a)[i - 1], "unordered axis");
  // Catmull-Rom weights assume uniform spacing in ln T and ln rho: verify rather than assume.
  for (const auto* a : {&lt_, &lr_}) {
    const double h = ((*a).back() - (*a).front()) / static_cast<double>(a->size() - 1);
    for (std::size_t i = 1; i < a->size(); ++i) require(std::abs((*a)[i] - (*a)[i - 1] - h) <= 1e-9 * std::abs(h), "nonuniform ln T or ln rho axis");
  }
}
EnvelopeSource::State EnvelopeSource::eval(double lnT, double lnrho, double X, double Y3) const {
  require(X >= X_.front() && X <= X_.back() && Y3 >= Y3_.front() && Y3 <= Y3_.back(), "composition outside declared source interval");
  const auto ix = std::min<std::size_t>(std::upper_bound(X_.begin(), X_.end(), X) - X_.begin() - 1, X_.size() - 2);
  const auto iy = std::min<std::size_t>(std::upper_bound(Y3_.begin(), Y3_.end(), Y3) - Y3_.begin() - 1, Y3_.size() - 2);
  const double fx = (X - X_[ix]) / (X_[ix + 1] - X_[ix]), fy = (Y3 - Y3_[iy]) / (Y3_[iy + 1] - Y3_[iy]);
  const auto it=cell(lt_,lnT,"temperature outside source"),ir=cell(lr_,lnrho,"density outside source");
  const auto wt=cr((lnT-lt_[it])/(lt_[it+1]-lt_[it])),wr=cr((lnrho-lr_[ir])/(lr_[ir+1]-lr_[ir]));
  const std::size_t nr=lr_.size();
  std::array<double,5> v{};
  for (int dx = 0; dx < 2; ++dx) for (int dy = 0; dy < 2; ++dy) {
    const double w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy); if (w == 0) continue;
    const std::size_t plane = (ix + dx) * Y3_.size() + (iy + dy);
    // Explicit support: never interpolate across a nonpositive source response.
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)for(std::size_t q=1;q<5;++q)
      require(planes_[plane][q][(it-1+a)*lr_.size()+ir-1+b]>0,
              "unsupported envelope-source response stencil");
    for (std::size_t q=0;q<5;++q) {
      const auto& values=planes_[plane][q];double value=0;
      for(int a=0;a<4;++a)for(int b=0;b<4;++b)
        value+=wt[a]*wr[b]*values[(it-1+a)*nr+ir-1+b];
      v[q]+=w*value;
    }
  }
  return {v[0], v[1], v[2], v[3], v[4]};
}
double EnvelopeSource::lnrho_from(double lnT, double lnP, double X, double Y3, double guess) const {
  double x = guess;
  for (int it = 0; it < 60; ++it) {
    const auto s = eval(lnT, x, X, Y3);
    require(s.chiRho > 0, "non-positive chiRho in source");
    const double dx = -(s.lnP - lnP) / s.chiRho; x += std::clamp(dx, -1.0, 1.0);
    if (std::abs(dx) < 1e-12) return x;
  }
  throw std::domain_error("envelope atmosphere: source density inversion did not converge");
}

EnvelopeAtmosphere::EnvelopeAtmosphere(const Atmosphere& top, const Opacity& opacity, const EnvelopeSource& source,
    double alpha_mlt, double total_mass, double envelope_mass, double steps_per_unit_lnP, EnvelopeMetals metals)
    : top_(top), opacity_(opacity), source_(&source), metals_(metals), alpha_(alpha_mlt), M_(total_mass), dM_(envelope_mass), per_unit_(steps_per_unit_lnP) {
  require(alpha_ > 0 && M_ > 0 && dM_ > 0 && dM_ < 0.05 * M_ && per_unit_ >= 10, "invalid envelope parameters");
}

EnvelopeAtmosphere::EnvelopeAtmosphere(const Atmosphere& top, const Opacity& opacity, const Eos& eos,
    double alpha_mlt,double total_mass,double envelope_mass,double steps_per_unit_lnP)
    :top_(top),opacity_(opacity),eos_(&eos),alpha_(alpha_mlt),M_(total_mass),dM_(envelope_mass),
     per_unit_(steps_per_unit_lnP) {
  require(alpha_>0 && M_>0 && dM_>0 && dM_<.05*M_ && per_unit_>=10,"invalid envelope parameters");
}

EnvelopeSource::PressureState EnvelopeAtmosphere::at_pressure(double lnT,double lnP,
    const Composition& c,double& guess) const {
  if(source_)return source_->at_pressure(lnT,lnP,c,guess,metals_);
  const double T=std::exp(lnT),P=std::exp(lnP);
  const double rho=eos_->rho_from_PT(T,P,c,std::exp(guess));guess=std::log(rho);
  const auto state=eos_->eval(T,rho,c);
  return {rho,state.cp,state.delta,state.grad_ad,state.chiRho};
}

EnvelopeSource::PressureState EnvelopeSource::at_pressure(double lnT, double lnP,
    const Composition& c, double& guess, EnvelopeMetals metals) const {
  using namespace constants;
  const double Z=c.Z(), D=c[Species::H2];
  require(std::isfinite(c.sum()) && std::abs(c.sum()-1)<1e-10
      && Z>=0 && Z<=.04 && D>=0 && D<=1e-4, "unsupported layer composition");
  for(double x:c.X)require(std::isfinite(x) && x>=0, "invalid layer abundance");
  require(metals!=EnvelopeMetals::reject || Z<=EnvelopeAtmosphere::max_Z,
      "metal fraction too large for the H/He-only layer source");
  // Preserve the previously selected trace-metal approximation in reject mode.
  // D uses the same total-H mass proxy as the top atmosphere; fuel is untouched.
  const double weight=metals==EnvelopeMetals::reject?1:1-Z;
  guess=lnrho_from(lnT,lnP,(c.X[0]+D)/weight,c.X[1]/weight,guess);
  const auto hh=eval(lnT,guess,(c.X[0]+D)/weight,c.X[1]/weight);
  const double rh=std::exp(guess);
  if(metals==EnvelopeMetals::reject || Z==0)return {rh,hh.cp,hh.delta,hh.grad_ad,hh.chiRho};
  const double T=std::exp(lnT),P=std::exp(lnP),Pr=a_rad*std::pow(T,4)/3,Pg=P-Pr;
  require(Pg>0 && hh.cp>0 && hh.delta>0 && hh.chiRho>0, "invalid layer thermodynamics");
  double particles=0;
  if(c.metal_inventory==MetalInventory::gs98)
    particles=c.metal_ion_moment(0)+(metals==EnvelopeMetals::ionized?c.metal_ion_moment(1):0);
  else for(std::size_t k=METAL_BEGIN;k<METAL_END;++k)
    particles+=c.X[k]/Z/c.abundance_weight(k)*(1+(metals==EnvelopeMetals::ionized?nuclides[k].Z:0));
  const double Rz=R_gas*particles, vz=Rz*T/Pg;
  const double dz=1+4*Pr/Pg, cz=2.5*Rz+4*Pr*vz/T*(4+dz);
  // Add Gibbs free energies at common TOTAL P,T. Each component includes
  // radiation, so their mass-weighted volume carries radiation exactly once.
  const double vh=weight/rh, vm=Z*vz, v=vh+vm;
  const double cp=weight*hh.cp+Z*cz, delta=(vh*hh.delta+vm*dz)/v;
  return {1/v,cp,delta,P*v*delta/(T*cp),v/(vh/hh.chiRho+vm*P/Pg)};
}

double EnvelopeSourceDensity::rho_from_PT(double T,double P,const Composition& c,double rho_guess) const {
  require(std::isfinite(T) && T>0 && std::isfinite(P) && P>0,"invalid density query");
  const double gasP=P-constants::a_rad*std::pow(T,4)/3;
  require(gasP>0,"nonpositive gas pressure in density query");
  // A neutral atomic gas provides a starting guess only. The returned value
  // satisfies the source pressure inversion, including its molecular physics.
  if(!(std::isfinite(rho_guess) && rho_guess>0))
    rho_guess=gasP/(constants::R_gas*T*c.mu_ions_inv());
  const auto range=source_.lnrho_interval();
  double guess=std::clamp(std::log(rho_guess),range[0],range[1]);
  return source_.at_pressure(std::log(T),std::log(P),c,guess,metals_).rho;
}

EnvelopeSurface EnvelopeAtmosphere::integrate(double Teff, double R, const Composition& c) const {
  using namespace constants;
  const double L = 4 * M_PI * R * R * sigma_SB * std::pow(Teff, 4), g0 = G * M_ / (R * R);
  const auto a = [&] {
    try{return top_.eval(Teff,g0,c);}
    catch(const std::exception& error) {
      throw std::domain_error("envelope surface Teff="+std::to_string(Teff)
          +" logg="+std::to_string(std::log10(g0))+": "+error.what());
    }
  }();
  // Integrate mass measured inward from the surface. Subtracting tiny shell
  // masses repeatedly from M loses precision in the final mass-depth match.
  const double target = dM_;
  double lnP = std::log(a.P), lnT = std::log(a.T), r = R, lnrho = std::log(a.rho);
  EnvelopeSurface out{R, Teff, L, a.T, a.P};
  auto rhs = [&](double lp, const std::array<double,3>& y, double& guess) {
    const double T = std::exp(y[0]), P = std::exp(lp);
    const auto s = at_pressure(y[0], lp, c, guess);
    const double rho = s.rho;
    const double enclosed_mass=M_-y[2];
    const double kappa = opacity_.eval(T, rho, c).kappa, g = G * enclosed_mass / (y[1] * y[1]);
    const double grad_rad = 3 * kappa * L * P / (16 * M_PI * a_rad * c_light * G * enclosed_mass * T * T * T * T);
    EosState st{}; st.P = P; st.cp = s.cp; st.delta = s.delta;
    const double U = mixing_length_U(T, rho, kappa, g, st, alpha_);
    const double grad = ledoux_mixing_length_gradient(grad_rad, s.grad_ad, 0., U).grad;
    const double drdlp = -P / (rho * g);
    return std::array<double,3>{grad, drdlp, -4 * M_PI * y[1] * y[1] * rho * drdlp};
  };
  auto step = [&](double lp, const std::array<double,3>& y, double h, double& guess) {
    auto add = [](const std::array<double,3>& u, const std::array<double,3>& k, double s) { return std::array<double,3>{u[0] + s * k[0], u[1] + s * k[1], u[2] + s * k[2]}; };
    const auto k1 = rhs(lp, y, guess), k2 = rhs(lp + h / 2, add(y, k1, h / 2), guess);
    const auto k3 = rhs(lp + h / 2, add(y, k2, h / 2), guess), k4 = rhs(lp + h, add(y, k3, h), guess);
    std::array<double,3> n{}; for (int i = 0; i < 3; ++i) n[i] = y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return n;
  };
  ++integrations_;
  const double h = 1.0 / per_unit_; std::array<double,3> y{lnT, r, 0};
  for (std::size_t n = 0; n < 200000; ++n) {
    double guess = lnrho; auto yn = step(lnP, y, h, guess);
    if (yn[2] >= target) {   // shrink the last step so the base lies exactly at M - dM
      double lo = 0, hi = h, gg = lnrho; std::array<double,3> ym{};
      for (int it = 0; it < 60; ++it) { const double hm = .5 * (lo + hi); gg = lnrho; ym = step(lnP, y, hm, gg); (ym[2] < target ? lo : hi) = hm; if (hi - lo < 1e-14) break; }
      out.P_base = std::exp(lnP + .5 * (lo + hi)); out.T_base = std::exp(ym[0]); out.r_base = ym[1];
      out.rho_base = at_pressure(ym[0], std::log(out.P_base), c, gg).rho;
      out.steps = n + 1; return out;
    }
    y = yn; lnP += h; lnrho = guess;
  }
  throw std::domain_error("envelope atmosphere: step limit before reaching the envelope base");
}

EnvelopeSurface EnvelopeAtmosphere::solve(double L, double r_b, const Composition& c) const {
  using namespace constants;
  auto teff = [&](double R) { return std::pow(L / (4 * M_PI * sigma_SB * R * R), .25); };
  EnvelopeSurface last{};
  auto f = [&](double R) { last = integrate(teff(R), R, c); return last.r_base - r_b; };
  const double q0 = ratio_.load();
  if (q0 > 1 && std::isfinite(q0)) {
    auto usable = [&](double radius) { return std::isfinite(radius) && radius > r_b * (1 + 1e-12); };
    try {
      double x0 = r_b * q0;
      if (usable(x0)) {
        double f0 = f(x0);
        if (std::abs(f0) < 1e-11 * r_b) { ratio_.store(x0 / r_b); return last; }
        double x1 = x0 - f0 * q0;
        for (int it = 0; it < 8 && usable(x1); ++it) {
          const double f1 = f(x1);
          if (std::abs(f1) < 1e-11 * r_b) { ratio_.store(x1 / r_b); return last; }
          if (!std::isfinite(f1) || f1 == f0) break;
          const double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
          x0 = x1; f0 = f1; x1 = x2;
        }
      }
    } catch (const std::domain_error&) {
      // An unsupported predicted radius falls back to the existing bracket.
    }
  }
  // r_base increases with R. Bracket around the last solution's R/r_b when available, widening geometrically.
  const double q = ratio_.load();
  double a, b, fa, fb;
  if (q > 1) {
    double w = 2e-6; a = r_b * (q - w); b = r_b * (q + w); fa = f(a); fb = f(b);
    for (int k = 0; (fa > 0 || fb < 0) && k < 40; ++k) {
      w *= 8;
      if (fa > 0) { b = a; fb = fa; a = r_b * std::max(1 + 1e-9, q - w); fa = f(a); }
      else { a = b; fa = fb; b = r_b * (q + w); fb = f(b); }
    }
  } else {
    // On restart the photospheric radius is unknown. r_b can correspond to
    // a temperature outside a narrow atmosphere table even when the physical
    // surface is supported. Find supported trials before asking for a bracket.
    // Failed source queries are never used as boundary values.
    std::string last_error;
    auto trial=[&](double radius)->std::optional<double> {
      try{return f(radius);}
      catch(const std::domain_error& error){last_error=error.what();return {};}
    };
    std::optional<std::pair<double,double>> low;
    bool bracketed=false;double previous=r_b;
    for(int k=0;k<512 && !bracketed;++k) {
      const double radius=r_b*std::exp(1e-7+.005*k);
      const auto value=trial(radius);
      if(value) {
        if(*value<=0)low=std::pair{radius,*value};
        else {
          b=radius;fb=*value;
          if(low){a=low->first;fa=low->second;bracketed=true;break;}
          // The first supported trial may lie just above the root. Approach
          // the domain edge from its supported side rather than extrapolate.
          double outside=previous,inside=radius;
          for(int j=0;j<40;++j) {
            const double mid=.5*(outside+inside);const auto v=trial(mid);
            if(v && *v<=0){a=mid;fa=*v;bracketed=true;break;}
            if(v){inside=mid;b=mid;fb=*v;}else outside=mid;
          }
          if(!bracketed)break;
        }
      }else if(low) {
        double inside=low->first,outside=radius;
        for(int j=0;j<40;++j) {
          const double mid=.5*(inside+outside);const auto v=trial(mid);
          if(v && *v>=0){a=low->first;fa=low->second;b=mid;fb=*v;bracketed=true;break;}
          if(v){inside=mid;low=std::pair{mid,*v};}else outside=mid;
        }
        if(!bracketed)low.reset();
      }
      previous=radius;
    }
    if(!bracketed)throw std::domain_error("envelope atmosphere: no supported photospheric-radius bracket; "+last_error);
  }
  require(fa <= 0 && fb >= 0, "photospheric radius not bracketed");
  for (int it = 0; it < 100; ++it) {   // Illinois regula falsi
    const double c1 = b - fb * (b - a) / (fb - fa), fc = f(c1);
    if (std::abs(fc) < 1e-11 * r_b) { ratio_.store(c1 / r_b); return last; }
    if (fc * fb < 0) { a = b; fa = fb; } else fa *= .5;
    b = c1; fb = fc;
  }
  throw std::domain_error("envelope atmosphere: photospheric radius did not converge");
}

EnvelopeSurface EnvelopeAtmosphere::photosphere(double Teff_b, double g_b, const Composition& c) const {
  require(eos_ || metals_!=EnvelopeMetals::reject || c.Z() <= max_Z, "metal fraction too large for the H/He-only layer source");
  using namespace constants;
  const double r_b = std::sqrt(G * M_ / g_b), L = 4 * M_PI * r_b * r_b * sigma_SB * std::pow(Teff_b, 4);
  return solve(L, r_b, c);
}

void EnvelopeAtmosphere::evaluation_threads(std::size_t threads) {
  if(threads==0 || threads>64)throw std::invalid_argument("envelope evaluation threads must be between 1 and 64");
  threads_=threads;
}

void EnvelopeAtmosphere::jacobian_reuse(double radius) {
  if(!std::isfinite(radius) || radius<0 || radius>.01)
    throw std::invalid_argument("envelope Jacobian reuse radius must lie in [0,0.01]");
  std::lock_guard lock(jacobian_mutex_);
  jacobian_radius_=radius;
  jacobian_.valid=false;
}

AtmosphereState EnvelopeAtmosphere::eval_value(double Teff_b, double g_b, const Composition& c) const {
  const auto s=photosphere(Teff_b,g_b,c);
  AtmosphereState out{};
  out.T=s.T_base;out.P=s.P_base;out.rho=s.rho_base;
  out.Pgas=s.P_base-constants::a_rad*std::pow(s.T_base,4)/3;
  out.tau=tau_sentinel;
  return out;
}

AtmosphereState EnvelopeAtmosphere::eval(double Teff_b, double g_b, const Composition& c) const {
  require(eos_ || metals_!=EnvelopeMetals::reject || c.Z() <= max_Z, "metal fraction too large for the H/He-only layer source");
  const auto s = photosphere(Teff_b, g_b, c);
  AtmosphereState out{};
  out.T=s.T_base;out.P=s.P_base;out.rho=s.rho_base;
  out.Pgas=s.P_base-constants::a_rad*std::pow(s.T_base,4)/3;
  out.tau=tau_sentinel;
  const double logTeff=std::log(Teff_b),logg=std::log(g_b);
  if(jacobian_radius_>0) {
    std::lock_guard lock(jacobian_mutex_);
    bool nearby=jacobian_.valid && jacobian_.uses<7
      && std::abs(logTeff-jacobian_.logTeff)<=jacobian_radius_
      && std::abs(logg-jacobian_.logg)<=jacobian_radius_
      && c.basis==jacobian_.composition.basis;
    for(std::size_t i=0;i<c.X.size();++i)
      nearby=nearby && std::abs(c.X[i]-jacobian_.composition.X[i])<=.01*jacobian_radius_;
    if(nearby) {
      ++jacobian_.uses;++jacobian_reused_;
      out.dlnT_dlnTeff=jacobian_.derivatives[0];out.dlnT_dlng=jacobian_.derivatives[1];
      out.dlnP_dlnTeff=jacobian_.derivatives[2];out.dlnP_dlng=jacobian_.derivatives[3];
      return out;
    }
  }
  constexpr double h = 1e-4;   // central differences of the full map; checked against h/2 in the controls
  const std::array<std::array<double,2>,4> offsets{{{h,0},{-h,0},{0,h},{0,-h}}};
  std::array<std::array<double,2>,4> values{};
  std::array<std::exception_ptr,4> failures{};
  const auto workers=detail::inside_parallel_evaluation?std::size_t{1}:std::min<std::size_t>(threads_,4);
  auto work=[&](std::size_t worker) {
    // Each derivative starts from the same radius guess and owns its root
    // state. Thread scheduling cannot change another derivative's guess.
    auto evaluate=[&](EnvelopeAtmosphere& local) {
      for(std::size_t i=worker;i<offsets.size();i+=workers) {
        try {
          local.ratio_.store(s.R/s.r_base);
          const auto e=local.photosphere(Teff_b*std::exp(offsets[i][0]),g_b*std::exp(offsets[i][1]),c);
          values[i]={std::log(e.T_base),std::log(e.P_base)};
        }catch(...) {failures[i]=std::current_exception();}
      }
      integrations_.fetch_add(local.integrations());
    };
    if(eos_) {EnvelopeAtmosphere local(top_,opacity_,*eos_,alpha_,M_,dM_,per_unit_);evaluate(local);}
    else {EnvelopeAtmosphere local(top_,opacity_,*source_,alpha_,M_,dM_,per_unit_,metals_);evaluate(local);}
  };
  if(workers==1)work(0);
  else detail::evaluation_workers().run(workers,[&](std::size_t worker) {
    detail::EvaluationScope scope;work(worker);
  });
  for(const auto& error:failures)if(error)std::rethrow_exception(error);
  const auto& tp=values[0];const auto& tm=values[1];const auto& gp=values[2];const auto& gm=values[3];
  out.dlnT_dlnTeff = (tp[0] - tm[0]) / (2 * h); out.dlnT_dlng = (gp[0] - gm[0]) / (2 * h);
  out.dlnP_dlnTeff = (tp[1] - tm[1]) / (2 * h); out.dlnP_dlng = (gp[1] - gm[1]) / (2 * h);
  ++jacobian_computed_;
  if(jacobian_radius_>0) {
    std::lock_guard lock(jacobian_mutex_);
    jacobian_={true,logTeff,logg,c,
      {out.dlnT_dlnTeff,out.dlnT_dlng,out.dlnP_dlnTeff,out.dlnP_dlng},0};
  }
  return out;
}

}  // namespace ember

namespace ember {
namespace {
void need(bool ok, const char* why) { if (!ok) throw std::domain_error(std::string("envelope map: ") + why); }
}
EnvelopeMapAtmosphere::EnvelopeMapAtmosphere(const std::string& path) {
  std::ifstream f(path); std::string tag; f >> tag; need(tag == "EMBER_ENVELOPE_MAP_V1", "unsupported format");
  std::size_t nx, ny, nt, ng; f >> M_ >> dM_ >> nx >> ny >> nt >> ng;
  need(f.good() && nx >= 1 && nx <= 100 && ny >= 1 && ny <= 100 && nt >= 4 && nt <= 1000
      && ng >= 4 && ng <= 1000 && nx*ny*nt*ng <= 1000000, "invalid dimensions");
  need(std::isfinite(M_) && std::isfinite(dM_) && dM_>0 && M_>dM_, "invalid envelope mass");
  auto read = [&](std::vector<double>& a, std::size_t n) { a.resize(n); for (auto& v : a) { f >> v; need(f.good() && std::isfinite(v), "invalid value"); } };
  read(X_, nx); read(Y3_, ny); read(lte_, nt); read(lg_, ng);
  for(const auto* a:{&X_,&Y3_,&lte_,&lg_})
    need(std::adjacent_find(a->begin(),a->end(),std::greater_equal<double>())==a->end(),"non-increasing axis");
  for (const auto* a : {&lte_, &lg_}) {
    const double h = (a->back() - a->front()) / static_cast<double>(a->size() - 1);
    for (std::size_t i = 1; i < a->size(); ++i) need(std::abs((*a)[i] - (*a)[i - 1] - h) <= 1e-9 * std::abs(h), "nonuniform axis");
  }
  planes_.resize(nx * ny);
  for (auto& p : planes_) {
    for (auto& q : p) {
      read(q, nt * ng);
      for(auto& v:q)if(v==1e300)v=std::numeric_limits<double>::quiet_NaN();
    }
    for(std::size_t i=0;i<nt*ng;++i)
      for(const auto& q:p)need(std::isnan(q[i])==std::isnan(p[0][i]),"partially masked node");
  }
}
std::array<double,9> EnvelopeMapAtmosphere::plane_hermite(std::size_t plane, double x, double y) const {
  const std::size_t nt = lte_.size(), ng = lg_.size();
  need(x >= lte_.front() && x <= lte_.back() && y >= lg_.front() && y <= lg_.back(), "(Teff_b, g_b) outside the tabulated corridor");
  const double hx = lte_[1] - lte_[0], hy = lg_[1] - lg_[0];
  const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>((x - lte_[0]) / hx), nt - 2);
  const std::size_t j = std::min<std::size_t>(static_cast<std::size_t>((y - lg_[0]) / hy), ng - 2);
  const double u = (x - lte_[i]) / hx, v = (y - lg_[j]) / hy;
  const auto& P = planes_[plane];
  auto at = [&](std::size_t q, std::size_t a, std::size_t b) { return P[q][a * ng + b]; };
  // cross derivative d2f/dxdy at node (a,b) from centred (one-sided at edges) differences of the stored x-slope along y
  auto cross = [&](std::size_t sx, std::size_t a, std::size_t b) {
    std::size_t b0 = b == 0 ? 0 : b - 1, b1 = b + 1 == ng ? b : b + 1;
    if(std::isnan(at(sx,a,b0)))b0=b;
    if(std::isnan(at(sx,a,b1)))b1=b;
    return b1==b0?0.:(at(sx, a, b1) - at(sx, a, b0)) / ((static_cast<double>(b1) - static_cast<double>(b0)) * hy);
  };
  for(int a=0;a<2;++a)for(int b=0;b<2;++b)
    need(!std::isnan(at(0,i+a,j+b)),"map cell touches an unsupported node");
  auto h0 = [](double t) { return (1 + 2 * t) * (1 - t) * (1 - t); };
  auto h1 = [](double t) { return t * (1 - t) * (1 - t); };
  auto d0 = [](double t) { return 6 * t * t - 6 * t; };           // d/dt of h0
  auto d1 = [](double t) { return (1 - t) * (1 - 3 * t); };        // d/dt of h1
  std::array<double,9> out{};
  for (int k = 0; k < 3; ++k) {                                      // lnT, lnP, lnR
    const std::size_t sx = 3 + 2 * k, sy = 4 + 2 * k;
    double f = 0, fx = 0, fy = 0;
    for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) {
      const double ta = a ? 1 - u : u, tb = b ? 1 - v : v, sa = a ? -1 : 1, sb = b ? -1 : 1;
      const double F = at(k, i + a, j + b), Fx = at(sx, i + a, j + b) * hx, Fy = at(sy, i + a, j + b) * hy, Fxy = cross(sx, i + a, j + b) * hx * hy;
      const double wx = h0(ta), wxd = h1(ta) * sa, wy = h0(tb), wyd = h1(tb) * sb;
      const double dwx = d0(ta) * sa, dwxd = d1(ta), dwy = d0(tb) * sb, dwyd = d1(tb);
      f += F * wx * wy + Fx * wxd * wy + Fy * wx * wyd + Fxy * wxd * wyd;
      fx += (F * dwx * wy + Fx * dwxd * wy + Fy * dwx * wyd + Fxy * dwxd * wyd) / hx;
      fy += (F * wx * dwy + Fx * wxd * dwy + Fy * wx * dwyd + Fxy * wxd * dwyd) / hy;
    }
    out[k] = f; out[3 + 2 * k] = fx; out[4 + 2 * k] = fy;
  }
  return out;
}
EnvelopeMapAtmosphere::Values EnvelopeMapAtmosphere::values(double Teff_b, double g_b, const Composition& c) const {
  need(c.Z() <= EnvelopeAtmosphere::max_Z, "metal fraction too large for the H/He-only layer source");
  const double X = c.X[0], Y3 = c.X[1], x = std::log(Teff_b), y = std::log(g_b);
  need(X >= X_.front() && X <= X_.back() && Y3 >= Y3_.front() && Y3 <= Y3_.back(), "composition outside declared interval");
  auto bracket = [](const std::vector<double>& a, double v, std::size_t& i, double& w) {
    if (a.size() == 1) { i = 0; w = 0; return; }
    i = std::min<std::size_t>(std::upper_bound(a.begin(), a.end(), v) - a.begin() - 1, a.size() - 2);
    w = (v - a[i]) / (a[i + 1] - a[i]);
  };
  std::size_t ix, iy; double fx, fy; bracket(X_, X, ix, fx); bracket(Y3_, Y3, iy, fy);
  std::array<double,9> s{};
  for (int dx = 0; dx < (X_.size() > 1 ? 2 : 1); ++dx) for (int dy = 0; dy < (Y3_.size() > 1 ? 2 : 1); ++dy) {
    const double w = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy); if (w == 0) continue;
    const auto p = plane_hermite((ix + dx) * Y3_.size() + iy + dy, x, y);
    for (int q = 0; q < 9; ++q) s[q] += w * p[q];
  }
  return {s[0], s[1], s[2], {s[3], s[4]}, {s[5], s[6]}, {s[7], s[8]}};
}
AtmosphereState EnvelopeMapAtmosphere::eval(double Teff_b, double g_b, const Composition& c) const {
  const auto v = values(Teff_b, g_b, c);
  AtmosphereState out{};
  out.T = std::exp(v.lnT); out.P = std::exp(v.lnP);
  out.Pgas = out.P - constants::a_rad * std::pow(out.T, 4) / 3;
  out.rho = std::numeric_limits<double>::quiet_NaN();   // not tabulated; the mesh EOS supplies rho at the base
  out.tau = EnvelopeAtmosphere::tau_sentinel;
  out.dlnT_dlnTeff = v.dT[0]; out.dlnT_dlng = v.dT[1]; out.dlnP_dlnTeff = v.dP[0]; out.dlnP_dlng = v.dP[1];
  return out;
}
std::array<double,2> EnvelopeMapAtmosphere::photosphere_R_Teff(double Teff_b, double g_b, const Composition& c) const {
  const auto v = values(Teff_b, g_b, c); const double R = std::exp(v.lnR);
  const double r_b = std::sqrt(constants::G * M_ / g_b);   // L = 4 pi r_b^2 sigma Teff_b^4 = 4 pi R^2 sigma Teff^4
  return {R, Teff_b * std::sqrt(r_b / R)};
}
}  // namespace ember

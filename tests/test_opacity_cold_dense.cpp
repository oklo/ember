#include "ember/opacity_cold_dense.hpp"
#include "ember/opacity_mixture.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace ember;
namespace {
void check(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
void near(double a,double b,const char* message) {
  if(std::abs(a-b)>2e-5*std::max({1.,std::abs(a),std::abs(b)})) {
    std::cerr<<message<<": "<<a<<" != "<<b<<'\n';throw std::runtime_error(message);
  }
}
struct Original : Opacity {
  OpacityState eval(double T,double rho,const Composition& c) const override {
    return {50*std::pow(T/5000.,-.8)*std::pow(rho/.1,.2)*std::exp(.3*c.X[0]+10*c.Z()),-.8,.2,.3,10};
  }
  const char* name() const override {return "analytic original";}
};
Composition comp(double X,double Z) {auto c=solar_scaled(X,Z);c.basis=AbundanceBasis::atomic_mass;return c;}
}
int main()try {
  const auto file=std::filesystem::temp_directory_path()/
      ("ember-cold-opacity-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".dat");
  struct Remove {std::filesystem::path p;~Remove(){std::error_code e;std::filesystem::remove(p,e);}} cleanup{file};
  {
    std::ofstream f(file);f<<std::setprecision(17)<<"2 4 7 synthetic\n-2 -1.5 -1 -.5 0 .5 1\n";
    for(double T:{3000.,5000.,12000.,20000.})f<<std::log10(T)<<' ';
    f<<'\n';
    for(double X:{.98,.9955}) {
      f<<X<<" 0\n";
      for(double T:{3000.,5000.,12000.,20000.}) {
        for(double r:{-2.,-1.5,-1.,-.5,0.,.5,1.})
          f<<std::log10(17.)-.6*std::log10(T/5000)+.4*(r+1)+.7*X/std::log(10.)<<' ';
        f<<'\n';
      }
    }
  }
  Original original;ColdDenseOpacity joined(original,file);
  for(auto [T,R,X,Z]:{std::array{5000.,5.75,.99,0.},std::array{3250.,5.95,.99,0.},
      std::array{11000.,5.8,.99,0.},std::array{18000.,5.8,.99,0.},std::array{6000.,6.1,.9825,0.},
      std::array{6000.,6.1,.99,5e-11},std::array{6000.,6.1,.99,0.}}) {
    const double rho=std::pow(10.,R+3*(std::log10(T)-6)),h=1e-6;
    auto c=comp(X,Z);const auto a=joined.eval(T,rho,c);
    const auto ln=[&](double t,double r,const Composition& q){return std::log(joined.eval(t,r,q).kappa);};
    near(a.dlnk_dlnT,(ln(T*std::exp(h),rho,c)-ln(T*std::exp(-h),rho,c))/(2*h),"T derivative");
    near(a.dlnk_dlnRho,(ln(T,rho*std::exp(h),c)-ln(T,rho*std::exp(-h),c))/(2*h),"rho derivative");
    near(a.dlnk_dX,(ln(T,rho,comp(X+h,Z))-ln(T,rho,comp(X-h,Z)))/(2*h),"H derivative");
    if(Z>0) {
      constexpr double hz=1e-14;
      near(a.dlnk_dZ,(ln(T,rho,comp(X,Z+hz))-ln(T,rho,comp(X,Z-hz)))/(2*hz),"metal join derivative");
    }
  }
  for(auto [T,rho,X,Z]:{std::array{2500.,.004,.99,0.},std::array{15000.,.1,.99,0.},
      std::array{6000.,.01,.99,0.},std::array{6000.,.5,.7,.02}}) {
    auto c=comp(X,Z);const auto a=joined.eval(T,rho,c),b=original.eval(T,rho,c);
    check(a.kappa==b.kappa && a.dlnk_dlnT==b.dlnk_dlnT && a.dlnk_dX==b.dlnk_dX,"inactive source changed");
  }
  auto c=comp(.99,0);
  for(auto [rho,X]:{std::array{11.,.99},std::array{.4,.999}}) {
    bool rejected=false;try{joined.eval(6000,rho,comp(X,0));}catch(const std::domain_error&){rejected=true;}
    check(rejected,"outside dense table accepted");
  }
  // The old 10–12 kK return queried the original source beyond its density
  // edge. The warmer overlap covers this actual cooling-envelope corner.
  struct BoundedOriginal : Original {
    OpacityState eval(double T,double rho,const Composition& c) const override {
      if(std::log10(rho)-3*(std::log10(T)-6)>6)
        throw std::domain_error("original density edge");
      return Original::eval(T,rho,c);
    }
  } bounded;
  ColdDenseOpacity warm(bounded,file);
  for(double T:{9999.,10000.,10010.,12000.,15999.,16000.}) {
    const double rho=std::pow(10.,6.01+3*(std::log10(T)-6));
    const auto v=warm.eval(T,rho,c);
    check(v.kappa>0,"warm dense-gas coverage gap");
  }
  for(double T:{16001.,18000.,19999.,20000.,20001.}) {
    const double rho=std::pow(10.,5.8+3*(std::log10(T)-6)),h=1e-6;
    auto f=[&](double t,double r){return std::log(warm.eval(t,r,c).kappa);};
    const auto v=warm.eval(T,rho,c);
    near(v.dlnk_dlnT,(f(T*std::exp(h),rho)-f(T*std::exp(-h),rho))/(2*h),"warm overlap T derivative");
    near(v.dlnk_dlnRho,(f(T,rho*std::exp(h))-f(T,rho*std::exp(-h)))/(2*h),"warm overlap density derivative");
  }
  // Below 3500 K the original ends at log R = 6 while cool envelopes reach 6.4: the density join must reach
  // the computed source before that edge at every T >= 3000 K (with smooth derivatives through it).
  for(double T:{3001.,3200.,3480.,3499.}) {
    const double rho=std::pow(10.,6.4+3*(std::log10(T)-6));
    check(warm.eval(T,rho,c).kappa>0,"cool dense-gas coverage gap");
  }
  for(double T:{3100.,3300.})for(double R:{5.91,5.94,5.97}) {
    const double rho=std::pow(10.,R+3*(std::log10(T)-6)),h=1e-6;
    auto f=[&](double t,double r){return std::log(warm.eval(t,r,c).kappa);};
    const auto v=warm.eval(T,rho,c);
    near(v.dlnk_dlnT,(f(T*std::exp(h),rho)-f(T*std::exp(-h),rho))/(2*h),"cool overlap T derivative");
    near(v.dlnk_dlnRho,(f(T,rho*std::exp(h))-f(T,rho*std::exp(-h)))/(2*h),"cool overlap density derivative");
  }
  bool corner_rejected=false;   // neither source covers T < 3000 K beyond log R 5.98
  try{warm.eval(2900.,std::pow(10.,6.4+3*(std::log10(2900.)-6)),c);}catch(const std::domain_error&){corner_rejected=true;}
  check(corner_rejected,"unsupported cool corner accepted");
  bool overlap_rejected=false;
  try{warm.eval(18000.,7.,c);}catch(const std::domain_error&){overlap_rejected=true;}
  check(overlap_rejected,"warm overlap skipped unsupported original source");
  ColdDenseOpacity small(original,file,.1),large(original,file,10);
  near(large.eval(6000,.4,c).kappa/small.eval(6000,.4,c).kappa,100,"scale control");
  ElementalOpacity mapped(joined);c.basis=AbundanceBasis::baryon_mass;
  const double source_k=17*std::pow(6000./5000,-.6)*std::pow(.4/.1,.4)*std::exp(.7*c.X[0]);
  check(std::abs(mapped.eval(6000,.4,c).kappa/source_k-1)<1e-11,"source mass convention changed opacity");
  c.X[1]=.003;c.X[2]-=.003;
  const auto a=mapped.eval(6000,.4,c);const double h=1e-6;
  auto plus=c,minus=c;plus.X[1]+=h;plus.X[2]-=h;minus.X[1]-=h;minus.X[2]+=h;
  near(a.dlnk_dY3,(std::log(mapped.eval(6000,.4,plus).kappa)-std::log(mapped.eval(6000,.4,minus).kappa))/(2*h),"isotope derivative");
  std::cout<<"cold dense opacity derivatives, unchanged overlap, source bounds and isotope mapping pass\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}

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
    std::ofstream f(file);f<<std::setprecision(17)<<"2 4 5 synthetic\n-2 -1.5 -1 -.5 0\n";
    for(double T:{3000.,5000.,12000.,20000.})f<<std::log10(T)<<' ';
    f<<'\n';
    for(double X:{.98,.9955}) {
      f<<X<<" 0\n";
      for(double T:{3000.,5000.,12000.,20000.}) {
        for(double r:{-2.,-1.5,-1.,-.5,0.})
          f<<std::log10(17.)-.6*std::log10(T/5000)+.4*(r+1)+.7*X/std::log(10.)<<' ';
        f<<'\n';
      }
    }
  }
  Original original;ColdDenseOpacity joined(original,file);
  for(auto [T,R,X,Z]:{std::array{5000.,5.75,.99,0.},std::array{3250.,5.95,.99,0.},
      std::array{11000.,5.8,.99,0.},std::array{6000.,6.1,.9825,0.},
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
  for(auto [T,rho,X,Z]:{std::array{2500.,.1,.99,0.},std::array{15000.,.1,.99,0.},
      std::array{6000.,.01,.99,0.},std::array{6000.,.5,.7,.02}}) {
    auto c=comp(X,Z);const auto a=joined.eval(T,rho,c),b=original.eval(T,rho,c);
    check(a.kappa==b.kappa && a.dlnk_dlnT==b.dlnk_dlnT && a.dlnk_dX==b.dlnk_dX,"inactive source changed");
  }
  auto c=comp(.99,0);
  for(auto [rho,X]:{std::array{1.1,.99},std::array{.4,.999}}) {
    bool rejected=false;try{joined.eval(6000,rho,comp(X,0));}catch(const std::domain_error&){rejected=true;}
    check(rejected,"outside dense table accepted");
  }
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

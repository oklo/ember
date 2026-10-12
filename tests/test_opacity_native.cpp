#include "ember/opacity_radiative.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace ember;

int main() {
  const auto directory=std::filesystem::temp_directory_path()/
      ("ember-native-opacity-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(directory);
  try {
    const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
    const std::vector<double> metals{0.,.02,.03,.06,.10,.16,.24};
    auto family=[&](const std::string& name,double offset,double maximum_share=1.) {
      const auto path=directory/(name+".dat");std::ofstream manifest(path);
      manifest<<std::setprecision(17)<<"EMBER_OPACITY_MIXTURE 1 "<<metals.size()<<" logRho analytic\n";
      for(std::size_t iz=0;iz<metals.size();++iz) {
        const double z=metals[iz];const auto file=name+"_"+std::to_string(iz)+".dat";
        manifest<<z<<' '<<std::quoted(file)<<'\n';
        std::ofstream table(directory/file);table<<std::setprecision(17);
        table<<"4 9 9 native composition control\n-5 -4 -3 -2 -1 0 1 2 3\n3 3.5 4 4.5 5 5.5 6 6.5 7\n";
        for(double fraction:{0.,.4,.8,1.}) {
          const double u=fraction*maximum_share;
          table<<u*(1-z)<<' '<<z<<'\n';
          for(int it=0;it<9;++it) {
            for(int ir=0;ir<9;++ir)
              table<<offset+.4*u+40*z*z+.3*(-5.+ir)-.2*(3.+.5*it)<<' ';
            table<<'\n';
          }
        }
      }
      return path;
    };
    const RadiativeOpacity::Tables files{family("low",0),family("warm",.2),
        family("bridge",.4),family("hot",.6)};
    RadiativeOpacity opacity(files,RadiativeOpacity::NativeComposition{});
    const auto c=[](double x,double z){return solar_scaled(x,z);};
    auto expected=[&](double T,double rho,double x,double z,double offset) {
      const auto i=interp::locate(metals,z);
      const double a=metals[i],b=metals[i+1],f=(z-a)/(b-a);
      return offset+.4*x/(1-z)+40*((1-f)*a*a+f*b*b)+.3*std::log10(rho)-.2*std::log10(T);
    };
    // The strong Z curvature makes ignoring the .06/.10/.16/.24 source planes detectable.
    for(const auto [lt,offset]:std::array<std::pair<double,double>,4>{{{4.,0},{4.7,.2},{5.3,.4},{6.,.6}}})
      for(double z:{0.,.06,.08241129899325579,.10,.22,.24}) {
        const double x=.6*(1-z),T=std::pow(10.,lt),rho=.1;
        const auto state=opacity.eval(T,rho,c(x,z));
        require(std::abs(std::log10(state.kappa)-expected(T,rho,x,z,offset))<1e-11,
            "native reader must use the supplied high-Z planes");
      }
    constexpr double h=1e-6,x=.6475887010067442,z=.08241129899325579;
    HydrogenShareMixtureOpacity source(files.low);
    const auto source_state=source.eval(1e4,.1,c(x,z));
    require(std::abs((std::log(source.eval(1e4,.1,c(x,z+h)).kappa)-
        std::log(source.eval(1e4,.1,c(x,z-h)).kappa))/(2*h)-source_state.dlnk_dZ)<2e-7,
        "source metal derivative at fixed hydrogen mismatch");
    for(double lt:{4.,4.43,4.7,5.075,5.3,5.65,6.}) {
      const double T=std::pow(10.,lt),rho=.1;
      const auto s=opacity.eval(T,rho,c(x,z));
      const auto ln=[&](double t,double r,double xx,double zz){return std::log(opacity.eval(t,r,c(xx,zz)).kappa);};
      require(std::abs((ln(T*std::exp(h),rho,x,z)-ln(T*std::exp(-h),rho,x,z))/(2*h)-s.dlnk_dlnT)<2e-7,
          "temperature derivative must include table joins");
      require(std::abs((ln(T,rho*std::exp(h),x,z)-ln(T,rho*std::exp(-h),x,z))/(2*h)-s.dlnk_dlnRho)<2e-7,
          "density derivative mismatch");
      require(std::abs((ln(T,rho,x+h,z)-ln(T,rho,x-h,z))/(2*h)-s.dlnk_dX)<2e-7,
          "hydrogen-share chain rule mismatch");
      // ElementalOpacity intentionally does not expose an arbitrary input-metal
      // derivative. Check its finite response against the independent source law.
      const double dz=std::log(10.)*(.4*x/((1-z)*(1-z))+40*(.06+.10));
      require(std::abs((ln(T,rho,x,z+h)-ln(T,rho,x,z-h))/(2*h)-dz)<2e-7,
          "assembled metal response mismatch");
      const auto range=opacity.density_range(T,c(x,z));
      require(range && range->min<=rho && range->max>=rho,"source density range missing");
    }
    int rejected=0;
    for(const auto point:std::array<std::array<double,4>,3>{{{1e4,.1,.5,.241},{1e4,1e4,.5,.08},{100.,.1,.5,.08}}})
      try{opacity.eval(point[0],point[1],c(point[2],point[3]));}catch(const std::domain_error&){++rejected;}
    require(rejected==3,"native mode must not extend Z, temperature or density");
    const auto limited=family("limited",0,.7);
    RadiativeOpacity limited_opacity({limited,limited,limited,limited},RadiativeOpacity::NativeComposition{});
    bool caught=false;
    try{limited_opacity.eval(1e4,.1,c(.8*(1-.08),.08));}catch(const std::domain_error&){caught=true;}
    require(caught,"native mode must not silently extend the hydrogen axis");
    // Selecting native composition must retain the independent cold-gas source.
    const auto cold_file=directory/"cold.dat";
    {
      std::ofstream cold(cold_file);cold<<std::setprecision(17);
      cold<<"2 4 7 native cold control\n-2 -1.5 -1 -.5 0 .5 1\n";
      for(double t:{3000.,5000.,12000.,20000.})cold<<std::log10(t)<<' ';
      cold<<'\n';
      for(double hydrogen:{.98,.9955}) {
        cold<<hydrogen<<" 0\n";
        for(double t:{3000.,5000.,12000.,20000.}) {
          for(double r:{-2.,-1.5,-1.,-.5,0.,.5,1.})
            cold<<std::log10(17.)-.6*std::log10(t/5000)+.4*(r+1)+.7*hydrogen/std::log(10.)<<' ';
          cold<<'\n';
        }
      }
    }
    auto cold_files=files;cold_files.cold_dense=cold_file;
    RadiativeOpacity cooling(cold_files,RadiativeOpacity::NativeComposition{});
    const double mass_scale=.99/nuclides[0].A+.01*4/nuclides[2].A;
    const double source_hydrogen=(.99/nuclides[0].A)/mass_scale;
    const double cold_expected=17*std::pow(4*mass_scale,.4)*std::exp(.7*source_hydrogen)*mass_scale;
    require(std::abs(cooling.eval(5000.,.4,c(.99,0)).kappa/cold_expected-1)<1e-10,
        "native composition lost the cold dense-gas source");
    require(cooling.eval(1e6,.1,c(x,z)).kappa==opacity.eval(1e6,.1,c(x,z)).kappa,
        "inactive cold source changed native high-Z opacity");
    std::filesystem::remove_all(directory);
    std::cout<<"PASS native opacity: high-Z source planes, joins, composition derivatives and coverage limits\n";
    return 0;
  }catch(const std::exception& e) {
    std::filesystem::remove_all(directory);std::cerr<<e.what()<<'\n';return 1;
  }
}

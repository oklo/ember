#include "ember/fingering.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace ember;

int main() {
  try {
    auto require=[](bool ok,const char* message) {
      if(!ok)throw std::runtime_error(message);
    };
    const auto stable=brown_fingering_flux(1e-6,1e-7,1e7);
    require(stable.regime==FingeringRegime::stable && stable.growth_rate==0
      && stable.chemical_nusselt_excess==0,"stable threshold must not mix");
    const auto stable_above=brown_fingering_flux(1e-6,1e-7,2e7);
    require(stable_above.regime==FingeringRegime::stable,"stable side of threshold");
    for(double R:{.5,1.}) {
      const auto s=brown_fingering_flux(1e-6,1e-7,R);
      require(s.regime==FingeringRegime::overturning && s.chemical_nusselt_excess==0,
        "ordinary convection must not be counted as fingering");
    }
    for(auto bad:std::array<std::array<double,3>,7>{{
      {0,1e-7,2},{1e-6,0,2},{1e-6,1,2},{1e-6,1e-7,0},
      {-1,1e-7,2},{1e-6,1e-7,std::numeric_limits<double>::infinity()},
      {std::numeric_limits<double>::quiet_NaN(),1e-7,2}}}) {
      bool caught=false;
      try {brown_fingering_flux(bad[0],bad[1],bad[2]);}
      catch(const std::invalid_argument&){caught=true;}
      require(caught,"invalid physical parameters must be rejected");
    }
    const auto poorer_fit=brown_fingering_flux(1e-6,2e-6,2.);
    require(!poorer_fit.within_calibrated_ratio,"tau>Pr must flag calibration limitation");

    std::ifstream data(std::string(EMBER_TEST_DATA_DIR)+"/fingering_reference.txt");
    require(bool(data),"missing independent eigenvalue reference");
    std::string line;int count=0;double maximum_error=0,maximum_derivative_error=0;
    const auto quantities=[](const FingeringFlux& f) {
      return std::array<double,4>{f.growth_rate,f.wavenumber_squared,
        f.thermal_nusselt_excess,f.chemical_nusselt_excess};
    };
    while(std::getline(data,line)) {
      if(line.empty() || line[0]=='#')continue;
      std::istringstream in(line);double pr,tau,R,lambda,q,thermal,chemical;
      require(bool(in>>pr>>tau>>R>>lambda>>q>>thermal>>chemical),"bad reference row");
      const auto s=brown_fingering_flux(pr,tau,R);
      require(s.regime==FingeringRegime::fingering,"unstable reference must mix");
      const std::array<double,4> actual{s.growth_rate,s.wavenumber_squared,
        s.thermal_nusselt_excess,s.chemical_nusselt_excess};
      const std::array<double,4> expected{lambda,q,thermal,chemical};
      for(std::size_t i=0;i<actual.size();++i) {
        const double error=std::abs(actual[i]/expected[i]-1);
        maximum_error=std::max(maximum_error,error);
        if(error>2e-5) {
          std::cerr<<"reference discrepancy: "<<pr<<' '<<tau<<' '<<R<<" field="<<i
            <<" actual="<<actual[i]<<" expected="<<expected[i]<<'\n';
          throw std::runtime_error("independent eigenvalue comparison failed");
        }
      }
      const auto response=brown_fingering_response(pr,tau,R);
      require(response.derivatives_defined && quantities(response.value)==actual,
        "response preserves the independently checked flux values");
      // Five-point differences differentiate independently evaluated fastest
      // modes. Perturbations remain strictly in the open unstable regime.
      const double step=std::min({2e-3,std::log(R)/32,std::log(1/(tau*R))/32});
      require(step>0,"derivative comparison must stay in its physical regime");
      for(std::size_t j=0;j<3;++j) {
        std::array<std::array<double,4>,4> nearby{};
        const std::array<double,4> offset{-2,-1,1,2};
        for(std::size_t k=0;k<4;++k) {
          std::array<double,3> input{pr,tau,R};input[j]*=std::exp(offset[k]*step);
          const auto evaluated=brown_fingering_flux(input[0],input[1],input[2]);
          require(evaluated.regime==FingeringRegime::fingering,"finite difference crossed stability boundary");
          nearby[k]=quantities(evaluated);
        }
        for(std::size_t i=0;i<4;++i) {
          const double numerical=(nearby[0][i]-8*nearby[1][i]+8*nearby[2][i]-nearby[3][i])/(12*step);
          const double scale=std::max(actual[i],std::abs(response.partials[i][j]));
          const double error=std::abs(numerical-response.partials[i][j])/scale;
          maximum_derivative_error=std::max(maximum_derivative_error,error);
          if(error>3e-6) {
            std::cerr<<"response discrepancy: "<<pr<<' '<<tau<<' '<<R<<" field="<<i
              <<" parameter="<<j<<" analytic="<<response.partials[i][j]<<" numerical="<<numerical<<'\n';
            throw std::runtime_error("independent derivative comparison failed");
          }
        }
      }
      ++count;
    }
    require(count>=90,"incomplete reference set");
    const auto near=brown_fingering_flux(1e-6,1e-7,.9999e7);
    const auto nearer=brown_fingering_flux(1e-6,1e-7,.99999e7);
    require(nearer.chemical_nusselt_excess<.011*near.chemical_nusselt_excess,
      "mixing must vanish quadratically at the stable threshold");
    require(!brown_fingering_response(1e-6,1e-7,1).derivatives_defined &&
      !brown_fingering_response(1e-6,1e-7,1e7).derivatives_defined,
      "do not claim open-regime derivatives at stability boundaries");
    for(double pr:{1e-6,.015})for(double tau:{1e-7,1e-4})for(double net:{.001,.1,.9}) {
      const auto a=two_composition_fingering(pr,{tau,tau},{1.,net-1.});
      const auto b=brown_fingering_flux(pr,tau,1/(1+(net-1.)));
      require(std::abs(a.growth_rate/b.growth_rate-1)<1e-12 &&
          std::abs(a.thermal_nusselt_excess/b.thermal_nusselt_excess-1)<1e-12 &&
          std::abs(a.mixing_over_thermal[0]/(tau*b.chemical_nusselt_excess)-1)<1e-12 &&
          a.mixing_over_thermal[0]==a.mixing_over_thermal[1],
          "equal composition diffusivities must recover Brown heat and mixing");
    }
    require(two_composition_fingering(.01,{1e-4,2e-4},{1e-4,0}).regime==FingeringRegime::stable,
        "two-component stationary threshold must not mix");
    require(two_composition_fingering(.01,{1e-4,2e-4},{-.1,-.1}).regime==FingeringRegime::stable,
        "two stabilizing fields must not mix");
    require(two_composition_fingering(.01,{1e-4,2e-4},{-.1,.01}).regime==FingeringRegime::stable,
        "Routh-Hurwitz test must recognize stable fast-field driving");
    for(std::size_t field=0;field<2;++field) {
      const std::array<double,2> tau{1e-5,1e-4};std::array<double,2> drive{};drive[field]=.03;
      const auto a=two_composition_fingering(.015,tau,drive);
      const auto b=brown_fingering_flux(.015,tau[field],1/.03);
      require(std::abs(a.growth_rate/b.growth_rate-1)<1e-10 &&
          std::abs(a.wavenumber_squared/b.wavenumber_squared-1)<1e-9,
          "a passive composition field must leave the active Brown mode unchanged");
    }
    const auto cancel=two_composition_fingering(.015,{1e-4,1.04e-4},{.05,-.05});
    require(cancel.regime==FingeringRegime::fingering && cancel.growth_rate>0 &&
        cancel.mixing_over_thermal[0]>cancel.mixing_over_thermal[1],
        "opposing diffusion modes can grow at zero net composition buoyancy");
    for(const auto& input:std::array<std::array<double,5>,2>{{
        {.0004,.00005,.0034,-1.1,.83},
        {.013569673496689305,.00011105469291715423,.00011544183077567232,
          -.18989186929200672,.19034483027947444}}}) {
      bool oscillatory=false;
      try{two_composition_fingering(input[0],{input[1],input[2]},{input[3],input[4]});}
      catch(const std::domain_error&){oscillatory=true;}
      require(oscillatory,"an actually dominant oscillatory mode must be reported");
    }
    std::ifstream oscillatory_data(std::string(EMBER_TEST_DATA_DIR)+"/fingering_oscillatory_reference.txt");
    require(bool(oscillatory_data),"missing independent oscillatory reference");
    int oscillatory_count=0;
    while(std::getline(oscillatory_data,line)) {
      if(line.empty() || line[0]=='#')continue;
      double pr,t0,t1,g0,g1,lambda,q,omega;std::istringstream in(line);
      require(bool(in>>pr>>t0>>t1>>g0>>g1>>lambda>>q>>omega),"bad oscillatory reference row");
      const auto a=two_composition_fingering(pr,{t0,t1},{g0,g1},OscillatoryMixing::growth_squared);
      const auto b=two_composition_fingering(pr,{t0,t1},{g0,g1},OscillatoryMixing::growth_frequency);
      require(a.regime==FingeringRegime::fingering &&
          std::abs(a.growth_rate/lambda-1)<3e-6 && std::abs(a.wavenumber_squared/q-1)<2e-4 &&
          std::abs(a.oscillation_frequency-omega)/std::max(lambda,omega)<3e-6,
          "oscillatory mode disagrees with independent eigenvalues");
      require(a.growth_rate==b.growth_rate && a.oscillation_frequency==b.oscillation_frequency &&
          a.wavenumber_squared==b.wavenumber_squared,"saturation choice changed the linear mode");
      const double ratio=std::hypot(a.growth_rate,a.oscillation_frequency)/a.growth_rate;
      require(a.velocity_squared>0 && std::abs(b.velocity_squared/a.velocity_squared/ratio-1)<1e-12 &&
          std::abs(b.thermal_nusselt_excess/a.thermal_nusselt_excess/ratio-1)<1e-12,
          "continuous oscillatory saturation has inconsistent heat amplitude");
      for(std::size_t i=0;i<2;++i)
        require(a.mixing_over_thermal[i]>0 &&
            std::abs(b.mixing_over_thermal[i]/a.mixing_over_thermal[i]/ratio-1)<1e-12,
            "oscillatory composition transport has inconsistent phase response");
      ++oscillatory_count;
    }
    require(oscillatory_count>=80,"incomplete oscillatory reference set");
    for(auto policy:{OscillatoryMixing::reject,OscillatoryMixing::growth_squared,OscillatoryMixing::growth_frequency}) {
      const auto f=two_composition_fingering(.015,{1e-4,1e-4},{.04,-.01},policy);
      require(f.oscillation_frequency==0 && f.velocity_squared>0,"stationary limit acquired an oscillation");
      require(two_composition_fingering(.015,{1e-4,2e-4},{-.1,-.1},policy).velocity_squared==0,
          "stable state acquired a saturation amplitude");
    }
    std::cout<<"PASS: "<<oscillatory_count<<" independent oscillatory cases and both explicit saturation choices\n";
    std::ifstream two_data(std::string(EMBER_TEST_DATA_DIR)+"/fingering_two_reference.txt");
    require(bool(two_data),"missing two-composition eigenvalue reference");
    int two_count=0;double two_error=0;
    while(std::getline(two_data,line)) {
      if(line.empty() || line[0]=='#')continue;
      double pr,t0,t1,g0,g1,lambda,q;std::istringstream in(line);
      require(bool(in>>pr>>t0>>t1>>g0>>g1>>lambda>>q),"bad two-composition reference row");
      const auto f=two_composition_fingering(pr,{t0,t1},{g0,g1});
      const double error=std::abs(f.growth_rate/lambda-1);two_error=std::max(two_error,error);
      require(f.regime==FingeringRegime::fingering && error<2e-6 &&
          std::abs(f.wavenumber_squared/q-1)<2e-4,"two-field mode disagrees with independent eigenvalues");
      require(f.thermal_nusselt_excess>0 && f.mixing_over_thermal[0]>0 && f.mixing_over_thermal[1]>0,
          "two-field saturation must preserve transport signs");
      ++two_count;
    }
    require(two_count>=25,"incomplete two-composition reference set");
    std::cout<<"PASS: "<<two_count<<" two-composition eigenvalue cases; maximum growth error "<<two_error<<'\n';
    std::cout<<"PASS: "<<count<<" independent spectral cases; maximum relative error "
      <<maximum_error<<"; maximum scaled derivative error "<<maximum_derivative_error
      <<"; instability boundaries, invalid inputs and calibration flag pass\n";
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<e.what()<<'\n';return 1;
  }
}

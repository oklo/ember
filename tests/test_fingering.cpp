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
    std::cout<<"PASS: "<<count<<" independent spectral cases; maximum relative error "
      <<maximum_error<<"; maximum scaled derivative error "<<maximum_derivative_error
      <<"; instability boundaries, invalid inputs and calibration flag pass\n";
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<e.what()<<'\n';return 1;
  }
}

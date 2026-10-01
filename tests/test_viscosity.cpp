#include "ember/viscosity.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

int main() try {
  auto require=[](bool ok,const char* message) {if(!ok)throw std::runtime_error(message);};
  std::ifstream input(std::string(EMBER_TEST_DATA_DIR)+"/liquid_viscosity_reference.txt");
  require(bool(input),"missing independent viscosity reference");
  std::string line;int count=0;double worst=0,derivative_error=0;
  while(std::getline(input,line)) {
    if(line.empty()||line[0]=='#')continue;
    std::array<double,4> x{};double e,i,log;
    std::istringstream row(line);require(bool(row>>x[0]>>x[1]>>x[2]>>x[3]>>e>>i>>log),"bad viscosity reference");
    const auto r=ember::ocp_liquid_viscosity(x[0],x[1],x[2],x[3]);
    for(double error:{std::abs(r.value.electron/e-1),std::abs(r.value.ion_classical/i-1),
                      std::abs(r.value.coulomb_log/log-1)}) {
      worst=std::max(worst,error);require(error<2e-6,"independent viscosity quadrature mismatch");
    }
    for(int j=0;j<(x[2]==x[3]?2:4);++j) {
      const double h=1e-4;auto lo=x,hi=x;lo[j]*=std::exp(-h);hi[j]*=std::exp(h);
      const auto a=ember::ocp_liquid_viscosity(lo[0],lo[1],lo[2],lo[3]).value;
      const auto b=ember::ocp_liquid_viscosity(hi[0],hi[1],hi[2],hi[3]).value;
      const std::array<double,2> numerical{(b.electron-a.electron)/(2*h),(b.ion_classical-a.ion_classical)/(2*h)};
      for(int k=0;k<2;++k) {
        const double error=std::abs(numerical[k]-r.partials[k][j])/(k?i:e);
        derivative_error=std::max(derivative_error,error);require(error<2e-6,"viscosity derivative mismatch");
      }
    }
    ++count;
  }
  require(count==8,"incomplete viscosity reference");
  for(auto x:std::array<std::array<double,4>,4>{{{0,1000,1,1},{1e6,1,1,1},{1e4,1e6,6,12},{1e6,1e8,1,1}}}) {
    bool caught=false;try {ember::ocp_liquid_viscosity(x[0],x[1],x[2],x[3]);}
    catch(const std::exception&){caught=true;}
    require(caught,"invalid or unsupported viscosity state accepted");
  }
  const auto h=ember::ocp_liquid_viscosity(159000,911.7,1,1).value;
  require(h.electron>40*h.ion_classical && h.electron_electron_frequency<.002*h.electron_ion_frequency,
          "H-layer viscosity must be dominated by electron-ion scattering");
  std::cout<<"quadrature error "<<worst<<", derivative error "<<derivative_error<<'\n';
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}

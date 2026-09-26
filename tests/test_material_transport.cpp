#include "ember/material_transport.hpp"
#include "material_test_eos.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace ember;
namespace {
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
template<class F> void rejects(F call){bool failed=false;try{call();}catch(const std::exception&){failed=true;}require(failed,"invalid input accepted");}
}
int main() {
  try {
    // Independent scalar thermal energy balance, solved by bisection.
    AnalyticMaterial eos{1,{4,0,0},{2,0,0},{},{1,1},0,1};
    std::vector<MaterialVector> old{{std::log(2.),0,0},{0,0,0}};
    std::vector<double> mass{.3,.7};MaterialMatrix heat{};heat[0][0]=.2;
    MaterialTransportOptions options;options.components=1;
    auto faces=[&](auto){return std::vector<MaterialMatrix>{heat};};
    const auto answer=transport_material(old,old,mass,.4,eos,faces,options);
    double lo=1.3,hi=2;
    for(int i=0;i<80;++i) {
      const double T0=.5*(lo+hi),T1=(1.3-.3*T0)/.7;
      const double residual=.3*.375*(T0-2)+.4*.2*(1/T1-1/T0);
      if(residual>0)hi=T0;else lo=T0;
    }
    require(std::abs(std::exp(answer.primitives[0][0])-.5*(lo+hi))<1e-10,"scalar heat root mismatch");
    require(std::abs(answer.integrated_conservation_error[0])<1e-11,"heat not conserved");
    require(answer.entropy_change+1e-11>=answer.backward_euler_entropy_bound,"entropy inequality failed");
    // Stiff uniform modes, non-diagonal unequal capacities. A uniform vector
    // has exactly zero face differences even at very large conductance.
    MaterialVector uniform{.3,-.8,1.1};
    for(double stiffness:{0.,1.,1e6,1e12}) {
      std::vector<MaterialMatrix> a(16),k(15);std::vector<MaterialVector> rhs(16);
      for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<3;++j)for(std::size_t v=0;v<3;++v) {
        a[i][j][v]=j==v?1.0+static_cast<double>(i+j):.1;
        rhs[i][j]+=a[i][j][v]*uniform[v];
      }
      for(auto& face:k)for(std::size_t j=0;j<3;++j)face[j][j]=stiffness;
      const auto result=solve_material_chain(a,k,rhs,3);
      for(const auto& cell:result)for(std::size_t j=0;j<3;++j)
        require(std::abs(cell[j]-uniform[j])<2e-10,"stiff uniform mode lost");
    }
    // A disconnected chain preserves each cell separately, including one cell.
    heat={};const auto sealed=transport_material(old,old,mass,1e12,eos,faces,options);
    require(sealed.primitives==old,"disconnected cells changed");
    rejects([&]{auto bad=mass;bad[0]=0;transport_material(old,old,bad,1,eos,faces,options);});
    heat[0][0]=-1;rejects([&]{transport_material(old,old,mass,1,eos,faces,options);});
    rejects([&]{solve_material_chain({}, {}, {},3);});
    std::cout<<"material transport controls passed\n";return 0;
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

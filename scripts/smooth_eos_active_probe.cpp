#include "ember/eos_smooth_mixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Input: X Y3 T rho selected_directions (bit 0 H1, bit 1 He3).
// Unselected derivatives are serialized as null; selected ones must be finite.
int main(int argc,char** argv) {
  if(argc!=2)return 2;
  try {
    ember::SmoothMetalHelmholtzEos eos(argv[1],
        ember::HelmholtzTableEos::Mixture::allow_documented_proxy);
    std::cout<<std::setprecision(17);
    double x,y,T,rho;int mask;
    while(std::cin>>x>>y>>T>>rho>>mask) {
      try {
        if(mask<0 || mask>3)throw std::invalid_argument("invalid active direction mask");
        auto c=ember::solar_scaled(x,.02);c.X[1]=y;c.X[2]-=y;
        c.basis=ember::AbundanceBasis::baryon_mass;c.metal_inventory=ember::MetalInventory::gs98;
        const std::array<bool,2> active{bool(mask&1),bool(mask&2)};
        const auto f=eos.active_composition_potential(T,rho,c,active);
        const auto range=eos.density_range(T,c).value();
        if(!std::isfinite(f.phi))throw std::runtime_error("nonfinite potential");
        auto value=[&](double v,bool selected) {
          if(selected) {
            if(!std::isfinite(v))throw std::runtime_error("nonfinite active derivative");
          } else if(!std::isnan(v))throw std::runtime_error("unselected derivative is not NaN");
        };
        for(unsigned k=0;k<2;++k) {
          value(f.gradient[k],active[k]);value(f.dgradient_dlnT[k],active[k]);
          value(f.dgradient_dlnRho[k],active[k]);
          for(unsigned l=0;l<2;++l)value(f.hessian[k][l],active[k]&&active[l]);
        }
        const double v[]={f.phi,f.gradient[0],f.gradient[1],f.hessian[0][0],
            f.hessian[0][1],f.hessian[1][0],f.hessian[1][1],
            f.dgradient_dlnT[0],f.dgradient_dlnT[1],
            f.dgradient_dlnRho[0],f.dgradient_dlnRho[1],range.min,range.max};
        std::cout<<"{\"ok\":true,\"values\":[";
        for(unsigned i=0;i<13;++i) {
          if(i)std::cout<<',';
          if(std::isnan(v[i]))std::cout<<"null";else std::cout<<v[i];
        }
        std::cout<<"]}\n";
      } catch(const std::exception& e) {
        std::cout<<"{\"ok\":false,\"error\":"<<std::quoted(e.what())<<"}\n";
      }
    }
    return std::cin.eof()?0:2;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

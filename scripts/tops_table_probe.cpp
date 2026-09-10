// Probe the actual runtime interpolation in a fixed-Z TOPS source family.
// Input: atomic X, T [K], rho [g/cm3]. No isotope mapping or AESOPUS blend.
#include "ember/opacity_table.hpp"
#include "ember/opacity_blend.hpp"
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {
  try {
    if(argc!=3)throw std::invalid_argument("tops_table_probe low_table high_table");
    ember::TabulatedOpacity low(argv[1],"TOPS low",ember::TabulatedOpacity::DensityAxis::logRho);
    ember::TabulatedOpacity high(argv[2],"TOPS high",ember::TabulatedOpacity::DensityAxis::logRho);
    if(low.metallicity()!=high.metallicity())throw std::invalid_argument("source metallicities differ");
    ember::BlendedOpacity source(low,high,5.6,5.7);
    double x,T,rho;
    std::cout<<std::setprecision(17);
    while(std::cin>>x>>T>>rho) {
      const auto c=ember::solar_scaled(x,low.metallicity());
      const auto state=source.eval(T,rho,c);
      std::cout<<state.kappa<<' '<<state.dlnk_dlnT<<' '<<state.dlnk_dlnRho<<'\n';
    }
    if(!std::cin.eof())throw std::invalid_argument("malformed probe input");
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

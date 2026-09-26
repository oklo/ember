// Compare two masked atmosphere tables while loading the selected EOS once.
#include "ember/atmosphere_grid.hpp"
#include "ember/eos_smooth_mixture.hpp"
#include <cmath>
#include <cstdio>

using namespace ember;

void emit(const CompositionAtmosphereGrid& grid, double t, double g, const Composition& c) {
  if (!grid.covers(t,g,c)) { std::printf("{\"covered\":false}"); return; }
  try {
    const auto s=grid.eval(t,g,c); const auto r=grid.composition_response(t,g,c);
    std::printf("{\"covered\":true,\"T\":%.17g,\"Pgas\":%.17g,\"P\":%.17g,\"rho\":%.17g,"
                "\"logarithmic_thermal_derivatives\":[%.17g,%.17g,%.17g,%.17g],"
                "\"composition_derivatives\":[%.17g,%.17g,%.17g,%.17g]}",
                s.T,s.Pgas,s.P,s.rho,s.dlnT_dlnTeff,s.dlnT_dlng,s.dlnP_dlnTeff,s.dlnP_dlng,
                r.dlnT_dXH,r.dlnT_dX3,r.dlnP_dXH,r.dlnP_dX3);
  } catch (const std::domain_error& e) {
    std::fprintf(stderr,"atmosphere/EOS domain: %s\n",e.what());
    std::printf("{\"covered\":false,\"atmosphere_covered\":true,\"eos_supported\":false}");
  }
}

int main(int argc,char** argv) {
  try {
    if(argc!=4) throw std::invalid_argument("EOS family and two atmosphere tables required");
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    CompositionAtmosphereGrid before(eos,argv[2],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    CompositionAtmosphereGrid after(eos,argv[3],CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    double x,y,t,lg;int fields;
    while((fields=std::scanf("%lf %lf %lf %lf",&x,&y,&t,&lg))==4) {
      if(!std::isfinite(x+y+t+lg)||t<=0)throw std::invalid_argument("invalid atmosphere query");
      auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
      c.X[1]=y;c.X[2]-=y;
      std::printf("{\"before\":");emit(before,t,std::pow(10.,lg),c);
      std::printf(",\"after\":");emit(after,t,std::pow(10.,lg),c);std::puts("}");
    }
    if(fields!=EOF)throw std::invalid_argument("expected four source coordinates");
  } catch(const std::exception& e) {std::fprintf(stderr,"%s\n",e.what());return 1;}
}

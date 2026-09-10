// Compare source-atmosphere density with Ember's EOS at the same T/Pgas.
// Input rows: XH, XHe3, T [K], Pgas [dyn/cm2], source rho [g/cm3].
// The fixed Z=.02 inventory is the one used by the 0.1 Msun grid.
#include "ember/constants.hpp"
#include "ember/eos_composition.hpp"
#include "ember/eos_mixture.hpp"
#include <cmath>
#include <cstdio>

int main(int argc,char** argv) {
  using namespace ember;
  try {
    const bool metal=argc==2 && std::string(argv[1])=="--gs98";
    if(argc!=1 && !metal)throw std::invalid_argument("optional --gs98 only");
    std::unique_ptr<Eos> selected;
    if(metal)selected=std::make_unique<MetalHelmholtzEos>("data/eos/freeeos300_gs98_z020.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
    else selected=std::make_unique<CompositionHelmholtzEos>("data/eos/freeeos300_hhe_extended.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
    const auto& eos=*selected;
    double x, y, t, pg, source;
    int fields;
    while ((fields = std::scanf("%lf %lf %lf %lf %lf", &x, &y, &t, &pg,
                               &source)) == 5) {
      auto composition = solar_scaled(x, .02);
      composition.basis = AbundanceBasis::baryon_mass;
      if(metal)composition.metal_inventory=MetalInventory::gs98;
      composition.X[1] = y;
      composition.X[2] -= y;
      const double rho = eos.rho_from_PT(
          t, pg + constants::a_rad * std::pow(t, 4) / 3, composition);
      std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", x, y, t, pg,
                  source, rho, rho / source - 1);
    }
    if (fields != EOF)
      throw std::invalid_argument("expected five numbers per source state");
  } catch (const std::exception &error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }
}

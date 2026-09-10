// Read full-profile rows from stdin. Emit local physical diagnostics as JSON arrays.
#include "ember/conduction_table.hpp"
#include "ember/constants.hpp"
#include "ember/convection.hpp"
#include "ember/eos_composition.hpp"
#include "ember/eos_mixture.hpp"
#include "ember/opacity_mixture.hpp"
#include <cmath>
#include <cstdio>
int main(int argc,char** argv) {
  using namespace ember;
  try {
    bool metal=false;
    std::string eos_family="data/eos/freeeos300_gs98_z020.dat";
    std::string opacity_directory="data/opacity";
    bool explicit_family=false;
    for(int i=1;i<argc;++i) {
      const std::string arg=argv[i];
      if(arg=="--gs98")metal=true;
      else if(arg=="--eos-family" && i+1<argc) {
        eos_family=argv[++i];explicit_family=true;
      } else if(arg=="--opacity-directory" && i+1<argc)opacity_directory=argv[++i];
      else throw std::invalid_argument("expected --gs98, --eos-family PATH or --opacity-directory PATH");
    }
    if(explicit_family && !metal)throw std::invalid_argument("explicit metal EOS family requires --gs98");
    std::unique_ptr<Eos> selected;
    if(metal)selected=std::make_unique<MetalHelmholtzEos>(eos_family,HelmholtzTableEos::Mixture::allow_documented_proxy);
    else selected=std::make_unique<CompositionHelmholtzEos>("data/eos/freeeos300_hhe_extended.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
    const auto& eos=*selected;
    StellarMixtureOpacity rad(opacity_directory);
    const std::string suffix=metal?"_metals.dat":".dat";
    TabulatedConduction cond("data/conduction/condtab21wd"+suffix), classic("data/conduction/condtab21_I"+suffix),
        undamped("data/conduction/condtab21nd"+suffix);
    HotConduction hot(cond);
    double m, r, rho, T, L, X, Y3, Y4;
    while (std::scanf("%lf %lf %lf %lf %lf %lf %lf %lf", &m, &r, &rho, &T, &L, &X, &Y3, &Y4) == 8) {
      auto c = solar_scaled(X, .02);
      c.basis = AbundanceBasis::baryon_mass;
      if(metal)c.metal_inventory=MetalInventory::gs98;
      c.X[1] = Y3;
      c.X[2] = Y4;
      const auto e = eos.eval(T, rho, c);
      const auto kr = rad.eval(T, rho, c);
      const auto kc = hot.eval(T, rho, c);
      const double k = 1 / (1 / kr.kappa + 1 / kc.kappa);
      const double gr = 3 * k * L * e.P /
                        (16 * M_PI * constants::a_rad * constants::c * constants::G * m * std::pow(T, 4)),
                   gravity = constants::G * m / (r * r);
      const auto conv = mixing_length_gradient(gr, e.grad_ad, mixing_length_U(T, rho, k, gravity, e, 1.9));
      const auto bounds = rad.density_range(T, c).value();
      std::printf("[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]\n",
                  e.P, e.E, e.S, e.grad_ad, kr.kappa, std::isfinite(kc.kappa) ? kc.kappa : 0., k, gr,
                  conv.grad, gr > 0 ? k / kc.kappa * conv.grad / gr : 0., bounds.min, bounds.max,
                  T >= 1e6 ? classic.eval(T, rho, c).kappa : 0.,
                  T >= 1e6 ? undamped.eval(T, rho, c).kappa : 0.);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}

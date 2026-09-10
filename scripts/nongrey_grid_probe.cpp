// Query the actual masked atmosphere interpolator and selected GS98 EOS.
// Input: XH X3 Teff[K] log10(g[cm/s2]). Output is one JSON object per input.
#include "ember/atmosphere_grid.hpp"
#include "ember/eos_mixture.hpp"
#include <cmath>
#include <cstdio>

int main(int argc, char **argv) {
  using namespace ember;
  try {
    if (argc != 3)
      throw std::invalid_argument("EOS manifest and atmosphere table required");
    MetalHelmholtzEos eos(argv[1], HelmholtzTableEos::Mixture::allow_documented_proxy);
    CompositionAtmosphereGrid grid(eos, argv[2],
        CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    double x, y, t, lg;
    int fields;
    while ((fields = std::scanf("%lf %lf %lf %lf", &x, &y, &t, &lg)) == 4) {
      auto c = solar_scaled(x, .02);
      c.basis = AbundanceBasis::baryon_mass;
      c.metal_inventory = MetalInventory::gs98;
      c.X[1] = y;
      c.X[2] -= y;
      const double g = std::pow(10., lg);
      if (!grid.covers(t, g, c)) {
        std::puts("{\"covered\":false}");
        continue;
      }
      const auto s = grid.eval(t, g, c);
      const auto r = grid.composition_response(t, g, c);
      std::printf("{\"covered\":true,\"T\":%.17g,\"Pgas\":%.17g,"
                  "\"P\":%.17g,\"rho\":%.17g,"
                  "\"logarithmic_thermal_derivatives\":[%.17g,%.17g,%.17g,%.17g],"
                  "\"composition_derivatives\":[%.17g,%.17g,%.17g,%.17g]}\n",
                  s.T, s.Pgas, s.P, s.rho, s.dlnT_dlnTeff, s.dlnT_dlng,
                  s.dlnP_dlnTeff, s.dlnP_dlng,
                  r.dlnT_dXH, r.dlnT_dX3, r.dlnP_dXH, r.dlnP_dX3);
    }
    if (fields != EOF)
      throw std::invalid_argument("expected four finite source coordinates");
  } catch (const std::exception &error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }
}

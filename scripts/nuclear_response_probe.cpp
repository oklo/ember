// Replay GS98 profile states and print the complete selected nuclear response.
#include "ember/nuclear.hpp"
#include <cstdio>
#include <exception>
int main() {
  using namespace ember;
  const PPChains nuclear(PPRates::solar_fusion_ii, PPScreening::salpeter_van_horn);
  double T, rho, X, Y3, Y4;
  try {
    while (std::scanf("%lf %lf %lf %lf %lf", &T, &rho, &X, &Y3, &Y4) == 5) {
      auto c = solar_scaled(X, .02);
      c.basis = AbundanceBasis::baryon_mass;
      c.metal_inventory = MetalInventory::gs98;
      c.X[1] = Y3;
      c.X[2] = Y4;
      const auto n = nuclear.composition_response(T, rho, c);
      std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g",
                  T, rho, X, Y3, Y4, n.state.eps, n.state.eps_neutrino,
                  n.state.dlneps_dlnT, n.state.dlneps_dlnRho);
      for (double value : n.state.dXdt) std::printf(" %.17g", value);
      for (double value : n.deps_dX) std::printf(" %.17g", value);
      for (const auto& row : n.d_dXdt_dX)
        for (double value : row) std::printf(" %.17g", value);
      std::putchar('\n');
    }
    return std::feof(stdin) ? 0 : 1;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }
}

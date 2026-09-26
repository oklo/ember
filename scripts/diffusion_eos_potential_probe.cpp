#include "ember/eos_mixture.hpp"
#include <iomanip>
#include <iostream>

// Read-only comparison of the production EOS's thermodynamic potential.
// Columns: X, Y3, T, rho, E/T-S, P, E, S, dP/dX, dE/dX.
int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("potential_probe family.dat");
    ember::MetalHelmholtzEos eos(
        argv[1], ember::HelmholtzTableEos::Mixture::allow_documented_proxy);
    double x, y, t, rho;
    std::cout << std::setprecision(17);
    while (std::cin >> x >> y >> t >> rho) {
      auto c = ember::solar_scaled(x, .02);
      c.basis = ember::AbundanceBasis::baryon_mass;
      c.metal_inventory = ember::MetalInventory::gs98;
      c.X[1] = y;
      c.X[2] -= y;
      const auto e = eos.eval_with_derivatives(t, rho, c);
      const auto d = eos.composition_response(t, rho, c);
      std::cout << x << ' ' << y << ' ' << t << ' ' << rho << ' '
                << e.state.E/t-e.state.S << ' ' << e.state.P << ' '
                << e.state.E << ' ' << e.state.S << ' '
                << d.dP[0] << ' ' << d.dE[0] << '\n';
    }
    if (!std::cin.eof()) throw std::invalid_argument("malformed input");
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

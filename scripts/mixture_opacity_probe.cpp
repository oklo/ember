#include "ember/opacity_mixture.hpp"
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("mixture_opacity_probe opacity_directory");
    ember::StellarMixtureOpacity opacity(argv[1]);
    double x, y, T, rho;
    std::cout << std::setprecision(17);
    while (std::cin >> x >> y >> T >> rho) {
      auto c = ember::solar_scaled(x, .02);
      c.basis = ember::AbundanceBasis::baryon_mass;
      c.metal_inventory = ember::MetalInventory::gs98;
      c.X[1] = y;
      c.X[2] -= y;
      const auto s = opacity.eval(T, rho, c);
      const auto range = opacity.density_range(T, c);
      if (!range) throw std::runtime_error("missing source density support");
      std::cout << x << ' ' << y << ' ' << T << ' ' << rho << ' ' << s.kappa << ' '
                << s.dlnk_dlnT << ' ' << s.dlnk_dlnRho << ' ' << s.dlnk_dX << ' '
                << s.dlnk_dY3 << ' ' << range->min << ' ' << range->max << '\n';
    }
    if (!std::cin.eof()) throw std::invalid_argument("malformed input");
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

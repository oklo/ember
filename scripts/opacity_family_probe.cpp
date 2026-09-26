// Query one source family, optionally mapping baryonic isotope abundances.
#include "ember/opacity_mixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("opacity_family_probe manifest");
    ember::MixtureOpacity source(argv[1]);
    ember::ElementalOpacity mapped(source);
    int baryonic;
    double x, y3, z, t, rho;
    std::cout << std::setprecision(17);
    while (std::cin >> baryonic >> x >> y3 >> z >> t >> rho) {
      if ((baryonic != 0 && baryonic != 1) || !(t>0) || !(rho>0))
        throw std::invalid_argument("invalid opacity query");
      auto c = ember::solar_scaled(x,z);
      c.basis = baryonic ? ember::AbundanceBasis::baryon_mass : ember::AbundanceBasis::atomic_mass;
      c.metal_inventory = ember::MetalInventory::gs98;
      c.X[1] = y3; c.X[2] -= y3;
      const ember::Opacity& opacity = baryonic ? static_cast<const ember::Opacity&>(mapped) : source;
      std::cout << "{\"query\":[" << baryonic << ',' << x << ',' << y3 << ',' << z << ',' << t << ',' << rho << ']';
      try {
        const auto s = opacity.eval(t,rho,c);
        const auto range = opacity.density_range(t,c).value();
        std::cout << ",\"covered\":true,\"kappa\":" << s.kappa
                  << ",\"dlnk_dlnT\":" << s.dlnk_dlnT
                  << ",\"dlnk_dlnrho\":" << s.dlnk_dlnRho
                  << ",\"dlnk_dX\":" << s.dlnk_dX
                  << ",\"density_range\":[" << range.min << ',' << range.max << "]}\n";
      } catch (const std::domain_error&) {
        std::cout << ",\"covered\":false}\n";
      }
    }
    if (!std::cin.eof()) throw std::invalid_argument("malformed opacity query");
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';return 1;
  }
}

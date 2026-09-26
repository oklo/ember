// Evaluate conductive opacity for an explicitly specified metal fraction.
// Input rows retain the saved-profile order: m r rho T L XH X3 X4.
#include "ember/conduction_table.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    if (argc != 4)
      throw std::invalid_argument("usage: conduction-composition TABLE PROFILE Z");
    const double Z = std::stod(argv[3]);
    if (!std::isfinite(Z) || Z < 0 || Z >= 1)
      throw std::invalid_argument("invalid metal fraction");
    ember::TabulatedConduction source(argv[1]);
    std::ifstream input(argv[2]);
    double m, r, rho, T, L, XH, X3, X4, checksum = 0;
    int count = 0;
    while (input >> m >> r >> rho >> T >> L >> XH >> X3 >> X4) {
      auto c = ember::solar_scaled(XH, Z);
      c.basis = ember::AbundanceBasis::baryon_mass;
      c.metal_inventory = ember::MetalInventory::gs98;
      c.X[1] = X3; c.X[2] = X4;
      const auto value = source.eval(T, rho, c);
      checksum += value.kappa;
      ++count;
      std::printf("%.17g %.17g %.17g %.17g %.17g\n", value.kappa,
          value.dlnk_dlnT, value.dlnk_dlnRho, value.dlnk_dX, value.dlnk_dY3);
    }
    if (!input.eof() || count == 0 || !std::isfinite(checksum))
      throw std::invalid_argument("invalid or empty conduction profile");
    std::printf("checksum %.17g\n", checksum);
  } catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what()); return 1;
  }
}

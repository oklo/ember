// Evaluate the selected conduction prescription at fixed material states.
#include "ember/conduction_table.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("usage: conduction-mixture TABLE < X Z Y3 T rho");
    ember::TabulatedConduction table(argv[1]);
    ember::HotConduction source(table);
    double x, z, y3, temperature, density;
    std::cout << std::setprecision(17);
    while (std::cin >> x >> z >> y3 >> temperature >> density) {
      auto composition = ember::solar_scaled(x, z);
      composition.basis = ember::AbundanceBasis::baryon_mass;
      composition.metal_inventory = ember::MetalInventory::gs98;
      composition.X[1] = y3;
      composition.X[2] = 1-x-z-y3;
      if (x < 0 || z < 0 || y3 < 0 || composition.X[2] < 0 || temperature <= 0 || density <= 0)
        throw std::invalid_argument("invalid material state");
      const auto value = source.eval(temperature, density, composition);
      if (!std::isfinite(value.kappa) || value.kappa <= 0)
        throw std::runtime_error("invalid conductive opacity");
      std::cout << value.kappa << '\n';
    }
    if (!std::cin.eof()) throw std::invalid_argument("malformed material query");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

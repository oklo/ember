// Query the actual fixed-Z opacity value, derivatives and supported density.
#include "ember/opacity_table.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("opacity_support_probe table");
    ember::TabulatedOpacity table(argv[1], "source opacity", ember::TabulatedOpacity::DensityAxis::logRho);
    std::cout << std::setprecision(17);
    double X, T, rho;
    while (std::cin >> X >> T >> rho) {
      if (!std::isfinite(X) || !(T > 0) || !(rho > 0))
        throw std::invalid_argument("invalid opacity query");
      const auto comp = ember::solar_scaled(X, table.metallicity());
      std::cout << "{\"query\":[" << X << ',' << T << ',' << rho << ']';
      try {
        const auto state = table.eval(T, rho, comp);
        const auto range = table.density_range(T, comp).value();
        std::cout << ",\"covered\":true,\"kappa\":" << state.kappa
                  << ",\"dlnk_dlnT\":" << state.dlnk_dlnT
                  << ",\"dlnk_dlnrho\":" << state.dlnk_dlnRho
                  << ",\"dlnk_dX\":" << state.dlnk_dX
                  << ",\"density_range\":[" << range.min << ',' << range.max << "]}\n";
      } catch (const std::domain_error&) {
        std::cout << ",\"covered\":false}\n";
      }
    }
    if (!std::cin.eof()) throw std::invalid_argument("malformed opacity query");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

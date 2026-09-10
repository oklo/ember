// SPDX-License-Identifier: GPL-3.0-or-later
// Optional offline LX-MIE source adapter; not linked into the Ember library.
#include "mie.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

int main() {
  try {
    double n{}, k{}, x{};
    std::cout << std::setprecision(17);
    while (std::cin >> n) {
      if (!(std::cin >> k >> x)) throw std::runtime_error("incomplete grain query");
      if (!std::isfinite(n) || !std::isfinite(k) || !std::isfinite(x)
          || n <= 0 || k < 0 || x < 1e-5 || x > 1e4)
        throw std::domain_error("invalid refractive index or size parameter");
      if (n == 1 && k == 0) {
        // No optical contrast. The phase asymmetry is immaterial when the
        // scattering cross section vanishes; represent it by zero.
        std::cout << "0 0 0 0\n" << std::flush;
        continue;
      }
      auto result = lxmie::Mie(std::complex<double>(n, -k), x);
      if (k == 0) {
        // A lossless material has exactly zero absorption. Check the raw
        // source cancellation before enforcing that analytic identity.
        const double tolerance = 256 * std::numeric_limits<double>::epsilon()
            * std::max(std::abs(result.q_ext), std::abs(result.q_sca));
        if (std::abs(result.q_abs) > tolerance)
          throw std::runtime_error("lossless Mie source violates energy conservation");
        result.q_abs = 0; result.q_ext = result.q_sca;
      }
      if (!std::isfinite(result.q_ext) || !std::isfinite(result.q_sca)
          || !std::isfinite(result.q_abs) || !std::isfinite(result.asymmetry_parameter)
          || result.q_ext < 0 || result.q_sca < 0 || result.q_abs < 0
          || std::abs(result.asymmetry_parameter) > 1)
        throw std::runtime_error("invalid grain optical response");
      std::cout << result.q_abs << ' ' << result.q_sca << ' '
                << result.q_ext << ' ' << result.asymmetry_parameter << '\n' << std::flush;
    }
    if (!std::cin.eof()) throw std::runtime_error("invalid grain query");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

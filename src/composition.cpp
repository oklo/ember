#include "ember/composition.hpp"
#include <cmath>

namespace ember {

// Asplund, Amarsi & Grevesse (2021), Table 2, photospheric log epsilon values
// (log n_X/n_H + 12).  Only the species we resolve are named; the rest of the
// metal budget is lumped into Zrest with an effective A = 20.
namespace {
constexpr double eps_C = 8.46, eps_N = 7.83, eps_O = 8.69;
// Isotopic splits: terrestrial/solar-wind C13/C12 and the CNO cycle's own
// starting point.  C13/C12 = 1/89 is the protosolar ratio.
constexpr double c13_over_c12 = 1.0 / 89.0;
}

Composition solar_scaled(double X, double Z) {
  Composition c{};
  // Number densities relative to hydrogen, converted to mass ratios.
  const double nC = std::pow(10.0, eps_C - 12.0);
  const double nN = std::pow(10.0, eps_N - 12.0);
  const double nO = std::pow(10.0, eps_O - 12.0);
  const double mC = nC * 12.011, mN = nN * 14.007, mO = nO * 15.999;
  // Total metal mass per hydrogen mass in the AAG21 mixture is dominated by
  // CNO plus Ne, Mg, Si, S, Fe; we resolve CNO and let the remainder ride in
  // Zrest, preserving the CNO *fraction* of Z rather than inventing one.
  constexpr double cno_fraction_of_Z = 0.72; // AAG21: C+N+O is ~72% of Z by mass
  const double mtot = mC + mN + mO;
  const double fC = mC / mtot, fN = mN / mtot, fO = mO / mtot;

  const double Zcno = Z * cno_fraction_of_Z;
  const double c12 = Zcno * fC / (1.0 + c13_over_c12);
  c[Species::C12]   = c12;
  c[Species::C13]   = c12 * c13_over_c12;
  c[Species::N14]   = Zcno * fN;
  c[Species::O16]   = Zcno * fO;
  c[Species::Zrest] = Z - Zcno;
  c[Species::H1]    = X;
  c[Species::He3]   = 0.0;
  c[Species::He4]   = 1.0 - X - Z;
  return c;
}

} // namespace ember

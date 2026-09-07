// Nuclear tests.  The pp chain has one number everybody knows - the Sun's
// central energy generation - and one behaviour that matters more for this
// code than that number: He3 must build up and then burn away, because that
// is what expands the core of a low-mass star over a trillion years.
#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <cstdio>
#include <string>

using namespace ember;
static int failures = 0;
static void check(bool ok, const std::string& what, double got, double want) {
  if (!ok) ++failures;
  std::printf("  [%s] %-50s got %-12.5g want %-12.5g\n",
              ok ? "PASS" : "FAIL", what.c_str(), got, want);
}

int main() {
  PPChains nuc;
  std::printf("ember nuclear - %s\n\n", nuc.name());

  // 1. Solar centre.  T = 1.55e7 K, rho = 150 g/cm^3, X = 0.35 gives an energy
  //    generation of order 10 erg/g/s; standard solar models put it near 17.
  {
    Composition c = solar_scaled(0.35, 0.014);
    c[Species::He4] = 1.0 - 0.35 - 0.014;
    const auto s = nuc.eval(1.55e7, 150.0, c);
    check(s.eps > 3.0 && s.eps < 60.0, "solar centre eps in the right decade", s.eps, 17.0);
    check(s.dlneps_dlnT > 3.5 && s.dlneps_dlnT < 6.0,
          "solar centre dln(eps)/dlnT ~ 4", s.dlneps_dlnT, 4.0);
  }

  // 2. The rate must switch off, exactly, when it is cold.
  {
    Composition c = solar_scaled(0.70, 0.014);
    const auto s = nuc.eval(1.0e5, 1.0, c);
    // Not exactly zero - the Gamow tail never is - but utterly negligible
    // against the 1e-7 erg/g/s that a cooling white dwarf radiates.
    check(s.eps < 1e-20, "burning is negligible at 1e5 K", s.eps, 0.0);
  }

  // 3. Baryon bookkeeping: the mass fractions must not drift.  Every reaction
  //    conserves nucleons, so the rates have to sum to zero.
  {
    Composition c = solar_scaled(0.50, 0.014);
    c[Species::He3] = 1e-3;
    const auto s = nuc.eval(1.2e7, 100.0, c);
    // Rest mass is *not* conserved: the mass fractions must fall at exactly
    // the rate the energy leaves, sum(dX/dt) = -eps/c^2.  Testing against zero
    // instead would hide a wrong branch ratio behind the mass defect.
    double sum = 0.0;
    for (double v : s.dXdt) sum += v;
    // eps is the *heating* rate, with the neutrino share already gone, so the
    // mass lost exceeds it.  The identity that must hold exactly is on the
    // total; here we check the weaker but still sharp statement that the two
    // agree to within the neutrino fraction, ~1%.
    const double want = -s.eps / (constants::c * constants::c);
    check(sum < 0.0 && std::abs(sum - want) <= 0.02 * std::abs(want),
          "mass defect: -sum dX/dt c^2 = eps + neutrinos", sum, want);
  }

  // 4. He3 grows where it is scarce and burns where it is abundant - the
  //    behaviour that drives the trillion-year core expansion.
  {
    Composition lo = solar_scaled(0.70, 0.014); lo[Species::He3] = 1e-8;
    Composition hi = solar_scaled(0.70, 0.014); hi[Species::He3] = 5e-2;
    const auto a = nuc.eval(8.0e6, 100.0, lo);
    const auto b = nuc.eval(8.0e6, 100.0, hi);
    const std::size_t i3 = static_cast<std::size_t>(Species::He3);
    check(a.dXdt[i3] > 0.0, "He3 accumulates when scarce", a.dXdt[i3], 1.0);
    check(b.dXdt[i3] < 0.0, "He3 is consumed when abundant", b.dXdt[i3], -1.0);
  }

  // 5. Hydrogen is consumed and helium made, always.
  {
    Composition c = solar_scaled(0.70, 0.014); c[Species::He3] = 1e-4;
    const auto s = nuc.eval(1.0e7, 100.0, c);
    check(s.dXdt[static_cast<std::size_t>(Species::H1)] < 0.0, "hydrogen decreases",
          s.dXdt[static_cast<std::size_t>(Species::H1)], -1.0);
    check(s.eps > 0.0, "energy generation is positive", s.eps, 1.0);
  }

  std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
              failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}

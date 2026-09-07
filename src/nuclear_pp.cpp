#include "ember/nuclear.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>

namespace ember {
using namespace constants;

namespace {

// Non-resonant thermonuclear reaction rate in the standard Gamow form.
// N_A <sigma v> = C * T9^{-2/3} exp(-tau) * (1 + corrections), with
// tau = 3 (E_G/4kT)^{1/3}.  Returning dln(rate)/dlnT alongside costs nothing
// and spares the solver a numerical derivative.
struct Rate { double v; double dlnv_dlnT; };

// Reduced mass factor and Gamow energy are folded into the coefficient and
// the exponent scale; both are taken from the Adelberger et al. (2011)
// compilation, expressed as fits in T9.
Rate pp(double T9) {                       // p(p,e+ nu)d  - the bottleneck
  const double t913 = std::cbrt(T9);
  const double t923 = t913 * t913;
  const double tau  = 3.381 / t913;
  // Adelberger+2011 eq. (2.32) form as used in standard compilations.
  const double f = 4.01e-15 / t923 * std::exp(-tau)
                 * (1.0 + 0.123 * t913 + 1.09 * t923 + 0.938 * T9);
  const double poly = 1.0 + 0.123 * t913 + 1.09 * t923 + 0.938 * T9;
  const double dpoly = (0.123 * t913 / 3.0 + 1.09 * t923 * 2.0 / 3.0 + 0.938 * T9) / poly;
  return {f, -2.0 / 3.0 + tau / 3.0 + dpoly};
}

Rate he3he3(double T9) {                   // He3(He3,2p)He4  - the ppI branch
  const double t913 = std::cbrt(T9);
  const double t923 = t913 * t913;
  const double tau  = 12.276 / t913;
  const double poly = 1.0 - 0.034 * t913 + 0.213 * t923 - 0.026 * T9;
  const double f = 6.04e10 / t923 * std::exp(-tau) * poly;
  const double dpoly = (-0.034 * t913 / 3.0 + 0.213 * t923 * 2.0 / 3.0 - 0.026 * T9) / poly;
  return {f, -2.0 / 3.0 + tau / 3.0 + dpoly};
}

Rate he3he4(double T9) {                   // He3(alpha,gamma)Be7 - ppII/ppIII
  const double t913 = std::cbrt(T9);
  const double t923 = t913 * t913;
  const double tau  = 12.826 / t913;
  const double f = 5.46e6 / t923 * std::exp(-tau);
  return {f, -2.0 / 3.0 + tau / 3.0};
}

// Salpeter weak screening.  Valid while the Coulomb coupling is small, which
// is the regime a hydrogen-burning low-mass star occupies; a cold dense
// remnant needs Chugunov et al. (2007) instead, which is why this lives behind
// the Nuclear interface rather than inside it.
double screen_weak(double T, double rho, const Composition& c, double z1, double z2) {
  const double zbar  = c.mu_elec_inv() / c.mu_ions_inv();
  const double z2bar = [&] {
    double s = 0.0;
    for (std::size_t i = 0; i < NSPEC; ++i)
      s += c.X[i] * nuclides[i].Z * nuclides[i].Z / nuclides[i].A;
    return s / c.mu_ions_inv();
  }();
  const double ne = c.mu_elec_inv() * NA * rho;
  const double lam = std::sqrt(4.0 * M_PI * (1.602176634e-10) * ne * (z2bar + zbar) / (kB * T));
  // H12 = z1 z2 e^2 / (kT) * kappa_D ; assembled in CGS below.
  const double e2 = 4.803204673e-10 * 4.803204673e-10;
  const double kD = std::sqrt(4.0 * M_PI * e2 * ne * (z2bar + zbar) / (kB * T));
  (void)lam;
  const double H = z1 * z2 * e2 * kD / (kB * T);
  return std::exp(std::min(H, 2.0));   // capped: beyond this the weak limit is void
}

} // namespace

NuclearState PPChains::eval(double T, double rho, const Composition& c) const {
  NuclearState s{};
  const double T9 = T * 1e-9;
  if (T9 < 1e-4) return s;                 // nothing happens; keep it exactly zero

  const double X  = c[Species::H1];
  const double Y3 = c[Species::He3];
  const double Y4 = c[Species::He4];
  const double A1 = nuclides[static_cast<std::size_t>(Species::H1)].A;
  const double A3 = nuclides[static_cast<std::size_t>(Species::He3)].A;
  const double A4 = nuclides[static_cast<std::size_t>(Species::He4)].A;

  const auto r_pp   = pp(T9);
  const auto r_33   = he3he3(T9);
  const auto r_34   = he3he4(T9);
  const double f_pp = screen_weak(T, rho, c, 1, 1);
  const double f_33 = screen_weak(T, rho, c, 2, 2);
  const double f_34 = screen_weak(T, rho, c, 2, 2);

  // Reactions per gram per second.  The 1/2 on identical-particle reactions is
  // the standard double-counting factor.
  const double n_pp = 0.5 * (X / A1) * (X / A1) * rho * r_pp.v * f_pp;
  const double n_33 = 0.5 * (Y3 / A3) * (Y3 / A3) * rho * r_33.v * f_33;
  const double n_34 =       (Y3 / A3) * (Y4 / A4) * rho * r_34.v * f_34;

  // Composition change first; the energy then follows from it.  Deriving the
  // release from the mass defect of the very nuclide masses the code carries,
  // rather than from a separately tabulated set of Q values, makes energy and
  // composition consistent by construction: they cannot drift apart, and a
  // mistaken branch ratio shows up as an energy error instead of hiding.
  auto& d = s.dXdt;
  const std::size_t iH1  = static_cast<std::size_t>(Species::H1);
  const std::size_t iHe3 = static_cast<std::size_t>(Species::He3);
  const std::size_t iHe4 = static_cast<std::size_t>(Species::He4);
  // Per p+p reaction three protons are consumed, not two: two make the
  // deuteron, and the fast d(p,gamma)He3 that follows takes a third.
  // ppI  : He3 + He3 -> He4 + 2p        returns two protons, makes one He4.
  // ppII : He3 + He4 + p -> 2 He4       consumes a proton, nets one He4.
  // Baryon number cancels identically; the mass defect does not, and leaves.
  d[iH1]  = (-3.0 * n_pp + 2.0 * n_33 - n_34) * A1;
  d[iHe3] = ( n_pp - 2.0 * n_33 - n_34) * A3;
  d[iHe4] = ( n_33 + n_34) * A4;

  double dm = 0.0;
  for (double v : d) dm += v;                 // negative: mass is lost
  const double eps_total = -dm * c_light * c_light;       // erg/g/s liberated

  // Neutrinos take their share straight out of the star.  pp emits 0.265 MeV
  // on average; the Be7 electron capture that opens ppII emits 0.861 MeV.
  constexpr double MeV = 1.602176634e-6;
  const double eps_nu = (n_pp * 0.265 + n_34 * 0.861) * MeV * NA;
  s.eps = eps_total - eps_nu;
  if (s.eps < 0.0) s.eps = 0.0;

  // Temperature and density derivatives.  Every channel is binary, so each
  // goes as rho^1; the temperature dependence is the Gamow exponent, weighted
  // by how much of the energy each channel is currently carrying.
  const double w_pp = n_pp * 6.67, w_33 = n_33 * 12.86, w_34 = n_34 * 19.0;
  const double wsum = w_pp + w_33 + w_34;
  if (wsum > 0.0) {
    s.dlneps_dlnT = (w_pp * r_pp.dlnv_dlnT + w_33 * r_33.dlnv_dlnT
                   + w_34 * r_34.dlnv_dlnT) / wsum;
    s.dlneps_dlnRho = 1.0;
  }
  return s;
}

} // namespace ember

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

// Classical Salpeter weak screening. Valid while the Coulomb coupling is
// small; intermediate coupling/degenerate electrons need further work. A cold dense
// remnant needs Chugunov et al. (2007) instead, which is why this lives behind
// the Nuclear interface rather than inside it.
struct Screening { double factor, dlnf_dlnT, dlnf_dlnRho; };
Screening screen_weak(double T, double rho, const Composition& comp, double z1, double z2) {
  // Classical electron + ion charge susceptibility: ne + sum(ni Zi^2).
  // Multiplying ion-averaged charges by ne instead of ni introduces an
  // erroneous extra mean ionic charge in mixtures containing helium.
  double charges=0;
  for(std::size_t i=0;i<NSPEC;++i) {
    const double z=nuclides[i].Z;
    charges+=comp.X[i]*(z*z+z)/comp.abundance_weight(i);
  }
  // H12 = z1 z2 e^2 / (kT) * kappa_D ; assembled in CGS below.
  const double e2 = 4.803204673e-10 * 4.803204673e-10;
  const double kD = std::sqrt(4.0 * M_PI * e2 * rho * NA * charges / (kB * T));
  const double H = z1 * z2 * e2 * kD / (kB * T);
  // Differentiate the implemented cap as well: above it the factor is
  // constant, while below it H is proportional to rho^(1/2) T^(-3/2).
  return {std::exp(std::min(H, 2.0)), H < 2.0 ? -1.5 * H : 0.0,
          H < 2.0 ? 0.5 * H : 0.0};
}

} // namespace

NuclearState PPChains::eval(double T, double rho, const Composition& c) const {
  NuclearState s{};
  const double T9 = T * 1e-9;
  if (T9 < 1e-4) return s;                 // nothing happens; keep it exactly zero

  const double X  = c[Species::H1];
  const double Y3 = c[Species::He3];
  const double Y4 = c[Species::He4];
  const double A1 = c.abundance_weight(0);
  const double A3 = c.abundance_weight(1);
  const double A4 = c.abundance_weight(2);

  const auto r_pp   = pp(T9);
  const auto r_33   = he3he3(T9);
  const auto r_34   = he3he4(T9);
  const auto f_pp = screen_weak(T, rho, c, 1, 1);
  const auto f_33 = screen_weak(T, rho, c, 2, 2);
  const auto f_34 = screen_weak(T, rho, c, 2, 2);

  // Reactions per gram per second.  The 1/2 on identical-particle reactions is
  // the standard double-counting factor.
  const double n_pp = 0.5 * (X / A1) * (X / A1) * rho * r_pp.v * f_pp.factor;
  const double n_33 = 0.5 * (Y3 / A3) * (Y3 / A3) * rho * r_33.v * f_33.factor;
  const double n_34 =       (Y3 / A3) * (Y4 / A4) * rho * r_34.v * f_34.factor;

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
  if(c.basis == AbundanceBasis::atomic_mass) {
    for(double v:d) dm+=v; // retain the legacy static convention exactly
  } else {
    for(std::size_t i=0;i<NSPEC;++i) dm+=d[i]*nuclides[i].A/mass_numbers[i];
  }
  const double eps_total = -dm * c_light * c_light;       // erg/g/s liberated

  // Neutrinos take their share straight out of the star.  pp emits 0.265 MeV
  // on average; the Be7 electron capture that opens ppII emits 0.861 MeV.
  constexpr double MeV = 1.602176634e-6;
  const double eps_nu = (n_pp * 0.265 + n_34 * 0.861) * MeV * NA;
  s.eps_neutrino = eps_nu;
  s.eps = eps_total - eps_nu;
  if (s.eps < 0.0) s.eps = 0.0;

  // Differentiate the same mass-defect heating rate, including the escaping
  // neutrinos and the density/temperature dependence of screening.
  const double c2 = c_light * c_light;
  const double m1=nuclides[0].A,m3=nuclides[1].A,m4=nuclides[2].A;
  const double w_pp = n_pp * ((3.0 * m1 - m3) * c2 - 0.265 * MeV * NA);
  const double w_33 = n_33 * ((2.0 * m3 - 2.0 * m1 - m4) * c2);
  const double w_34 = n_34 * ((m3 + m1 - m4) * c2 - 0.861 * MeV * NA);
  if (s.eps > 0.0) {
    s.dlneps_dlnT = (w_pp * (r_pp.dlnv_dlnT + f_pp.dlnf_dlnT)
                   + w_33 * (r_33.dlnv_dlnT + f_33.dlnf_dlnT)
                   + w_34 * (r_34.dlnv_dlnT + f_34.dlnf_dlnT)) / s.eps;
    s.dlneps_dlnRho = (w_pp * (1.0 + f_pp.dlnf_dlnRho)
                     + w_33 * (1.0 + f_33.dlnf_dlnRho)
                     + w_34 * (1.0 + f_34.dlnf_dlnRho)) / s.eps;
  }
  return s;
}

NuclearResponse PPChains::composition_response(double T,double rho,const Composition& comp) const {
  NuclearResponse result;result.state=eval(T,rho,comp);
  if(T<1e5) return result;
  const double w1=comp.abundance_weight(0),w3=comp.abundance_weight(1),w4=comp.abundance_weight(2);
  const double h1=comp.X[0]/w1,h3=comp.X[1]/w3,h4=comp.X[2]/w4;
  const auto s1=screen_weak(T,rho,comp,1,1),s2=screen_weak(T,rho,comp,2,2);
  const double c1=.5*rho*pp(T*1e-9).v*s1.factor;
  const double c2=.5*rho*he3he3(T*1e-9).v*s2.factor;
  const double c3=rho*he3he4(T*1e-9).v*s2.factor;
  double charges=0;
  for(std::size_t j=0;j<NSPEC;++j) {
    const double y=comp.X[j]/comp.abundance_weight(j),z=nuclides[j].Z;
    charges+=(z*z+z)*y;
  }
  const double m1=nuclides[0].A,m3=nuclides[1].A,m4=nuclides[2].A;
  constexpr double mev=1.602176634e-6;
  const double q1=(3*m1-m3)*c_light*c_light-.265*mev*NA;
  const double q2=(2*m3-2*m1-m4)*c_light*c_light;
  const double q3=(m3+m1-m4)*c_light*c_light-.861*mev*NA;
  for(std::size_t j=0;j<NSPEC;++j) {
    const double z=nuclides[j].Z,w=comp.abundance_weight(j);
    const double derivative=.5*(z*z+z)/(w*charges);
    const double df1=s1.dlnf_dlnRho>0?std::log(s1.factor)*derivative:0;
    const double df2=s2.dlnf_dlnRho>0?std::log(s2.factor)*derivative:0;
    const double n1=c1*h1*h1*df1+(j==0?2*c1*h1/w1:0);
    const double n2=c2*h3*h3*df2+(j==1?2*c2*h3/w3:0);
    const double n3=c3*h3*h4*df2+(j==1?c3*h4/w3:0)+(j==2?c3*h3/w4:0);
    result.d_dXdt_dX[0][j]=(-3*n1+2*n2-n3)*w1;
    result.d_dXdt_dX[1][j]=(n1-2*n2-n3)*w3;
    result.d_dXdt_dX[2][j]=(n2+n3)*w4;
    result.deps_dX[j]=q1*n1+q2*n2+q3*n3;
  }
  return result;
}

} // namespace ember

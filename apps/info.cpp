#include "ember/eos.hpp"
#include "ember/constants.hpp"
#include <cstdio>
int main() {
  using namespace ember;
  std::printf("ember - stellar evolution for the lowest-mass stars\n");
  std::printf("  Msun = %.6e g   Rsun = %.4e cm   Lsun = %.4e erg/s\n",
              constants::Msun, constants::Rsun, constants::Lsun);
  IdealEos eos;
  const auto comp = solar_scaled(0.70, 0.014);
  const auto s = eos.eval(1.5e7, 100.0, comp);
  std::printf("  eos(%s) at T=1.5e7 K, rho=100: P=%.4e  grad_ad=%.4f  Gamma1=%.4f\n",
              eos.name(), s.P, s.grad_ad, s.Gamma1);
  return 0;
}

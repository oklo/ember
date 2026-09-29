#pragma once
#include "ember/eos_helmholtz.hpp"
#include <array>
#include <filesystem>
#include <vector>

namespace ember {
// Material H/He free energy from pure, thermodynamically consistent potentials
// combined at common T and pressure. Nuclear mixing entropy and the helium
// isotope correction remain analytic in the caller. This component supplies
// neither metals nor a transition to another EOS; it is not an Eos selection.
// Input densities and energies use the conserved baryonic gram.
class AdditiveVolumePotential {
 public:
  explicit AdditiveVolumePotential(const std::filesystem::path&);
  // Residual Phi=F/T before R sum(n_i ln n_i) and isotope entropy are added.
  // Channels: value, H/He3 first derivatives, HH/H-He3/He3-He3 second
  // derivatives. Thermal jet order + composition order never exceeds three.
  std::array<HelmholtzJet,6> residual_jets(double T,double rho,double X,double Y3,
                                        bool composition=true)const;
  Eos::DensityRange density_range(double T,double X,double Y3)const;
  double minimum_temperature()const;
 private:
  struct Pure {
    std::vector<double> t,r;
    std::vector<std::array<double,16>> cells;
    std::vector<std::array<double,2>> density_bounds;
  };
  Pure h_,he_;
  static std::array<double,2> bounds(const Pure&,double logT);
  template<class J>static J phi(const Pure&,const J&,const J&,unsigned dr=0);
  static double density(const Pure&,double logT,double logP);
  template<class J>static J density_jet(const Pure&,const J&,const J&);
  template<class J>J residual(const J& logT,const J& logRho,const J& X,const J& Y3)const;
};
} // namespace ember

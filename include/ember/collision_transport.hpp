#pragma once
#include "ember/material_transport.hpp"
#include <memory>
#include <string>

namespace ember {
template<std::size_t Species> struct BasicCollisionTransportResponse {
  // Response to (-chemical_acceleration/T, -energy_scale grad(T)/T^2).
  // Last flux is reduced heat / energy_scale; EOS enthalpy is not included.
  std::array<std::array<double,Species+1>,Species+1> mobility{};
  std::array<double,Species> transport_enthalpy{};
  std::array<bool,Species> active_species{};
  double conductivity{}, energy_scale{}, eta{}, electron_density{}, b_thermal{};
  double maximum_solve_backward_error{};
};
template<std::size_t Species> struct BasicCollisionTransportPartial {
  std::array<std::array<double,Species+1>,Species+1> mobility{};
  std::array<double,Species> transport_enthalpy{};
  double conductivity{},energy_scale{},eta{},electron_density{},b_thermal{};
};
template<std::size_t Species> struct BasicCollisionTransportDerivatives {
  BasicCollisionTransportResponse<Species> value;
  // Coordinates: ln T, ln rho, X, Y3, Z, ln screening_length.
  std::array<BasicCollisionTransportPartial<Species>,6> partials{};
  // An absent species stays removed. Its composition derivative is undefined
  // at the boundary; callers must not use the zero storage as a derivative.
  std::array<bool,6> defined{};
};
using CollisionTransportResponse=BasicCollisionTransportResponse<2>;
using CollisionTransportPartial=BasicCollisionTransportPartial<2>;
using CollisionTransportDerivatives=BasicCollisionTransportDerivatives<2>;
using BulkMetalCollisionResponse=BasicCollisionTransportResponse<3>;
using BulkMetalCollisionDerivatives=BasicCollisionTransportDerivatives<3>;
struct CollisionScreeningDerivatives {
  double value{};
  // Coordinates: ln T, ln rho, X, Y3, Z, ln electron_stiffness.
  std::array<double,6> partials{};
  std::array<bool,6> defined{};
};

// Conditional kinetic model: classical repulsive screened ion collisions,
// nonrelativistic Born/Pauli electrons, ten electron energy modes, ion heat
// variables and leading Brownian ion heat relaxation. Metals are fully
// stripped, stationary GS98 carriers. Input fractions are baryonic.
// This is not an ionization or screening prescription for a cool atmosphere.
class ScreenedCollisionTransport {
 public:
  explicit ScreenedCollisionTransport(const std::string& exported_table);
  CollisionTransportResponse eval(double temperature, double density,
      double X, double Y3, double Z, double screening_length_cm) const;
  CollisionTransportDerivatives derivatives(double temperature,double density,
      double X,double Y3,double Z,double screening_length_cm) const;
  // Additional mass direction: all GS98 metals move at one shared velocity.
  // He4 and electron velocities enforce zero mass flow and electric current.
  // Metal heat variables and collision charges remain distinct. This is an
  // explicit grouping approximation, not independently settling elements.
  BulkMetalCollisionResponse bulk_metal_eval(double temperature,double density,
      double X,double Y3,double Z,double screening_length_cm) const;
  BulkMetalCollisionDerivatives bulk_metal_derivatives(double temperature,double density,
      double X,double Y3,double Z,double screening_length_cm) const;
  // Optional explicit effective static-screening calculation. The supplied
  // electron stiffness is (dP_e/dln n_e)/n_e at fixed T. No EOS is inferred.
  static double screening_length(double temperature, double density,
      double X, double Y3, double Z, double electron_stiffness_erg,
      bool include_ions);
  static CollisionScreeningDerivatives screening_derivatives(double temperature,double density,
      double X,double Y3,double Z,double electron_stiffness_erg,bool include_ions);
 private:
  struct Data;
  std::shared_ptr<const Data> data_;
};
} // namespace ember

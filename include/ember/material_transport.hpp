#pragma once
#include <array>
#include <functional>
#include <span>
#include <vector>

namespace ember {
// At most two independent mass fractions followed by scaled internal energy.
// Only the first `components` entries are used; unused entries remain zero.
using MaterialVector = std::array<double, 3>;
using MaterialMatrix = std::array<MaterialVector, 3>;

struct MaterialPoint {
  MaterialVector conserved{}, potential{};
  MaterialMatrix capacity{}, primitive_from_conserved{};
  double entropy{};
};
using MaterialThermodynamics = std::function<MaterialPoint(
    std::size_t, const MaterialVector&, bool derivatives)>;
using MaterialConductance = std::function<std::vector<MaterialMatrix>(
    std::span<const MaterialVector>)>;

struct MaterialTransportOptions {
  std::size_t components{3}, max_iterations{60}, max_backtracks{50};
  double tolerance{2e-11}, composition_sum_limit{1.0};
};
struct MaterialTransportResult {
  std::vector<MaterialVector> primitives, conserved, face_flux;
  std::vector<double> residual_history;
  MaterialVector integrated_conservation_error{};
  std::size_t iterations{};
  double residual{}, entropy_change{}, backward_euler_entropy_bound{};
};

// Positive block solve for mass capacity plus a closed-chain diffusion operator.
// Inputs are symmetric positive definite, except exactly zero face matrices.
// Work/storage are linear in cells. No diagonal floor or eigenvalue clipping.
MaterialMatrix material_positive_inverse(const MaterialMatrix&, std::size_t components);
std::vector<MaterialVector> solve_material_chain(std::span<const MaterialMatrix> capacities,
    std::span<const MaterialMatrix> conductances, std::span<const MaterialVector> rhs,
    std::size_t components);

// Backward Euler on fixed cell masses, with closed boundaries. Primitives are
// independent fractions and ln(T); U=(X...,E/e0), W=(mu/T...,-e0/T)/s0.
// Conductances already include s0, face area/distance, and any mass units.
// One shared face flux gives exact algebraic conservation. The coefficient
// callback is reevaluated in line searches; its derivative is frozen in the
// Newton direction. No general convergence guarantee for arbitrary callbacks.
// Old states may contain zeros: derivatives=false needs only U and entropy.
// Guesses must be interior; they never replace old conserved abundances.
// Remove globally absent species before calling. No burning, pressure work,
// convection, radiation flux, or hydrostatic update is performed here.
MaterialTransportResult transport_material(std::span<const MaterialVector> old_primitives,
    std::span<const MaterialVector> initial_guess, std::span<const double> mass,
    double dt, const MaterialThermodynamics&, const MaterialConductance&,
    const MaterialTransportOptions& = {});

class SmoothMetalHelmholtzEos;
struct Composition;
// Fixed-density native EOS callback for [XH,X3,ln(T)] with He4 as reference.
// Metal inventory remains fixed. E,S include local radiation thermodynamics.
// Chemical derivatives are not requested at old, possibly zero, abundances.
MaterialPoint native_material_point(const SmoothMetalHelmholtzEos&, double T, double rho,
    const Composition&, double energy_scale, double entropy_scale, bool derivatives);
} // namespace ember

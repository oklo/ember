#pragma once
#include "ember/evolution.hpp"
#include "ember/evolution_state.hpp"
#include <array>
#include <functional>
#include <limits>
#include <string_view>

namespace ember {
struct EvolutionAudit {
  bool pass{};
  double maximum_species_error{}, mass_error_surface{}, source_error{}, first_law{};
};

struct EvolutionControlOptions {
  EvolutionOptions step{};
  double target_age{};
  double maximum_dt{};
  double structure_tolerance{1e-4}, species_tolerance{1e-8}, energy_tolerance{.005};
  double safety_factor{.8}, minimum_growth{.5}, maximum_growth{1.5};
  double minimum_error{1e-8}, rejection_factor{.5};
  std::size_t maximum_consecutive_rejections{8};
  std::size_t maximum_steps{500};
  double maximum_cpu_seconds{std::numeric_limits<double>::infinity()};
  bool audit_failure_is_fatal{true};
};

struct EvolutionAttempt {
  double start_age{}, dt{}, error_norm{1e30};
  bool converged{}, audit_pass{}, accepted{};
  std::string message;
  // Full step, first half and second half, respectively.
  std::array<EvolutionAudit, 3> audits{};
};

// The controller owns time stepping. Physics selection, validity checks and
// output are supplied by the caller; none of them depends on a command line.
// Callbacks may reject an unsupported state by throwing. The returned state
// then remains the last accepted state, with the failure in the result.
struct EvolutionControlHooks {
  std::function<const Physics&(const Model&)> physics;
  // Optional state-dependent selection before a solve. The second half
  // interval uses its own starting model, including any exhausted isotope.
  std::function<void(const Model&, EvolutionOptions&)> configure_step;
  std::function<double(const Composition&, const Composition&)> species_difference;
  std::function<EvolutionAudit(const Model&, const EvolutionStep&, double)> audit;
  std::function<void(const Model&, std::span<const std::array<double, 3>>)> assess;
  std::function<double()> cpu_seconds;
  std::function<bool(std::string_view)> terminal_failure;
  std::function<void(const EvolutionAttempt&)> attempted;
  // Called at the starting state (dt=error=0), then after each acceptance.
  std::function<void(const EvolutionState&, double dt, double error)> accepted;
};

struct EvolutionControlResult {
  bool requested_age_reached{};
  std::size_t accepted_this_invocation{};
  std::string stop_reason;
};

EvolutionControlResult evolve(EvolutionState&, const Atmosphere&,
    const EvolutionControlOptions&, const EvolutionControlHooks&);
} // namespace ember

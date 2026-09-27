#include "ember/controller.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
void check_options(const EvolutionControlOptions& options, const EvolutionControlHooks& hooks) {
  if (!hooks.physics || !hooks.species_difference || !hooks.audit || !hooks.assess
      || !hooks.cpu_seconds)
    throw std::invalid_argument("evolution controller requires physics, error, audit, assessment and clock callbacks");
  if (!(std::isfinite(options.target_age) && options.target_age > 0
        && std::isfinite(options.maximum_dt) && options.maximum_dt > 0
        && options.structure_tolerance > 0 && options.species_tolerance > 0
        && options.energy_tolerance > 0 && options.safety_factor > 0
        && options.minimum_growth > 0 && options.maximum_growth >= options.minimum_growth
        && options.minimum_error > 0 && options.rejection_factor > 0 && options.rejection_factor < 1
        && options.maximum_consecutive_rejections > 0 && options.maximum_steps > 0
        && options.maximum_cpu_seconds > 0))
    throw std::invalid_argument("invalid evolution controller options");
}
} // namespace

EvolutionControlResult evolve(EvolutionState& state, const Atmosphere& atmosphere,
    const EvolutionControlOptions& options, const EvolutionControlHooks& hooks) {
  check_options(options, hooks);
  EvolutionControlResult result;
  std::size_t failed = 0;
  try {
    hooks.assess(state.model, state.metal_heat_rates);
    if (hooks.accepted) hooks.accepted(state, 0, 0);
    for (;;) {
      const double target = options.target_age;
      if (state.model.age >= target
          || target - state.model.age <= 8 * (std::nextafter(target, INFINITY) - target)) {
        result.stop_reason = "requested age";
        result.requested_age_reached = true;
        break;
      }
      if (result.accepted_this_invocation >= options.maximum_steps) {
        result.stop_reason = "accepted-step budget";
        break;
      }
      if (hooks.cpu_seconds() >= options.maximum_cpu_seconds) {
        result.stop_reason = "CPU budget";
        break;
      }
      const double ds = std::min({state.next_dt, options.maximum_dt, target - state.model.age});
      if (state.model.age + ds == state.model.age || ds <= 0)
        throw std::runtime_error("timestep is below clock resolution");
      auto selected = options.step;
      selected.previous_metal_heat_rates = state.metal_heat_rates;
      const auto& physics = hooks.physics(state.model);
      const auto full = evolve_step(state.model, physics, atmosphere, ds, selected);
      const auto h1 = full.converged ? evolve_step(state.model, physics, atmosphere, ds / 2, selected) : full;
      auto second_options = selected;
      second_options.previous_metal_heat_rates = h1.total_metal_species_rates;
      const auto h2 = h1.converged ? evolve_step(h1.model, physics, atmosphere, ds / 2, second_options) : h1;
      EvolutionAttempt attempt;
      attempt.start_age = state.model.age;
      attempt.dt = ds;
      attempt.converged = full.converged && h1.converged && h2.converged;
      attempt.audits = {hooks.audit(state.model, full, ds), hooks.audit(state.model, h1, ds / 2),
                        hooks.audit(h1.model, h2, ds / 2)};
      attempt.audit_pass = std::all_of(attempt.audits.begin(), attempt.audits.end(),
                                      [](const auto& audit) { return audit.pass; });
      if (attempt.converged) {
        double& error = attempt.error_norm;
        error = 0;
        for (std::size_t i = 0; i < state.model.size(); ++i) {
          for (const auto v : {Var::lnr, Var::lnrho, Var::lnT})
            error = std::max(error, std::abs(full.model.y[i][v] - h2.model.y[i][v]) / options.structure_tolerance);
          error = std::max(error, hooks.species_difference(full.model.comp[i], h2.model.comp[i])
                                 / options.species_tolerance);
        }
        const double half_power = .5 * (h1.nuclear_luminosity + h2.nuclear_luminosity);
        const double power_scale = std::max(std::abs(half_power),
            .5 * (std::abs(h1.model.y.back().L) + std::abs(h2.model.y.back().L)));
        error = std::max(error, std::abs(full.nuclear_luminosity - half_power) / (options.energy_tolerance * power_scale));
        error = std::max(error, std::abs(full.model.y.back().L / h2.model.y.back().L - 1) / (4 * options.structure_tolerance));
        hooks.assess(full.model, full.total_metal_species_rates);
        hooks.assess(h1.model, h1.total_metal_species_rates);
        hooks.assess(h2.model, h2.total_metal_species_rates);
      }
      attempt.accepted = attempt.converged && attempt.audit_pass && attempt.error_norm <= 1;
      attempt.message = !full.converged ? full.message : !h1.converged ? h1.message : h2.message;
      if (attempt.converged && !attempt.audit_pass)
        attempt.message = "physical inventory/energy audit failed";
      if (hooks.attempted) hooks.attempted(attempt);
      if (attempt.accepted) {
        state.model = h2.model;
        state.metal_heat_rates = h2.total_metal_species_rates;
        ++state.accepted;
        ++result.accepted_this_invocation;
        failed = 0;
        state.next_dt = std::min(options.maximum_dt,
            ds * std::clamp(options.safety_factor / std::sqrt(std::max(attempt.error_norm, options.minimum_error)),
                            options.minimum_growth, options.maximum_growth));
        if (hooks.accepted) hooks.accepted(state, ds, attempt.error_norm);
      } else {
        ++state.rejected;
        state.next_dt = ds * options.rejection_factor;
        if (attempt.converged && !attempt.audit_pass && options.audit_failure_is_fatal)
          throw std::runtime_error("physical inventory/energy audit failed");
        if (hooks.terminal_failure && hooks.terminal_failure(attempt.message))
          throw std::runtime_error(attempt.message);
        if (++failed >= options.maximum_consecutive_rejections)
          throw std::runtime_error("repeated step rejection: " + attempt.message);
      }
    }
  } catch (const std::exception& error) {
    result.stop_reason = error.what();
  }
  return result;
}
} // namespace ember

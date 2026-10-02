#include "ember/controller.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace ember {
namespace {
void check_options(const EvolutionControlOptions& options, const EvolutionControlHooks& hooks) {
  if (!hooks.physics || !hooks.species_difference || !hooks.audit || !hooks.assess
      || !hooks.cpu_seconds)
    throw std::invalid_argument("evolution controller requires physics, error, audit, assessment and clock callbacks");
  if (options.richardson_extrapolation && !hooks.assess_extrapolated)
    throw std::invalid_argument("Richardson extrapolation requires a physical assessment callback");
  if (!(std::isfinite(options.target_age) && options.target_age > 0
        && std::isfinite(options.maximum_dt) && options.maximum_dt > 0
        && options.structure_tolerance > 0 && options.species_tolerance > 0
        && options.energy_tolerance > 0 && options.safety_factor > 0
        && options.minimum_growth > 0 && options.maximum_growth >= options.minimum_growth
        && options.minimum_error > 0 && options.rejection_factor > 0 && options.rejection_factor < 1
        && options.maximum_consecutive_rejections > 0 && options.maximum_steps > 0
        && options.maximum_cpu_seconds > 0))
    throw std::invalid_argument("invalid evolution controller options");
  if (options.relaxed_intervals > 16 || !std::isfinite(options.relaxed_accuracy_factor)
      || options.relaxed_accuracy_factor < 1 || options.relaxed_accuracy_factor > 10
      || (!options.relaxed_intervals && options.relaxed_accuracy_factor != 1))
    throw std::invalid_argument("invalid periodic accuracy settings");
}
EvolutionOptions solve_accuracy(EvolutionOptions selected, double factor) {
  selected.relaxation.residual_tolerance *= factor;
  selected.relaxation.correction_tolerance *= factor;
  selected.coupling_stop_tolerance *= factor;
  selected.verification_residual_tolerance *= factor;
  selected.verification_correction_tolerance *= factor;
  return selected;
}
bool rapid_change(const Model& start, const EvolutionStep& trial,
                  const std::optional<double>& previous_power) {
  if (!trial.converged) return true;
  const double luminosity = start.y.back().L;
  if (!(luminosity > 0 && trial.model.y.back().L > 0)
      || std::abs(std::log(trial.model.y.back().L / luminosity)) > .10) return true;
  if (previous_power && std::abs(trial.nuclear_luminosity - *previous_power)
      > .10 * std::max(std::abs(*previous_power), luminosity)) return true;
  for (std::size_t i=0;i<start.size();++i)
    if (std::abs(trial.model.y[i].lnT-start.y[i].lnT) > .05
        || std::abs(trial.model.y[i].lnrho-start.y[i].lnrho) > .10) return true;
  return false;
}
// 2*h2 - full, or nothing when the combination is not admissible. The cheap
// count/mass prefilter is only a shortcut; the physical partition comparison
// is made by the assessment callback.
std::optional<Model> richardson(const EvolutionStep& full, const EvolutionStep& h1, const EvolutionStep& h2) {
  if (full.mixed_regions != h2.mixed_regions || h1.mixed_regions != h2.mixed_regions
      || full.convective_mass_fraction != h2.convective_mass_fraction
      || h1.convective_mass_fraction != h2.convective_mass_fraction)
    return std::nullopt;
  if (full.model.m != h2.model.m || full.model.M != h2.model.M || full.model.envelope_mass != h2.model.envelope_mass
      || full.model.comp.size() != h2.model.comp.size()) return std::nullopt;
  Model out = h2.model;
  for (std::size_t i = 0; i < out.size(); ++i) {
    for (const auto v : {Var::lnr, Var::lnrho, Var::lnT, Var::L})
      out.y[i][v] = h2.model.y[i][v] + (h2.model.y[i][v] - full.model.y[i][v]);
    for (const auto v : {Var::lnr, Var::lnrho, Var::lnT, Var::L})
      if (!std::isfinite(out.y[i][v])) return std::nullopt;
    auto& c = out.comp[i];const auto& a = h2.model.comp[i];const auto& b = full.model.comp[i];
    if (a.cn_molality.has_value() != b.cn_molality.has_value()) return std::nullopt;
    for (std::size_t j = 0; j < NSPEC; ++j) {
      c.X[j] = 2 * a.X[j] - b.X[j];
      if (!(c.X[j] >= 0 && c.X[j] <= 1) || (a.X[j] > 0) != (c.X[j] > 0)) return std::nullopt;
    }
    if (c.cn_molality)
      for (std::size_t j = 0; j < 3; ++j) {
        (*c.cn_molality)[j] = 2 * (*a.cn_molality)[j] - (*b.cn_molality)[j];
        if (!std::isfinite((*c.cn_molality)[j]) || !((*c.cn_molality)[j] >= 0)) return std::nullopt;
      }
  }
  return out;
}
} // namespace

EvolutionControlResult evolve(EvolutionState& state, const Atmosphere& atmosphere,
    const EvolutionControlOptions& options, const EvolutionControlHooks& hooks) {
  check_options(options, hooks);
  EvolutionControlResult result;
  std::size_t failed = 0;
  std::optional<Model> preceding;
  std::optional<double> previous_power;
  double preceding_dt = 0;
  const auto solve = [&](const Model& start, const Physics& physics, double dt,
                         EvolutionOptions selected, const Model* before, double ratio) {
    std::optional<Model> guess;
    selected.initial_structure_guess = nullptr;
    if (options.predict_structure && before && ratio > 0 && std::isfinite(ratio)
        && before->M == start.M && before->envelope_mass == start.envelope_mass && before->m == start.m && before->size() == start.size()
        && before->luminosity_grid == start.luminosity_grid) {
      bool smooth_history = true;
      for (std::size_t i = 0; i < start.size(); ++i)
        if (hooks.species_difference(start.comp[i], before->comp[i]) > selected.max_abundance_change)
          smooth_history = false;
      // Instantaneous mixing changes are not proportional to elapsed time.
      // Do not extrapolate their structural response into the next step.
      if (smooth_history) guess = start;
      bool finite = smooth_history;
      for (std::size_t i = 0; smooth_history && i < start.size(); ++i)
        for (std::size_t k = 0; k < NVAR; ++k) {
          const auto v = static_cast<Var>(k);
          guess->y[i][v] += ratio * (start.y[i][v] - before->y[i][v]);
          finite = finite && std::isfinite(guess->y[i][v]);
        }
      if (finite) selected.initial_structure_guess = &*guess;
    }
    auto step = evolve_step(start, physics, atmosphere, dt, selected);
    if (!step.converged && selected.initial_structure_guess) {
      // A poor prediction must not force a smaller physical timestep.
      selected.initial_structure_guess = nullptr;
      step = evolve_step(start, physics, atmosphere, dt, selected);
    }
    return step;
  };
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
      double accuracy_factor = options.relaxed_intervals && result.accepted_this_invocation
          && !failed && state.accepted % (options.relaxed_intervals + 1) != 0
          && ds < target-state.model.age ? options.relaxed_accuracy_factor : 1.;
      auto selected = solve_accuracy(options.step,accuracy_factor);
      if (hooks.configure_step) hooks.configure_step(state.model, selected);
      selected.previous_metal_heat_rates = state.metal_heat_rates;
      const auto& physics = hooks.physics(state.model);
      const Model* before = preceding ? &*preceding : nullptr;
      const double ratio = preceding_dt > 0 ? ds / preceding_dt : 0;
      auto full = solve(state.model, physics, ds, selected, before, ratio);
      bool tightened_after_trial = false;
      if (accuracy_factor > 1 && rapid_change(state.model,full,previous_power)) {
        accuracy_factor = 1; tightened_after_trial = true;
        selected = options.step;
        if (hooks.configure_step) hooks.configure_step(state.model,selected);
        selected.previous_metal_heat_rates = state.metal_heat_rates;
        full = solve(state.model,physics,ds,selected,before,ratio);
      }
      const auto h1 = full.converged ? solve(state.model, physics, ds / 2, selected, before, ratio / 2) : full;
      auto second_options = solve_accuracy(options.step,accuracy_factor);
      if (h1.converged && hooks.configure_step) hooks.configure_step(h1.model, second_options);
      second_options.previous_metal_heat_rates = h1.total_metal_species_rates;
      const auto h2 = h1.converged ? solve(h1.model, hooks.physics(h1.model), ds / 2,
                                         second_options, &state.model, 1.) : h1;
      EvolutionAttempt attempt;
      attempt.accuracy_factor = accuracy_factor;
      attempt.tightened_after_trial = tightened_after_trial;
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
            error = std::max(error, std::abs(full.model.y[i][v] - h2.model.y[i][v]) / (accuracy_factor * options.structure_tolerance));
          error = std::max(error, hooks.species_difference(full.model.comp[i], h2.model.comp[i])
                                 / (accuracy_factor * options.species_tolerance));
        }
        const double half_power = .5 * (h1.nuclear_luminosity + h2.nuclear_luminosity);
        const double power_scale = std::max(std::abs(half_power),
            .5 * (std::abs(h1.model.y.back().L) + std::abs(h2.model.y.back().L)));
        error = std::max(error, std::abs(full.nuclear_luminosity - half_power) / (options.energy_tolerance * power_scale));
        error = std::max(error, std::abs(full.model.y.back().L / h2.model.y.back().L - 1) / (4 * accuracy_factor * options.structure_tolerance));
        hooks.assess(full.model, full.total_metal_species_rates);
        hooks.assess(h1.model, h1.total_metal_species_rates);
        hooks.assess(h2.model, h2.total_metal_species_rates);
      }
      std::optional<Model> candidate;
      std::vector<std::array<double,3>> rates;
      // Assess only intervals already accepted by the unchanged time error.
      // Failure of this optional improvement must preserve the valid h2 path.
      if (options.richardson_extrapolation && attempt.converged && attempt.audit_pass && attempt.error_norm <= 1) {
        candidate = richardson(full,h1,h2);
        if (candidate && full.total_metal_species_rates.size() != h2.total_metal_species_rates.size()) candidate.reset();
        if (candidate) {
          rates = h2.total_metal_species_rates;
          for (std::size_t i=0;i<rates.size();++i) for(std::size_t k=0;k<3;++k) {
            rates[i][k] += h2.total_metal_species_rates[i][k]-full.total_metal_species_rates[i][k];
            if (!std::isfinite(rates[i][k])) candidate.reset();
          }
        }
        if (candidate) {
          try {
            if (!hooks.assess_extrapolated(state,full,h1,h2,*candidate,rates,ds).empty()) candidate.reset();
          } catch (const std::exception&) { candidate.reset(); }
        }
        if (!candidate) hooks.assess(h2.model,h2.total_metal_species_rates);
      }
      attempt.accepted = attempt.converged && attempt.audit_pass && attempt.error_norm <= 1;
      attempt.message = !full.converged ? full.message : !h1.converged ? h1.message : h2.message;
      if (attempt.converged && !attempt.audit_pass)
        attempt.message = "physical inventory/energy audit failed";
      if (hooks.attempted) hooks.attempted(attempt);
      if (attempt.accepted) {
        previous_power = h2.nuclear_luminosity;
        if (options.predict_structure) {
          preceding = state.model;
          preceding_dt = ds;
        }
        state.model = h2.model;
        state.metal_heat_rates = h2.total_metal_species_rates;
        if (options.richardson_extrapolation) {
          if (candidate) {
            state.model = std::move(*candidate);
            state.metal_heat_rates = std::move(rates);
            ++result.richardson_accepted;
          } else ++result.richardson_declined;
        }
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

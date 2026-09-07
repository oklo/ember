#include "ember/boundary.hpp"
#include "ember/eos_composite.hpp"
#include "ember/relaxation.hpp"
#include "../examples/radiative_polytrope.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0.0) {
  if (!ok) ++failures;
  std::printf("  [%s] %-66s %.6g\n", ok ? "PASS" : "FAIL", label, value);
}
template<class F> static void rejects(F&& f, const char* label) {
  bool threw = false;
  try { f(); } catch (const std::exception&) { threw = true; }
  check(threw, label);
}
static double distance(const Model& a, const Model& b) {
  double error = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
    for (std::size_t v = 0; v < NVAR; ++v) {
      const double scale = v == 3 ? std::max(std::abs(b.y[i].L), 1.0) : 1.0;
      error = std::max(error, std::abs(a.y[i][static_cast<Var>(v)] - b.y[i][static_cast<Var>(v)]) / scale);
    }
  return error;
}

int main() {
  std::printf("ember central boundary and stellar relaxation\n");
  // The leading regular-center error in an n=3 polytrope is -xi^2/15.
  double previous_error = 0.0;
  for (double xi : {0.1, 0.05, 0.025}) {
    example::RadiativePolytrope benchmark(16, xi);
    Physics p{&benchmark.eos, &benchmark.opacity, &benchmark.heating, 1.9};
    const auto c = central_residual(benchmark.reference, p);
    check(std::abs(c.f[0] / (-xi * xi / 15.0) - 1.0) < 0.003,
          "central sphere approaches the regular Lane-Emden expansion", c.f[0]);
    if (previous_error != 0.0)
      check(std::abs(c.f[0] / previous_error - 0.25) < 0.001,
            "halving inner radius quarters the central-boundary error", c.f[0] / previous_error);
    previous_error = c.f[0];
    check(c.f[1] == 0.0, "constant heating satisfies L_inner = epsilon*m_inner");
  }

  // Central derivatives include pp heating, electron degeneracy, and the
  // same gravitational energy convention as the zone equations.
  {
    CompositeEos eos; PPChains nuclear;
    example::PolytropeOpacity opacity; opacity.kappa = 1.0;
    Physics p{&eos, &opacity, &nuclear, 1.9};
    Model m; m.M = 0.1 * constants::Msun; m.m = {1e-8 * m.M, m.M};
    m.comp = {solar_scaled(0.7, 0.014), solar_scaled(0.7, 0.014)};
    m.y = {{std::log(1e7), std::log(1e5), std::log(3e6), 1e25},
           {std::log(1e10), std::log(1.0), std::log(1e5), 1e30}};
    Model old = m; old.y[0].lnT -= 0.02; old.y[0].lnrho -= 0.03;
    for (double dt : {0.0, 1e11}) {
      const auto c = central_residual(m, p, dt, dt > 0.0 ? &old : nullptr);
      double error = 0.0;
      for (std::size_t k = 0; k < 2; ++k) {
        double scale = 0.0, row_error = 0.0;
        for (std::size_t v = 0; v < NVAR; ++v) {
          const auto var = static_cast<Var>(v);
          const double unit = v == 3 ? m.y[0].L : 1.0;
          Model plus = m, minus = m;
          plus.y[0][var] += 1e-5 * unit; minus.y[0][var] -= 1e-5 * unit;
          const double num = (central_residual(plus, p, dt, dt > 0.0 ? &old : nullptr).f[k]
                             - central_residual(minus, p, dt, dt > 0.0 ? &old : nullptr).f[k]) / (2e-5);
          const double analytic = c.dfdy[k][v] * unit;
          row_error = std::max(row_error, std::abs(num - analytic));
          scale = std::max({scale, std::abs(num), std::abs(analytic)});
        }
        error = std::max(error, row_error / scale);
      }
      check(error < 2e-6, "central Jacobian agrees with fixed-previous-state differences", error);
    }
    m.m[0] = 0.0;
    rejects([&] { (void)central_residual(m, p); }, "zero inner mass rejected for logarithmic radius");
  }

  // Recover an independently integrated n=3 polytrope from a perturbed guess.
  // Convergence to the continuum solution must improve with mesh resolution,
  // independently of how small the Newton residual becomes.
  previous_error = 0.0;
  for (std::size_t points : {64UL, 128UL, 256UL}) {
    example::RadiativePolytrope benchmark(points);
    Physics p{&benchmark.eos, &benchmark.opacity, &benchmark.heating, 1.9};
    const auto initial = benchmark.initial();
    const Model untouched = initial;
    const auto r = relax(initial, p, benchmark.atmosphere);
    std::printf("       %zu points: %s; %zu iterations; residual %.3e\n",
                points, r.message.c_str(), r.iterations, r.residual);
    check(r.converged && r.residual < 1e-9 && r.correction < 1e-8,
          "complete model satisfies residual and correction tolerances", r.residual);
    check(distance(initial, untouched) == 0.0, "relaxation leaves its input unchanged");
    if (!r.converged) continue;
    const double error = std::abs(r.model.r(points - 1) / benchmark.reference.r(points - 1) - 1.0);
    check(error < 0.04, "radius approaches the independent Lane-Emden value", error);
    if (previous_error > 0.0)
      check(error / previous_error < 0.3 && error / previous_error > 0.2,
            "doubling the mesh gives second-order radius convergence", error / previous_error);
    previous_error = error;
    if (points == 256) check(error < 0.003, "256-point radius error is below 0.3 percent", error);
    double energy_error = 0.0;
    for (std::size_t i = 0; i < points; ++i)
      energy_error = std::max(energy_error, std::abs(r.model.y[i].L / (benchmark.heating.epsilon * r.model.m[i]) - 1.0));
    check(energy_error < 1e-9, "integrated luminosity equals the enclosed heat generation", energy_error);
    bool decreasing = true;
    for (std::size_t i = 1; i < r.history.size(); ++i)
      decreasing = decreasing && r.history[i].residual < r.history[i - 1].residual;
    check(decreasing, "accepted Newton updates decrease the fixed-scale residual");
    if (points == 64) {
      Model perturbed = r.model;
      for (auto& point : perturbed.y) { point.lnT += 0.01; point.lnrho -= 0.02; }
      const auto dynamic = relax(perturbed, p, benchmark.atmosphere, {}, 1e12, &r.model);
      check(dynamic.converged, "fixed-previous-model thermal solve converges", dynamic.residual);
      check(distance(dynamic.model, r.model) < 1e-7,
            "an equilibrium stays fixed under a thermal step", distance(dynamic.model, r.model));
      check(dynamic.model.age == perturbed.age, "relaxation does not advance age or claim a time integrator");
    }
  }

  // Failure contracts: never return a failed solve as converged or change the
  // caller's model when a step limit or a physics-domain rejection intervenes.
  {
    example::RadiativePolytrope b(32);
    Physics p{&b.eos, &b.opacity, &b.heating, 1.9};
    const Model initial = b.initial();
    RelaxationOptions options; options.max_iterations = 0;
    const auto limited = relax(initial, p, b.atmosphere, options);
    check(!limited.converged && limited.iterations == 0 && distance(limited.model, initial) == 0.0,
          "iteration exhaustion returns failure and the last accepted model");
    options.max_iterations = 1; options.max_log_step = 1e-18;
    const auto stagnant = relax(initial, p, b.atmosphere, options);
    check(!stagnant.converged && stagnant.correction > options.correction_tolerance,
          "an unrepresentably small damped step cannot certify convergence");
    options = {};
    class RejectTrials final : public Atmosphere {
      const Atmosphere& base;
      mutable bool used{};
    public:
      explicit RejectTrials(const Atmosphere& a) : base(a) {}
      AtmosphereState eval(double T, double g, const Composition& c) const override {
        if (used) throw std::domain_error("test atmosphere table edge");
        used = true; return base.eval(T, g, c);
      }
      const char* name() const override { return "trial-rejecting test atmosphere"; }
    } rejecting(b.atmosphere);
    options.max_iterations = 10; options.max_backtracks = 3;
    const auto failed = relax(initial, p, rejecting, options);
    check(!failed.converged && failed.message.find("test atmosphere table edge") != std::string::npos
          && distance(failed.model, initial) == 0.0,
          "line-search failure preserves the model and the physics error");
    Model invalid = initial; invalid.m.back() *= 0.99;
    rejects([&] { (void)relax(invalid, p, b.atmosphere); }, "surface mass must equal the model total mass");
    rejects([&] { (void)relax(initial, p, b.atmosphere, {}, 1.0); }, "thermal solve requires a previous model");
  }
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "ALL PASS", failures);
  return failures ? 1 : 0;
}

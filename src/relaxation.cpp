#include "ember/relaxation.hpp"
#include "ember/boundary.hpp"
#include "ember/henyey.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
namespace {
struct System {
  BoundaryBlock inner, outer;
  std::vector<ZoneResidual> zones;
  double norm{};
};

void validate(const Model& m, const Physics& p, const RelaxationOptions& o, double dt, const Model* prev) {
  if (m.size() < 2 || m.m.size() != m.size() || m.comp.size() != m.size()
      || !p.eos || !p.opacity || !p.nuclear)
    throw std::invalid_argument("relax: invalid model arrays or missing physics");
  if (!std::isfinite(m.M) || !(m.M > 0.0) || m.m.back() != m.M || !(m.m.front() > 0.0))
    throw std::invalid_argument("relax: mesh must start above zero and end at total mass");
  for (std::size_t i = 0; i < m.size(); ++i) {
    if (!std::isfinite(m.m[i]) || (i > 0 && !(m.m[i] > m.m[i - 1])))
      throw std::invalid_argument("relax: mass mesh must be finite and strictly increasing");
    for (std::size_t v = 0; v < NVAR; ++v)
      if (!std::isfinite(m.y[i][static_cast<Var>(v)]))
        throw std::invalid_argument("relax: non-finite initial state");
  }
  if (!(m.y.back().L > 0.0) || !std::isfinite(dt))
    throw std::invalid_argument("relax: positive surface luminosity and finite dt required");
  if (dt > 0.0 && (!prev || prev->m != m.m || prev->size() != m.size() || prev->comp.size() != m.size()))
    throw std::invalid_argument("relax: positive dt requires previous model on the same mesh");
  for (double v : {o.residual_tolerance, o.correction_tolerance, o.max_log_step, o.max_luminosity_step})
    if (!(v > 0.0) || !std::isfinite(v)) throw std::invalid_argument("relax: invalid tolerances or step limits");
  if (o.max_backtracks == 0) throw std::invalid_argument("relax: at least one line-search trial is required");
}

System assemble(const Model& m, const Physics& p, const Atmosphere& atmosphere,
                const std::vector<double>& Lunit, double dt, const Model* prev, bool jacobian) {
  System s{};
  const auto inner = central_residual(m, p, dt, prev);
  const auto outer = surface_residual(m.y.back(), m.M, m.comp.back(), *p.eos, atmosphere);
  s.inner = {inner.f, inner.dfdy}; s.outer = {outer.f, outer.dfdy};
  s.inner.f[1] /= Lunit.front();
  for (double& v : s.inner.dfdy[1]) v /= Lunit.front();
  for (std::size_t k = 0; k < 2; ++k) {
    s.inner.dfdy[k][3] *= Lunit.front(); s.outer.dfdy[k][3] *= Lunit.back();
    s.norm = std::max({s.norm, std::abs(s.inner.f[k]), std::abs(s.outer.f[k])});
  }
  if (jacobian) s.zones.reserve(m.size() - 1);
  for (std::size_t i = 0; i + 1 < m.size(); ++i) {
    ZoneResidual z{};
    if (jacobian) z = zone_residual(m, i, p, dt, prev);
    else z.f = zone_equations(m, i, p, dt, prev);
    const double dm = m.m[i + 1] - m.m[i];
    for (std::size_t k = 0; k < NVAR; ++k) {
      const double scale = k == 2 ? dm / std::max(Lunit[i], Lunit[i + 1]) : dm;
      z.f[k] *= scale;
      if (!std::isfinite(z.f[k])) throw std::domain_error("relax: non-finite scaled zone residual");
      s.norm = std::max(s.norm, std::abs(z.f[k]));
      if (jacobian) {
        for (std::size_t v = 0; v < NVAR; ++v) {
          z.dfdy_lo[k][v] *= scale * (v == 3 ? Lunit[i] : 1.0);
          z.dfdy_hi[k][v] *= scale * (v == 3 ? Lunit[i + 1] : 1.0);
        }
      }
    }
    if (jacobian) s.zones.push_back(z);
  }
  if (!std::isfinite(s.norm)) throw std::domain_error("relax: non-finite scaled residual");
  return s;
}
} // namespace

RelaxationResult relax(const Model& initial, const Physics& p, const Atmosphere& atmosphere,
                       const RelaxationOptions& options, double dt, const Model* prev) {
  validate(initial, p, options, dt, prev);
  // Fix units for the whole solve, including every line-search trial. Small
  // central luminosities get a local enclosed-mass scale, without a solar
  // luminosity floor or division by luminosity that can cross zero.
  double reference_L = 0.0;
  for (const auto& point : initial.y) reference_L = std::max(reference_L, std::abs(point.L));
  std::vector<double> Lunit(initial.size());
  for (std::size_t i = 0; i < initial.size(); ++i) {
    Lunit[i] = std::max(std::abs(initial.y[i].L), reference_L * (initial.m[i] / initial.M));
    if (!(Lunit[i] > 0.0) || !std::isfinite(Lunit[i]))
      throw std::domain_error("relax: luminosity units are not representable");
  }
  RelaxationResult result{};
  result.model = initial;
  for (;;) {
    const auto system = assemble(result.model, p, atmosphere, Lunit, dt, prev, true);
    result.residual = system.norm;
    result.correction = std::numeric_limits<double>::infinity();
    HenyeyCorrection correction;
    try { correction = solve_henyey(system.inner, system.zones, system.outer); }
    catch (const std::runtime_error& e) {
      result.message = e.what();
      return result;
    }
    double log_step = 0.0, L_step = 0.0;
    for (const auto& dy : correction.dy) {
      for (std::size_t v = 0; v < 3; ++v) log_step = std::max(log_step, std::abs(dy[v]));
      L_step = std::max(L_step, std::abs(dy[3]));
    }
    result.correction = std::max(log_step, L_step);
    if (result.residual <= options.residual_tolerance && result.correction <= options.correction_tolerance) {
      result.converged = true; result.message = "converged";
      result.history.push_back({result.residual, result.correction, 0.0, correction.backward_error});
      return result;
    }
    if (result.iterations >= options.max_iterations) {
      result.message = "Newton iteration limit reached";
      return result;
    }
    double damping = 1.0;
    if (log_step > options.max_log_step) damping = std::min(damping, options.max_log_step / log_step);
    if (L_step > options.max_luminosity_step) damping = std::min(damping, options.max_luminosity_step / L_step);
    const double surface_dL = correction.dy.back()[3] * Lunit.back();
    if (surface_dL < 0.0)
      damping = std::min(damping, -0.8 * result.model.y.back().L / surface_dL);
    bool accepted = false;
    std::string last_rejection = "residual did not decrease";
    for (std::size_t trial = 0; trial < options.max_backtracks; ++trial) {
      Model candidate = result.model;
      for (std::size_t i = 0; i < candidate.size(); ++i)
        for (std::size_t v = 0; v < NVAR; ++v)
          candidate.y[i][static_cast<Var>(v)] += damping * correction.dy[i][v] * (v == 3 ? Lunit[i] : 1.0);
      try {
        const double norm = assemble(candidate, p, atmosphere, Lunit, dt, prev, false).norm;
        if (norm < system.norm && norm <= (1.0 - 1e-4 * damping) * system.norm) {
          result.model = std::move(candidate); accepted = true; break;
        }
        last_rejection = "residual did not decrease";
      } catch (const std::domain_error& e) {
        last_rejection = e.what();
      } catch (const std::out_of_range& e) {
        last_rejection = e.what();
      } catch (const std::runtime_error& e) {
        last_rejection = e.what();
      }
      damping *= 0.5;
    }
    if (!accepted) {
      result.message = "line search failed: " + last_rejection;
      return result;
    }
    result.history.push_back({result.residual, result.correction, damping, correction.backward_error});
    ++result.iterations;
  }
}

} // namespace ember

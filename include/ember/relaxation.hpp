#pragma once
#include "ember/atmosphere.hpp"
#include "ember/structure.hpp"
#include <string>
#include <vector>
#include <limits>

namespace ember {

struct RelaxationOptions {
  std::size_t max_iterations{50};
  std::size_t max_backtracks{24};
  double residual_tolerance{1e-9};
  double correction_tolerance{1e-8};
  double max_log_step{0.2};
  double max_luminosity_step{0.5};  // fraction of each point's fixed luminosity unit
};

struct RelaxationIteration {
  double residual{}, correction{}, damping{}, linear_error{};
};

struct RelaxationResult {
  Model model;  // last accepted state, including on failure; input is never modified
  bool converged{};
  std::size_t iterations{};  // accepted Newton updates
  double residual{}, correction{std::numeric_limits<double>::infinity()};
  std::string message;
  std::vector<RelaxationIteration> history;
};

// Damped Newton solve on a fixed mass mesh and fixed composition. dt<=0 is
// static; positive dt uses a fixed previous model on exactly the same mesh.
// Does not advance age, burn composition, remesh, or choose a time step.
// Invalid inputs/initial physics throw. Numerical nonconvergence is returned
// explicitly, with diagnostics and the last accepted model.
RelaxationResult relax(const Model& initial, const Physics&, const Atmosphere&,
                       const RelaxationOptions& = {}, double dt = 0.0,
                       const Model* prev = nullptr);

} // namespace ember

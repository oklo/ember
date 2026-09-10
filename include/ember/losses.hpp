#pragma once
#include "ember/composition.hpp"
#include <cmath>
#include <stdexcept>

namespace ember {

// Specific thermal-neutrino energy sinks. Their relevance must be assessed along
// the actual trajectory, including remnant contraction and cooling; initial
// stellar mass alone does not establish that every channel is negligible.
struct LossState {
  double eps{};          // erg/g/s, positive is a sink
  double dlneps_dlnT{};
  double dlneps_dlnRho{};
};

class NeutrinoLosses {
public:
  virtual ~NeutrinoLosses() = default;
  virtual LossState eval(double T, double rho, const Composition&) const = 0;
  virtual const char* name() const = 0;
};

// Explicit zero-loss implementation for configured controls.
class NoNeutrinoLosses final : public NeutrinoLosses {
public:
  LossState eval(double, double, const Composition&) const override { return {}; }
  const char* name() const override { return "none"; }
};

// Haft, Raffelt & Weiss (1994), equations 23--27. Plasma decay only, using
// the fully ionized electron count and the source weak mixing angle .23.
// The published few-percent accuracy concerns the total emission rate where
// plasma decay dominates; it is not an all-domain error bound on this channel.
// Pair/photo/bremsstrahlung/recombination and partial ionization are omitted.
class PlasmaNeutrinoLosses final : public NeutrinoLosses {
public:
  LossState eval(double T,double rho,const Composition&) const override;
  const char* name() const override { return "Haft-Raffelt-Weiss plasma neutrinos; fully ionized approximation only"; }
};

inline LossState evaluate_losses(const NeutrinoLosses* source,double T,double rho,const Composition& comp) {
  const auto result=source?source->eval(T,rho,comp):LossState{};
  if(!std::isfinite(result.eps) || result.eps<0 || !std::isfinite(result.dlneps_dlnT)
      || !std::isfinite(result.dlneps_dlnRho))
    throw std::domain_error("NeutrinoLosses: invalid sink or derivative");
  return result;
}

} // namespace ember

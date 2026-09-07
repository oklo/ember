#pragma once
#include "ember/eos.hpp"
#include "ember/eos_component.hpp"

namespace ember {

// An equation of state assembled from additive components.
//
// The default assembly - ions, radiation, electrons - is the analytic one.
// Extending to cold, dense matter means pushing components on:
//
//     eos.add<CoulombOCP>();        // ion lattice energy, Gamma up to melting
//     eos.add<Crystallisation>();   // latent heat and the Debye solid
//     eos.add<PhaseSeparation>();   // C/O separation on freezing
//
// and nothing else in the code changes.  That is the point of the shape: the
// white dwarfs this code will be asked about next - a solar remnant, or
// something near the Chandrasekhar mass - differ from the present target in
// exactly these terms and in no other way.
class CompositeEos final : public Eos {
public:
  CompositeEos();                       // ions + radiation + electrons
  explicit CompositeEos(std::vector<std::unique_ptr<EosComponent>> parts)
      : parts_(std::move(parts)) {}

  template <class C, class... Args>
  CompositeEos& add(Args&&... args) {
    parts_.push_back(std::make_unique<C>(std::forward<Args>(args)...));
    return *this;
  }

  EosState eval(double T, double rho, const Composition&) const override;
  EosResponse eval_with_derivatives(double T, double rho, const Composition&) const override;
  const char* name() const override { return name_.c_str(); }
  std::size_t size() const { return parts_.size(); }

private:
  std::vector<std::unique_ptr<EosComponent>> parts_;
  std::string name_{"composite"};
};

} // namespace ember

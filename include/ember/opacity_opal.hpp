#pragma once
#include "ember/opacity_table.hpp"

namespace ember {
// Fixed-Z table; errors at every edge. See data/opacity/README.md for coverage.
class OpalOpacity final : public TabulatedOpacity {
public:
  explicit OpalOpacity(const std::filesystem::path& file)
      : TabulatedOpacity(file, "OPAL (GS98, rectangular subset)") {}
};
} // namespace ember

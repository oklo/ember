#pragma once
#include "ember/opacity_table.hpp"

namespace ember {
// AESOPUS 2.1 (Marigo et al. 2024), pressure-broadened gas opacity.
// The shipped GS98 table has logT=2..4.5 and logR=-8..6 at Z=.020.
// No grain opacity: the source archive is explicitly the gas distribution.
class AesopusOpacity final : public TabulatedOpacity {
public:
  explicit AesopusOpacity(const std::filesystem::path& file)
      : TabulatedOpacity(file, "AESOPUS 2.1 gas (GS98)") {}
};
} // namespace ember

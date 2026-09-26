// Apply the same exact source-retention check used by stellar restarts.
#include "ember/eos_mixture.hpp"
#include <iostream>
#include <string_view>

int main(int argc, char** argv) {
  try {
    if (argc != 4)
      throw std::invalid_argument("check_eos_extension density|temperature original.dat candidate.dat");
    const std::string_view kind(argv[1]);
    if (kind != "density" && kind != "temperature")
      throw std::invalid_argument("extension must be density or temperature");
    const auto mixture = ember::HelmholtzTableEos::Mixture::allow_documented_proxy;
    ember::MetalHelmholtzEos original(argv[2], mixture), candidate(argv[3], mixture);
    const auto added = kind == "density"
        ? candidate.check_density_extension(original)
        : candidate.check_temperature_extension(original);
    std::cout << "{\"retained_source_values_and_masks_exact\":true,\"added_states\":"
              << added << "}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

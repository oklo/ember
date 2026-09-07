#include "ember/relaxation.hpp"
#include "../examples/radiative_polytrope.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <exception>

int main(int argc, char** argv) {
  std::size_t points = 128;
  if (argc > 2) { std::fprintf(stderr, "usage: ember-polytrope [mesh_points]\n"); return 2; }
  if (argc == 2) {
    const char* end = argv[1] + std::strlen(argv[1]);
    const auto parsed = std::from_chars(argv[1], end, points);
    if (parsed.ec != std::errc{} || parsed.ptr != end) {
      std::fprintf(stderr, "mesh_points must be an integer\n"); return 2;
    }
  }
  try {
    ember::example::RadiativePolytrope benchmark(points);
    ember::Physics physics{&benchmark.eos, &benchmark.opacity, &benchmark.heating, 1.9};
    const auto result = ember::relax(benchmark.initial(), physics, benchmark.atmosphere);
    if (!result.converged) {
      std::fprintf(stderr, "polytrope solve failed: %s (residual %.6g, correction %.6g)\n",
                   result.message.c_str(), result.residual, result.correction);
      return 1;
    }
    const auto& m = result.model;
    std::printf("{\n  \"benchmark\": \"n=3 radiative polytrope; controlled physics, artificial atmosphere\",\n");
    std::printf("  \"converged\": true, \"points\": %zu, \"iterations\": %zu,\n", points, result.iterations);
    std::printf("  \"residual\": %.16g, \"correction\": %.16g,\n", result.residual, result.correction);
    std::printf("  \"mass_g\": %.16g, \"radius_cm\": %.16g, \"luminosity_erg_s\": %.16g,\n",
                m.M, m.r(points - 1), m.y.back().L);
    std::printf("  \"radius_relative_error\": %.16g,\n", m.r(points - 1) / benchmark.reference.r(points - 1) - 1.0);
    std::printf("  \"history\": [\n");
    for (std::size_t i = 0; i < result.history.size(); ++i) {
      const auto& h = result.history[i];
      std::printf("    {\"residual\": %.16g, \"correction\": %.16g, \"damping\": %.16g, \"linear_error\": %.16g}%s\n",
                  h.residual, h.correction, h.damping, h.linear_error, i + 1 == result.history.size() ? "" : ",");
    }
    std::printf("  ],\n  \"columns\": [\"mass_g\", \"radius_cm\", \"density_g_cm3\", \"temperature_K\", \"luminosity_erg_s\"],\n  \"profile\": [\n");
    for (std::size_t i = 0; i < m.size(); ++i)
      std::printf("    [%.16g, %.16g, %.16g, %.16g, %.16g]%s\n", m.m[i], m.r(i), m.rho(i), m.T(i),
                  m.y[i].L, i + 1 == m.size() ? "" : ",");
    std::printf("  ]\n}\n");
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what()); return 1;
  }
}

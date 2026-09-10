#pragma once
#include <array>
namespace ember {
// Generated from the pinned GS98 TOPS request and synple-elements.json.
// Metal fractions and representative integer isotopes match the atmosphere
// generator. Atomic weights support source-mass conversion; this is not a
// complete stellar isotope inventory.
struct Gs98Metal { double charge, mass_number, atomic_weight, fraction; };
inline constexpr std::array<Gs98Metal, 19> gs98_metals{{
  {6, 12, 12.0107, 0.17183599999999999}, // C
  {7, 14, 14.006740000000001, 0.050334999999999998}, // N
  {8, 16, 15.9994, 0.46735599999999999}, // O
  {10, 20, 20.1797, 0.10483099999999998}, // Ne
  {11, 23, 22.98977, 0.0020899999999999998}, // Na
  {12, 24, 24.305, 0.039924000000000001}, // Mg
  {13, 27, 26.981539999999999, 0.0036029999999999999}, // Al
  {14, 28, 28.0855, 0.044056999999999999}, // Si
  {15, 31, 30.973759999999999, 0.00042299999999999998}, // P
  {16, 32, 32.066000000000003, 0.023513000000000003}, // S
  {17, 35, 35.4527, 0.000292}, // Cl
  {18, 40, 39.948, 0.0043349999999999994}, // Ar
  {19, 39, 39.098300000000002, 0.00022800000000000007}, // K
  {20, 40, 40.078000000000003, 0.0038959999999999993}, // Ca
  {22, 48, 47.866999999999997, 0.00019500000000000002}, // Ti
  {24, 52, 51.996099999999998, 0.0011169999999999999}, // Cr
  {25, 55, 54.938049999999997, 0.00077899999999999996}, // Mn
  {26, 56, 55.844999999999999, 0.076433000000000001}, // Fe
  {28, 59, 58.693399999999997, 0.0047569999999999991}, // Ni
}};
constexpr double gs98_ion_moment(int power) {
  double sum=0;
  for(const auto& e:gs98_metals) {
    double weight=e.fraction/e.mass_number;
    for(int i=0;i<power;++i)weight*=e.charge;
    sum+=weight;
  }
  return sum;
}
constexpr double gs98_atomic_mass_scale() {
  double sum=0;for(const auto& e:gs98_metals)sum+=e.fraction*e.atomic_weight/e.mass_number;
  return sum;
}
} // namespace ember

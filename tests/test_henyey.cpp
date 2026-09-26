#include "ember/henyey.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace ember;
static int failures = 0;
static void check(bool ok, const char* label, double value = 0.0) {
  if (!ok) ++failures;
  std::printf("  [%s] %-60s %.6g\n", ok ? "PASS" : "FAIL", label, value);
}

// Independent dense Gauss-Jordan reference, with no use of the Henyey blocks
// after the complete matrix is assembled. Used only for small systems.
static std::vector<double> dense_solve(std::vector<std::vector<double>> a) {
  const std::size_t n = a.size();
  for (auto& row : a) {
    double norm = 0.0;
    for (std::size_t j = 0; j < n; ++j) norm = std::max(norm, std::abs(row[j]));
    for (double& v : row) v /= norm;
  }
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
    std::swap(a[k], a[pivot]);
    const double value = a[k][k];
    if (value == 0.0) throw std::runtime_error("singular test reference");
    for (double& v : a[k]) v /= value;
    for (std::size_t i = 0; i < n; ++i) {
      if (i == k) continue;
      const double factor = a[i][k];
      for (std::size_t j = 0; j <= n; ++j) a[i][j] -= factor * a[k][j];
    }
  }
  std::vector<double> result(n);
  for (std::size_t i = 0; i < n; ++i) result[i] = a[i][n];
  return result;
}

// This fixed, nonsingular three-point system reproduced backward error one
// in the original solver solely from roundoff in a zero boundary value.
static void held_boundary_regression() {
  BoundaryBlock inner{},outer{};std::vector<ZoneResidual> zones(2);
  inner.f={-1.2779163432084399,-0.11940608868942114};
  inner.dfdy[0]={1,0.33333333333333331,0,0};
  inner.dfdy[1]={0,0,0.13,1};
  outer.f={0,0};
  outer.dfdy[0]={0,0,1,0};
  outer.dfdy[1]={0,0.71999999999999997,0.28000000000000003,0};
  zones[0].f={2.3563395610015059,1.123428513923683,-0.55017329508732082,-0.30517315159161673};
  zones[0].dfdy_lo[0]={-1.5382176593881898,0.54522231372417906,0.10856665149336919,0.39159182649142332};
  zones[0].dfdy_lo[1]={-0.64248210125975025,-1.1739488988255629,0.26390841455676828,0.055575145480202322};
  zones[0].dfdy_lo[2]={0.18031010816552162,-0.095430368616358383,-1.1862380077168004,0.037062601488677724};
  zones[0].dfdy_lo[3]={0.24039261934761999,-0.45076236651786239,0.061641369458058377,-0.33160628300267619};
  zones[0].dfdy_hi[0]={1.390428510554351,0.62971031212648099,-0.61020673306863915,0.19671440544479879};
  zones[0].dfdy_hi[1]={-0.53914449678784193,0.57014376080613205,-0.66178249789920685,0.23841534569003939};
  zones[0].dfdy_hi[2]={0.44493649770537425,0.20769490853089526,0.58422398656583219,0.56303014137405472};
  zones[0].dfdy_hi[3]={0.03516366230628476,0.43351904472773489,0.12204579924888487,1.5033235131791358};
  zones[1].f={-2.7844092500016369,-2.1900305773201452,-0.6254838950552899,0.90312051393906112};
  zones[1].dfdy_lo[0]={-1.3618429931813192,-0.66750003552760484,0.25988269098284705,0.64687021554581192};
  zones[1].dfdy_lo[1]={-0.57248983913682061,-1.5096309101748906,-0.29026003695756358,0.048772291386540378};
  zones[1].dfdy_lo[2]={0.30247600152946891,0.53234268713103938,-1.6821802345891921,0.13578025734554192};
  zones[1].dfdy_lo[3]={-0.63674843203938658,0.094441340285491351,0.30809149784635781,-0.72544862492907536};
  zones[1].dfdy_hi[0]={0.92650986567488025,0.30917900205441001,-0.015873718656732105,0.61645288540519594};
  zones[1].dfdy_hi[1]={0.23009527410308361,1.4747410051274241,0.36188043578117018,-0.53606428060598632};
  zones[1].dfdy_hi[2]={-0.031427452115107522,0.15012246888368275,0.7821875087755642,-0.61216184955044994};
  zones[1].dfdy_hi[3]={-0.63062033872723444,0.50845898250040522,-0.54393411806913394,1.4600128820926681};
  const std::vector<std::array<double,4>> expected{{0.99998672587241755,0.83378885200806718,-0.88754262715998233,0.23478663022021884},{-0.97388115740795855,-0.74497031964948113,-0.69174889007088658,0.67505065332405856},{0.85570261060684305,0,0,-0.14412326835876843}};
  for(bool reverse:{false,true}) {
    if(reverse) {std::swap(outer.f[0],outer.f[1]);std::swap(outer.dfdy[0],outer.dfdy[1]);}
    const auto result=solve_henyey(inner,zones,outer);double error=0;
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<4;++j)
      error=std::max(error,std::abs(result.dy[i][j]-expected[i][j]));
    check(error<1e-12,"held-boundary regression recovers the known solution",error);
    check(result.dy.back()[1]==0 && result.dy.back()[2]==0,"both held boundary values are exactly zero");
    check(result.backward_error<1e-12,"held-boundary correction satisfies every original equation",result.backward_error);
  }
}

int main() {
  held_boundary_regression();
  std::printf("ember Henyey block elimination\n");
  for(bool held_boundary : {false,true})for (std::size_t points : {1UL, 2UL, 7UL, 31UL, 2048UL}) {
    BoundaryBlock inner{}, outer{};
    // These inner rows have zero leading columns, forcing elimination to
    // pivot with the zone equations instead of assuming an invertible 2x2.
    inner.dfdy[0][held_boundary?0:2] = 1.0; inner.dfdy[1][3] = 1.0;
    outer.dfdy[0][held_boundary?2:0] = 1.0; outer.dfdy[1][1] = 1.0;
    if(held_boundary)outer.dfdy[1][2]=.2;
    std::vector<ZoneResidual> zones(points - 1);
    const double h = 1.0 / static_cast<double>(points);
    for (std::size_t i = 0; i < zones.size(); ++i) {
      auto& z = zones[i];
      for (std::size_t k = 0; k < NVAR; ++k) {
        const double unit = std::pow(10.0, 80.0 * std::sin(static_cast<double>(4 * i + k)));
        for (std::size_t v = 0; v < NVAR; ++v) {
          const double angle = static_cast<double>(7 * i + 3 * k + v);
          z.dfdy_lo[k][v] = unit * (-(k == v ? 1.0 : 0.0) + 0.1 * h * std::sin(angle));
          z.dfdy_hi[k][v] = unit * ((k == v ? 1.0 : 0.0) + 0.2 * h * std::cos(angle));
        }
      }
    }
    std::vector<std::array<double, NVAR>> expected(points);
    for (std::size_t i = 0; i < points; ++i)
      for (std::size_t v = 0; v < NVAR; ++v)
        expected[i][v] = std::sin(0.31 * static_cast<double>(4 * i + v)) + 0.2;
    // Exactly fixed outer temperature and pressure. Relative backward error
    // must not reject an arbitrarily small spurious correction to a variable
    // whose defining equation determines an exact zero.
    if(held_boundary)expected.back()[1]=expected.back()[2]=0;
    for (std::size_t k = 0; k < 2; ++k)
      for (std::size_t v = 0; v < NVAR; ++v) {
        inner.f[k] -= inner.dfdy[k][v] * expected.front()[v];
        outer.f[k] -= outer.dfdy[k][v] * expected.back()[v];
      }
    for (std::size_t i = 0; i < zones.size(); ++i)
      for (std::size_t k = 0; k < NVAR; ++k)
        for (std::size_t v = 0; v < NVAR; ++v)
          zones[i].f[k] -= zones[i].dfdy_lo[k][v] * expected[i][v]
                           + zones[i].dfdy_hi[k][v] * expected[i + 1][v];
    HenyeyCorrection solved;
    try {solved=solve_henyey(inner,zones,outer);}
    catch(const std::runtime_error& e) {
      std::printf("       %zu points, held boundary %d: %s\n",points,held_boundary,e.what());
      check(false,"linear system with a known solution must converge");continue;
    }
    double error = 0.0;
    for (std::size_t i = 0; i < points; ++i)
      for (std::size_t v = 0; v < NVAR; ++v)
        error = std::max(error, std::abs(solved.dy[i][v] - expected[i][v]));
    std::printf("       %zu mesh points, held boundary %d\n", points,held_boundary);
    if(held_boundary)check(solved.dy.back()[1]==0 && solved.dy.back()[2]==0,
      "zero outer corrections obey their defining boundary equations exactly");
    check(error < 1e-10, "known solution recovered despite row scales spanning 160 decades", error);
    check(solved.backward_error < 1e-12, "correction satisfies the original linear equations", solved.backward_error);
    if (points <= 31) {
      const std::size_t n = NVAR * points;
      std::vector<std::vector<double>> a(n, std::vector<double>(n + 1));
      for (std::size_t k = 0; k < 2; ++k) {
        for (std::size_t v = 0; v < NVAR; ++v) {
          a[k][v] = inner.dfdy[k][v]; a[n - 2 + k][n - NVAR + v] = outer.dfdy[k][v];
        }
        a[k][n] = -inner.f[k]; a[n - 2 + k][n] = -outer.f[k];
      }
      for (std::size_t i = 0; i < zones.size(); ++i)
        for (std::size_t k = 0; k < NVAR; ++k) {
          const auto row = 2 + NVAR * i + k;
          for (std::size_t v = 0; v < NVAR; ++v) {
            a[row][NVAR * i + v] = zones[i].dfdy_lo[k][v];
            a[row][NVAR * (i + 1) + v] = zones[i].dfdy_hi[k][v];
          }
          a[row][n] = -zones[i].f[k];
        }
      const auto dense = dense_solve(a);
      error = 0.0;
      for (std::size_t i = 0; i < points; ++i)
        for (std::size_t v = 0; v < NVAR; ++v)
          error = std::max(error, std::abs(solved.dy[i][v] - dense[NVAR * i + v]));
      check(error < 1e-11, "block result agrees with independent dense elimination", error);
    }
  }
  bool rejected = false;
  try { (void)solve_henyey({}, {}, {}); } catch (const std::runtime_error&) { rejected = true; }
  check(rejected, "singular systems are rejected");
  BoundaryBlock bad{}; bad.dfdy[0][0] = std::numeric_limits<double>::quiet_NaN();
  rejected = false;
  try { (void)solve_henyey(bad, {}, {}); } catch (const std::domain_error&) { rejected = true; }
  check(rejected, "non-finite coefficients are rejected");
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "ALL PASS", failures);
  return failures ? 1 : 0;
}

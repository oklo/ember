#pragma once
#include "ember/detail/differential.hpp"
#include "ember/detail/taylor3.hpp"
namespace ember::detail {
inline double scalar(double x){return x;}
template<class D> double scalar(const D& x) {
  if constexpr(requires {x.value();})return x.value();
  else return x.value;
}
// Potekhin & Chabrier (2000) one-component-plasma Coulomb free energy per ion
// (CD09 eq. 24). A small-Gamma series keeps the Debye--Hueckel limit exact.
template<class D> D classical_ocp_free_energy(const D& g) {
  using std::sqrt;using std::log;using std::log1p;using std::atan;
  constexpr double A1=-.907,A2=.62954,B1=.00456,B2=211.6,B3=-1e-4,B4=.00462;
  const double A3=-std::sqrt(3.)/2-A1/std::sqrt(A2);
  const D s=sqrt(g),sa=sqrt(g/A2);
  if(scalar(g)<1e-3) {
    // Series of the same fit avoid subtracting nearly equal square roots,
    // logarithms and arctangents in the dilute limit.
    const D x=g/A2,z=g/B2;
    return A1*A2*sa*x*(2./3+x*(-1./5+x*(3./28+x*(-5./72+x*35./704))))
      +2*A3*s*g*(1./3+g*(-1./5+g*(1./7+g*(-1./9+g/11))))
      +B1*g*z*(.5+z*(-1./3+z*(.25-z/5)))+.5*B3*log1p(g*g/B4);
  }
  return A1*(sqrt(g*(A2+g))-A2*log(sa+sqrt(1+g/A2)))+2*A3*(s-atan(s))
        +B1*(g-B2*log1p(g/B2))+.5*B3*log1p(g*g/B4);
}
} // namespace ember::detail

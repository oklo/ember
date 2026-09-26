#pragma once
#include "ember/microscopic_transport.hpp"
#include "../src/differential.hpp"
#include <cmath>
#include <stdexcept>

// Prescribed analytic coefficients for integration tests, not a stellar
// diffusion prescription. Candidate-state dependence and all geometric,
// abundance and thermal cross terms are deliberately nonzero.
struct AnalyticMicroscopicTransport final:ember::MicroscopicTransport {
  double diffusion{1.},heat_conductivity{1e11},thermal_force{.01},enthalpy{3e14};
  ember::MicroscopicFaceResponse eval(std::size_t,double mlo,double mhi,
      const ember::Point& lo,const ember::Composition& a,
      const ember::Point& hi,const ember::Composition& b,bool derivatives)const override {
    using namespace ember;using D=detail::Differential<8>;using detail::exp;using detail::log;
    const D r=.5*(exp(D::variable(lo.lnr,0))+exp(D::variable(hi.lnr,4)));
    const D rho=.5*(exp(D::variable(lo.lnrho,1))+exp(D::variable(hi.lnrho,5)));
    const D T=.5*(exp(D::variable(lo.lnT,2))+exp(D::variable(hi.lnT,6)));
    const D logratio=D::variable(hi.lnT,6)-D::variable(lo.lnT,2);
    const D area=4*M_PI*r*r;
    const D g=area*area*rho*rho/(mhi-mlo)*diffusion*exp(.3*log(T/1e7)-.2*log(rho/1e3));
    // This simple composition-dependent mobility gives nonzero d(mobility)/dX.
    const double factor=1+.2*(a.X[0]+b.X[0]);
    std::array<D,2> rates;
    rates[0]=g*factor*(a.X[0]-b.X[0]-thermal_force*logratio);
    rates[1]=.7*g*factor*(a.X[1]-b.X[1]+.2*thermal_force*logratio);
    const D heat=heat_conductivity*exp(.7*log(T/1e7)+.4*log(rho/1e3));
    const D carried=enthalpy*T/1e7*(rates[0]-.4*rates[1]);
    MicroscopicFaceResponse out;out.carried_luminosity=carried.value;out.conductivity=heat.value;
    for(std::size_t j=0;j<2;++j) {
      out.species.rate[j]=rates[j].value;
      if(derivatives) {
        out.species.dleft[j][j]=g.value*factor*(j?.7:1.);
        out.species.dright[j][j]=-out.species.dleft[j][j];
        out.species.dleft[j][0]+=.2*rates[j].value/factor;
        out.species.dright[j][0]+=.2*rates[j].value/factor;
      }
    }
    if(derivatives)for(std::size_t v=0;v<NVAR;++v) {
      out.dcarried_lo[v]=carried.d[v];out.dcarried_hi[v]=carried.d[v+4];
      out.dconductivity_lo[v]=heat.d[v];out.dconductivity_hi[v]=heat.d[v+4];
    }
    return out;
  }
  const char* name()const override{return "analytic microscopic integration test";}
};

// Exercise a thermal response defined before nuclear production of He3,
// while species evaluation requires a positive new-time He3 candidate.
struct ThermalSubspaceControl final:ember::MicroscopicTransport {
  AnalyticMicroscopicTransport value;
  bool total_species_heat{false};
  mutable std::size_t absent_heat_calls{},species_calls{};
  mutable std::size_t total_heat_calls{};
  ember::MicroscopicFaceResponse eval(std::size_t face,double mlo,double mhi,
      const ember::Point& lo,const ember::Composition& a,
      const ember::Point& hi,const ember::Composition& b,bool derivatives)const override {
    ++species_calls;
    if(a.X[1]==0 || b.X[1]==0)throw std::domain_error("species derivatives require newly produced He3");
    return value.eval(face,mlo,mhi,lo,a,hi,b,derivatives);
  }
  ember::MicroscopicHeatResponse heat(std::size_t face,double mlo,double mhi,
      const ember::Point& lo,const ember::Composition& a,
      const ember::Point& hi,const ember::Composition& b,bool derivatives)const override {
    if(a.X[1]==0 && b.X[1]==0)++absent_heat_calls;
    return value.eval(face,mlo,mhi,lo,a,hi,b,derivatives);
  }
  bool requires_positive_species_guess()const override{return true;}
  bool uses_total_species_heat()const override{return total_species_heat;}
  ember::MicroscopicHeatResponse heat_with_total_species_rate(std::size_t face,double mlo,double mhi,
      const ember::Point& lo,const ember::Composition& a,const ember::Point& hi,
      const ember::Composition& b,const ember::SpeciesVector& total,bool derivatives)const override {
    ++total_heat_calls;
    auto r=heat(face,mlo,mhi,lo,a,hi,b,derivatives);
    // Synthetic split: half the original coefficient is material enthalpy,
    // half kinetic transported heat. Only the former follows total mixing.
    const double Tlo=std::exp(lo.lnT),Thi=std::exp(hi.lnT),T=.5*(Tlo+Thi);
    const double h=.5*value.enthalpy*T/1e7;
    r.total_rate_enthalpy={h,-.4*h};
    const double extra=h*(total[0]-.4*total[1]);
    r.carried_luminosity=.5*r.carried_luminosity+extra;
    for(std::size_t i=0;i<ember::NVAR;++i) {
      r.dcarried_lo[i]*=.5;r.dcarried_hi[i]*=.5;
    }
    if(derivatives){r.dcarried_lo[2]+=extra*.5*Tlo/T;r.dcarried_hi[2]+=extra*.5*Thi/T;}
    return r;
  }
  const char* name()const override{return "separate thermal-domain integration control";}
};

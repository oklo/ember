#pragma once
#include "ember/stellar_seed.hpp"
#include "ember/boundary.hpp"
#include "ember/structure.hpp"
#include "ember/relaxation.hpp"

namespace ember::driver {
// Entropy loss constructs only the age-zero star. Actual evolution uses
// gravothermal energy differences and the physical nuclear source.
struct ContractingSource final:Nuclear {
 const Nuclear& source;double entropy_loss;
 ContractingSource(const Nuclear& n,double rate):source(n),entropy_loss(rate){}
 NuclearState eval(double T,double rho,const Composition& c)const override {
  auto s=source.eval(T,rho,c);const double cooling=T*entropy_loss,total=s.eps+cooling;
  s.dlneps_dlnT=(s.eps*s.dlneps_dlnT+cooling)/total;
  s.dlneps_dlnRho*=s.eps/total;s.eps=total;return s;
 }
 const char* name()const override{return "initial contracting model only";}
};
inline Model contracting_guess(std::size_t points,double mass,double radius,double teff,
    const Composition& c,const Physics& physics,const Atmosphere& atmosphere,
    double envelope_mass=0) {
 const auto& eos=*physics.eos;const auto& seed_source=*physics.nuclear;
 auto guess=example::stellar_seed(points,mass,radius,c,seed_source,atmosphere,1.5,teff,
                                  LuminosityGrid::volume_faces,envelope_mass);
    // Initial guess only: integrate the selected discrete stellar equations
    // inward through the envelope, then blend into the polytropic interior.
    // Global relaxation and the physical acceptance criteria remain unchanged.
    const Model polytrope = guess;
    const auto surface = atmosphere.eval(teff, constants::G*guess.M/
        std::pow(guess.r(guess.size()-1),2), c);
    guess.y.back().lnT=std::log(surface.T);
    guess.y.back().lnrho=std::log(surface.rho);
    const double initial_L = guess.y.back().L;
    auto norm4=[](const std::array<double,4>& a) {
      double n=0.;for(auto v:a)n=std::max(n,std::abs(v));return n;
    };
    for(std::size_t hi=guess.size()-1;hi>0 && guess.m[hi-1]/guess.M>.90;--hi) {
      const std::size_t lo=hi-1;const double dm=guess.m[hi]-guess.m[lo];
      const auto e=eos.eval(guess.T(hi),guess.rho(hi),c);
      const double Plo=e.P+constants::G*.5*(guess.m[hi]+guess.m[lo])*dm/
        (4*M_PI*std::pow(guess.r(hi),4));
      const double Tlo=guess.T(hi)*std::pow(Plo/e.P,e.grad_ad);
      guess.y[lo]={guess.y[hi].lnr-dm/(4*M_PI*std::pow(guess.r(hi),3)*guess.rho(hi)),
        std::log(eos.rho_from_PT(Tlo,Plo,c,guess.rho(hi))),std::log(Tlo),
        guess.y[hi].L-energy_interval_mass(guess,lo)*seed_source.eval(guess.T(hi),guess.rho(hi),c).eps};
      bool converged=false;
      for(std::size_t it=0;it<40;++it) {
        const auto z=zone_residual(guess,lo,physics,0.);
        std::array<double,4> f{};std::array<std::array<double,5>,4> mat{};
        for(std::size_t k=0;k<4;++k) {
          const double scale=k==2?dm/initial_L:dm;f[k]=z.f[k]*scale;
          for(std::size_t v=0;v<4;++v)
            mat[k][v]=z.dfdy_lo[k][v]*scale*(v==3?initial_L:1.);
          mat[k][4]=-f[k];
        }
        const double before=norm4(f);
        if(before<1e-10) {
          converged=true;break;
        }
        for(std::size_t k=0;k<4;++k) {
          std::size_t pivot=k;
          for(std::size_t j=k+1;j<4;++j)if(std::abs(mat[j][k])>std::abs(mat[pivot][k]))pivot=j;
          if(std::abs(mat[pivot][k])<1e-15)throw std::runtime_error("singular envelope initializer");
          std::swap(mat[pivot],mat[k]);
          const double scale=mat[k][k];for(std::size_t j=k;j<5;++j)mat[k][j]/=scale;
          for(std::size_t i=0;i<4;++i)if(i!=k) {
            const double factor=mat[i][k];for(std::size_t j=k;j<5;++j)mat[i][j]-=factor*mat[k][j];
          }
        }
        double damping=1.;for(std::size_t v=0;v<4;++v)
          if(std::abs(mat[v][4])>.2)damping=std::min(damping,.2/std::abs(mat[v][4]));
        const auto base=guess.y[lo];bool accepted=false;
        for(std::size_t trial=0;trial<24;++trial) {
          for(std::size_t v=0;v<4;++v)guess.y[lo][static_cast<Var>(v)]=
            base[static_cast<Var>(v)]+damping*mat[v][4]*(v==3?initial_L:1.);
          try {
            auto candidate=zone_equations(guess,lo,physics,0.);
            for(std::size_t k=0;k<4;++k)candidate[k]*=k==2?dm/initial_L:dm;
            if(norm4(candidate)<before){accepted=true;break;}
          }catch(const std::exception&){}
          damping*=.5;
        }
        if(!accepted)throw std::runtime_error("envelope initializer line search failed at "+std::to_string(lo));
      }
      if(!converged)throw std::runtime_error("envelope initializer iteration limit at "+std::to_string(lo));
    }
    for(std::size_t i=0;i<guess.size();++i) {
      double f=std::clamp((guess.m[i]/guess.M-.90)/.08,0.,1.);f=f*f*(3-2*f);
      for(std::size_t v=0;v<4;++v)guess.y[i][static_cast<Var>(v)]=
        f*guess.y[i][static_cast<Var>(v)]+(1-f)*polytrope.y[i][static_cast<Var>(v)];
    }

 return guess;
}
} // namespace ember::driver

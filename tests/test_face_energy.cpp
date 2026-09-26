#include "ember/energy_grid.hpp"
#include "ember/evolution.hpp"
#include "ember/boundary.hpp"
#include "ember/constants.hpp"
#include "ember/evolution_checkpoint.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>

using namespace ember;
namespace {
int checks=0;
void require(bool ok,const char* message){++checks;if(!ok)throw std::runtime_error(message);}
template<class F> bool rejects(F action){try{action();}catch(const std::exception&){return true;}return false;}
constexpr double R=constants::R_gas/.6,chemical_energy=1e14;
struct Gas final:Eos {
  EosState eval(double T,double rho,const Composition& c)const override {
    EosState s;s.P=R*rho*T;s.E=1.5*R*T+chemical_energy*c.X[0];
    s.cv=1.5*R;s.cp=2.5*R;s.delta=s.chiT=s.chiRho=1;s.grad_ad=.4;return s;
  }
  EosResponse eval_with_derivatives(double T,double rho,const Composition& c)const override {
    EosResponse s;s.state=eval(T,rho,c);return s;
  }
  const char* name()const override{return "ideal gas with constant composition energy";}
};
struct Radiation final:Opacity {
  OpacityState eval(double,double,const Composition&)const override{return {1,0,0};}
  const char* name()const override{return "constant radiative opacity";}
};
struct NoBurn final:Nuclear {
  NuclearState eval(double,double,const Composition&)const override{return {};}
  const char* name()const override{return "no reactions";}
};
struct ConstantHeatFlux final:MicroscopicTransport {
  std::vector<double> K;
  MicroscopicFaceResponse eval(std::size_t i,double,double,const Point&,const Composition&,
      const Point&,const Composition&,bool)const override {
    MicroscopicFaceResponse s;s.carried_luminosity=K.at(i);return s;
  }
  const char* name()const override{return "prescribed constant-enthalpy material exchange";}
};
Model fixture() {
  Model m;m.M=5e31;m.m={1e31,2e31,3e31,4e31,5e31};
  auto c=solar_scaled(.5,.02);c.basis=AbundanceBasis::baryon_mass;
  for(std::size_t i=0;i<m.m.size();++i) {
    const double rho=1e3*std::exp(-.2*static_cast<double>(i));
    m.y.push_back({std::log(2e9*(1+.1*static_cast<double>(i))),std::log(rho),std::log(1e6),1e30});
    m.comp.push_back(c);
  }
  m.y[0].lnr=std::log(3*m.m[0]/(4*M_PI*m.rho(0)))/3;
  return m;
}
double derivative_check(const Model& m,const Model& old,const Physics& p,double dt) {
  double worst=0;
  for(std::size_t i=0;i+1<m.size();++i) {
    const auto a=zone_residual(m,i,p,dt,&old),b=zone_residual_numerical(m,i,p,dt,&old,2e-5);
    const auto values=zone_equations(m,i,p,dt,&old);
    require(a.f==values,"analytic and values-only residuals differ");
    for(std::size_t row=0;row<NVAR;++row) {
      double scale=0,error=0;
      for(std::size_t side=0;side<2;++side)for(std::size_t col=0;col<NVAR;++col) {
        const double unit=col==3?1e32:1.;
        const double x=(side?a.dfdy_hi:a.dfdy_lo)[row][col]*unit;
        const double y=(side?b.dfdy_hi:b.dfdy_lo)[row][col]*unit;
        scale=std::max({scale,std::abs(x),std::abs(y)});error=std::max(error,std::abs(x-y));
      }
      worst=std::max(worst,error/std::max(scale,1e-100));
    }
  }
  const auto a=central_residual(m,p,dt,&old);
  for(std::size_t row=0;row<2;++row) {
    double error=0,scale=0;
    for(std::size_t col=0;col<NVAR;++col) {
      const double unit=col==3?1e32:1.,h=2e-5*unit;
      auto plus=m,minus=m;plus.y[0][static_cast<Var>(col)]+=h;minus.y[0][static_cast<Var>(col)]-=h;
      const double expected=(central_residual(plus,p,dt,&old).f[row]-central_residual(minus,p,dt,&old).f[row])/(2*h)*unit;
      const double actual=a.dfdy[row][col]*unit;
      error=std::max(error,std::abs(actual-expected));scale=std::max({scale,std::abs(actual),std::abs(expected)});
    }
    worst=std::max(worst,error/scale);
  }
  require(worst<2e-6,"energy/transport Jacobian differs from finite differences");return worst;
}
}
int main() {
  try {
    Gas eos;Radiation opacity;NoBurn nuclear;Physics p{&eos,&opacity,&nuclear,1.9};
    auto old=fixture(),m=old;m.luminosity_grid=old.luminosity_grid=LuminosityGrid::volume_faces;
    ConstantHeatFlux exchange;exchange.K={0,1e32,0,0};p.microscopic=&exchange;
    constexpr double dt=1e11;
    const auto w=nodal_mass_weights(m);
    for(std::size_t i=0;i<m.size();++i) {
      const double left=i?exchange.K[i-1]:0,right=i+1<m.size()?exchange.K[i]:0;
      const double change=-dt*(right-left)/(chemical_energy*w[i]);
      m.comp[i].X[0]+=change;m.comp[i].X[2]-=change;m.y[i].L=right;
    }
    double maximum_energy=std::abs(central_residual(m,p,dt,&old).f[1])/1e32;
    double maximum_thermal=0;
    for(std::size_t i=0;i+1<m.size();++i) {
      const auto z=zone_equations(m,i,p,dt,&old);
      maximum_energy=std::max(maximum_energy,std::abs(z[2])*w[i+1]/1e32);
      maximum_thermal=std::max(maximum_thermal,std::abs(z[3])*(m.m[i+1]-m.m[i]));
    }
    require(maximum_energy<2e-12,"constant-enthalpy exchange fails local energy balance");
    require(maximum_thermal<1e-14,"material heat generates spurious thermal flux");
    require(convective_mixing_regions(m,p).size()==m.size(),"pure material heat falsely triggers convection");
    // Independently integrate the older nodal equation for exactly the same
    // divergence. Its interpolation leaves half the heat flux on this face.
    auto nodal=m;nodal.luminosity_grid=LuminosityGrid::mass_nodes;
    std::vector<double> S(m.size());
    for(std::size_t i=0;i<m.size();++i)
      S[i]=((i+1<m.size()?exchange.K[i]:0)-(i?exchange.K[i-1]:0))/w[i];
    nodal.y[0].L=nodal.m[0]*S[0];
    for(std::size_t i=1;i<m.size();++i)
      nodal.y[i].L=nodal.y[i-1].L+.5*(S[i-1]+S[i])*(nodal.m[i]-nodal.m[i-1]);
    const double nodal_artifact=std::abs(.5*(nodal.y[1].L+nodal.y[2].L)-exchange.K[1])/exchange.K[1];
    require(std::abs(nodal_artifact-.5)<1e-14,"manufactured old-grid counterexample changed");

    // Uneven spacing, changing thermal state and neutrino losses exercise
    // every energy volume, including both boundaries, rather than only a
    // uniform mesh where an incorrect interval width could happen to agree.
    PlasmaNeutrinoLosses losses;p.neutrino_losses=&losses;p.microscopic=nullptr;
    old=fixture();old.m={1e29,3e30,3.2e30,4.9e31,5e31};old.luminosity_grid=LuminosityGrid::volume_faces;
    m=old;
    for(std::size_t i=0;i<m.size();++i) {
      m.y[i].lnT+=.03*(1+static_cast<double>(i));m.y[i].lnrho+=.02;
      m.y[i].L=1e31*(1+static_cast<double>(i));
    }
    const auto weights=nodal_mass_weights(m);
    double residual=central_residual(m,p,dt,&old).f[1];
    long double direct=m.y.back().L;
    for(std::size_t i=0;i<m.size();++i) {
      require(nodal_volume_mass(m,i)==weights[i],"energy volume differs from composition volume");
      const auto e=eos.eval(m.T(i),m.rho(i),m.comp[i]),e0=eos.eval(old.T(i),old.rho(i),old.comp[i]);
      const long double heat=-((static_cast<long double>(e.E)-e0.E)+
          static_cast<long double>(e.P)*(1.L/m.rho(i)-1.L/old.rho(i)))/dt;
      direct-=weights[i]*(heat-losses.eval(m.T(i),m.rho(i),m.comp[i]).eps);
      if(i)residual+=weights[i]*zone_equations(m,i-1,p,dt,&old)[2];
    }
    const double identity_error=std::abs(static_cast<double>((residual-direct)/direct));
    require(identity_error<2e-12,"whole-star energy does not telescope on uneven mass mesh");
    const double jacobian_error=derivative_check(m,old,p,dt);
    // A zero material-heat provider cannot change the temperature equation,
    // convection boundaries or mixing coefficients of a volume-face model.
    // A nonzero temperature contrast exposes the old arithmetic/log mean
    // mismatch; an isothermal test would miss it.
    auto zero=p;ConstantHeatFlux absent;absent.K.assign(m.size()-1,0);
    p.microscopic=nullptr;zero.microscopic=&absent;
    for(std::size_t i=0;i+1<m.size();++i) {
      require(zone_equations(m,i,p,dt,&old)==zone_equations(m,i,zero,dt,&old),
          "zero material heat changes volume-face structure equations");
    }
    require(convective_mixing_regions(m,p)==convective_mixing_regions(m,zero),
        "zero material heat changes convection boundaries");
    const auto mix0=convective_mixing_faces(m,p),mix1=convective_mixing_faces(m,zero);
    for(std::size_t i=0;i<mix0.size();++i)
      require(mix0[i].diffusivity==mix1[i].diffusivity,
          "zero material heat changes convective diffusivity");
    (void)derivative_check(m,old,p,dt);
    auto incompatible=old;incompatible.luminosity_grid=LuminosityGrid::mass_nodes;
    require(rejects([&]{zone_residual(m,0,p,dt,&incompatible);}),"changed luminosity placement accepted by energy row");
    require(rejects([&]{central_residual(m,p,dt,&incompatible);}),"changed luminosity placement accepted at center");

    // Checkpoint carries the convention and isotope inventory, even after D
    // has disappeared. Old callers must explicitly opt into the new format.
    m=fixture();m.luminosity_grid=LuminosityGrid::volume_faces;
    for(auto& c:m.comp) {
      c.metal_inventory=MetalInventory::gs98;c.cn_molality=initial_gs98_cn(c);c=explicit_cn_material(c);
    }
    const auto stamp=std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto path=std::filesystem::temp_directory_path()/("ember-face-"+std::to_string(stamp)+".checkpoint");
    const driver::Selections selection{"test","test","test","test","test"};
    const driver::Identities identities{{"executable","numerical-test"}};
    for(double D:{2e-5,std::numeric_limits<double>::denorm_min(),0.}) {
      for(auto& c:m.comp){c[Species::H2]=D;c.X[0]=.5-D;}
      driver::write_checkpoint(path,{m,dt,7,2},selection,1e-12,identities);
      const auto saved=driver::read_checkpoint(path,m.size(),m.M,m.comp[0],selection,1e-12,identities,LuminosityGrid::volume_faces);
      require(saved.model.luminosity_grid==m.luminosity_grid && saved.model.comp==m.comp,"face checkpoint loses grid or D/CN inventory");
      for(std::size_t i=0;i<m.size();++i)require(saved.model.y[i].L==m.y[i].L,"checkpoint alters luminosity");
      require(rejects([&]{driver::read_checkpoint(path,m.size(),m.M,m.comp[0],selection,1e-12,identities);}),"nodal-only reader accepts face checkpoint");
      std::filesystem::remove(path);
    }
    std::cout<<std::setprecision(17)<<"{\"outcome\":\"passed\",\"checks\":"<<checks
      <<",\"exchange_energy_error\":"<<maximum_energy<<",\"spurious_thermal_gradient\":"<<maximum_thermal
      <<",\"old_nodal_flux_error\":"<<nodal_artifact<<",\"uneven_mesh_energy_identity\":"<<identity_error
      <<",\"jacobian_error\":"<<jacobian_error<<"}\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<"after "<<checks<<" checks: "<<e.what()<<'\n';return 1;}
}

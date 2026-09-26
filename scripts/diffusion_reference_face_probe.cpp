// Read-only comparison of native transport and its ideal-mixture force limit.
// This does not change or evolve a checkpoint.
#include "ember/screened_microscopic_transport.hpp"
#include "ember/eos_component.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace ember;

int main(int argc,char** argv) {
  if(argc!=3)return 2;
  try {
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);
    ScreenedMicroscopicTransport transport(eos,collisions,true,2e6);
    ElectronGas electrons;
    std::string line;
    while(std::getline(std::cin,line)) {
      std::istringstream input(line);
      std::array<double,2> m,r,rho,T,X,Y;
      for(int i=0;i<2;++i)if(!(input>>m[i]>>r[i]>>rho[i]>>T[i]>>X[i]>>Y[i]))return 2;
      try {
        std::array<Composition,2> c;
        std::array<Point,2> point;
        std::array<CompositionPotentialResponse,2> potential;
        for(int i=0;i<2;++i) {
          c[i]=solar_scaled(X[i],.02);c[i].basis=AbundanceBasis::baryon_mass;
          c[i].metal_inventory=MetalInventory::gs98;c[i].X[1]=Y[i];c[i].X[2]=.98-X[i]-Y[i];
          point[i]={std::log(r[i]),std::log(rho[i]),std::log(T[i]),1.};
          potential[i]=eos.active_composition_potential(T[i],rho[i],c[i],{true,true});
        }
        const double rb=.5*(r[0]+r[1]),rhob=.5*(rho[0]+rho[1]),Tb=std::sqrt(T[0]*T[1]);
        auto cb=c[0];for(std::size_t j=0;j<NSPEC;++j)cb.X[j]=.5*(c[0].X[j]+c[1].X[j]);
        const double ne=cb.mu_elec_inv()*constants::NA*rhob;
        const auto electron=electrons.eval(Tb,rhob,cb);
        const double stiffness=electron.dP_dlnRho/ne;
        const double length=ScreenedCollisionTransport::screening_length(Tb,rhob,cb.X[0],cb.X[1],.02,stiffness,true);
        const auto kinetic=collisions.eval(Tb,rhob,cb.X[0],cb.X[1],.02,length);
        const auto heat=eos.composition_heat(Tb,rhob,cb,{true,true},false);
        const auto native=transport.eval(0,m[0],m[1],point[0],c[0],point[1],c[1],false);
        const double area=4*M_PI*rb*rb,dr=(m[1]-m[0])/(area*rhob),dtemp=(T[1]-T[0])/(T[0]*T[1]);
        const double dlogrho=std::log(rho[1]/rho[0]),dlogT=std::log(T[1]/T[0]);
        const double dlogne=dlogrho+std::log(c[1].mu_elec_inv()/c[0].mu_elec_inv());
        const double dlogHe4=dlogrho+std::log(c[1].X[2]/c[0].X[2]);
        std::array<double,2> ideal_phi{
          constants::R_gas*(dlogrho+std::log(X[1]/X[0])-.25*dlogHe4+.5*dlogne-1.875*dlogT),
          constants::R_gas*((dlogrho+std::log(Y[1]/Y[0]))/3-.25*dlogHe4+dlogne/6-.375*dlogT)};
        std::array<double,2> ideal_h{3.125*constants::R_gas*Tb,.625*constants::R_gas*Tb};
        std::array<double,2> ideal_force{},native_force{},ideal_rate{},recomputed{};
        for(int j=0;j<2;++j) {
          ideal_force[j]=ideal_phi[j]+(ideal_h[j]+kinetic.transport_enthalpy[j])*dtemp;
          native_force[j]=potential[1].gradient[j]-potential[0].gradient[j]
              +(heat.exchange_enthalpy[j]+kinetic.transport_enthalpy[j])*dtemp;
        }
        for(int i=0;i<2;++i)for(int j=0;j<2;++j) {
          ideal_rate[i]-=area/dr*kinetic.mobility[i][j]*ideal_force[j];
          recomputed[i]-=area/dr*kinetic.mobility[i][j]*native_force[j];
        }
        std::cout<<std::setprecision(17)<<"{\"T\":"<<Tb<<",\"rho\":"<<rhob
          <<",\"X\":"<<cb.X[0]<<",\"Y3\":"<<cb.X[1]<<",\"radius\":"<<rb
          <<",\"dr_mass\":"<<dr<<",\"dr_geometric\":"<<r[1]-r[0]
          <<",\"eta\":"<<kinetic.eta<<",\"electron_stiffness_over_kT\":"<<stiffness/(constants::R_gas/constants::NA*Tb)
          <<",\"screening_length\":"<<length<<",\"native_rate\":["<<native.species.rate[0]<<','<<native.species.rate[1]
          <<"],\"recomputed_rate\":["<<recomputed[0]<<','<<recomputed[1]
          <<"],\"ideal_force_rate\":["<<ideal_rate[0]<<','<<ideal_rate[1]
          <<"],\"native_H_velocity\":"<<native.species.rate[0]/(area*rhob*cb.X[0])
          <<",\"ideal_force_H_velocity\":"<<ideal_rate[0]/(area*rhob*cb.X[0])
          <<",\"native_force\":["<<native_force[0]<<','<<native_force[1]
          <<"],\"ideal_force\":["<<ideal_force[0]<<','<<ideal_force[1]<<"]}\n";
      }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n";}
      std::cout.flush();
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

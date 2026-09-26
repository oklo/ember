#include "ember/collision_transport.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

int main(int argc,char** argv) {
  if(argc<2 || argc>3) return 2;
  const bool derivatives=argc==3 && std::string(argv[2])=="derivatives";
  if(argc==3 && !derivatives) return 2;
  try {
    ember::ScreenedCollisionTransport physics(argv[1]);
    std::string line;
    while(std::getline(std::cin,line)) {
      double T,rho,X,Y3,Z,length;std::istringstream input(line);
      if(!(input>>T>>rho>>X>>Y3>>Z>>length)) return 2;
      try {
        double stiffness=0;int include_ions=0;
        const bool screened=static_cast<bool>(input>>stiffness>>include_ions);
        if(screened)
          length=ember::ScreenedCollisionTransport::screening_length(T,rho,X,Y3,Z,stiffness,include_ions!=0);
        const auto start=std::chrono::steady_clock::now();
        ember::CollisionTransportDerivatives d;
        const auto r=derivatives ? (d=physics.derivatives(T,rho,X,Y3,Z,length),d.value) : physics.eval(T,rho,X,Y3,Z,length);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<std::setprecision(17)<<"{\"screening_length\":"<<length<<",\"eta\":"<<r.eta<<",\"ne\":"<<r.electron_density
          <<",\"b\":"<<r.b_thermal<<",\"energy_scale\":"<<r.energy_scale<<",\"conductivity\":"<<r.conductivity
          <<",\"enthalpy\":["<<r.transport_enthalpy[0]<<','<<r.transport_enthalpy[1]<<"],\"mobility\":[";
        for(std::size_t i=0;i<3;++i) {if(i)std::cout<<',';std::cout<<'[';
          for(std::size_t j=0;j<3;++j) {if(j)std::cout<<',';std::cout<<r.mobility[i][j];}std::cout<<']';}
        std::cout<<"],\"active\":["<<r.active_species[0]<<','<<r.active_species[1]<<"],\"backward_error\":"
          <<r.maximum_solve_backward_error<<",\"seconds\":"<<seconds;
        if(derivatives) {
          std::cout<<",\"defined\":[";
          for(std::size_t k=0;k<6;++k) {if(k)std::cout<<',';std::cout<<d.defined[k];}
          std::cout<<"],\"partials\":[";
          for(std::size_t k=0;k<6;++k) {
            if(k)std::cout<<',';const auto& a=d.partials[k];
            std::cout<<"{\"eta\":"<<a.eta<<",\"ne\":"<<a.electron_density<<",\"b\":"<<a.b_thermal
              <<",\"energy_scale\":"<<a.energy_scale<<",\"conductivity\":"<<a.conductivity<<",\"enthalpy\":["
              <<a.transport_enthalpy[0]<<','<<a.transport_enthalpy[1]<<"],\"mobility\":[";
            for(std::size_t i=0;i<3;++i) {if(i)std::cout<<',';std::cout<<'[';
              for(std::size_t j=0;j<3;++j) {if(j)std::cout<<',';std::cout<<a.mobility[i][j];}std::cout<<']';}
            std::cout<<"]}";
          }
          std::cout<<']';
          if(screened) {
            const auto s=ember::ScreenedCollisionTransport::screening_derivatives(T,rho,X,Y3,Z,stiffness,include_ions!=0);
            std::cout<<",\"screening_partials\":[";
            for(std::size_t k=0;k<6;++k) {if(k)std::cout<<',';std::cout<<s.partials[k];}
            std::cout<<']';
          }
        }
        std::cout<<"}\n";
      } catch(const std::exception& e) {std::cout<<"{\"error\":\""<<e.what()<<"\"}\n";}
      std::cout.flush();
    }
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

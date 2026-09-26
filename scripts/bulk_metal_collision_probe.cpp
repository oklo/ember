#include "ember/collision_transport.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <type_traits>

namespace {
template<class T,std::size_t N> void array(const std::array<T,N>& a) {
  std::cout<<'[';
  for(std::size_t i=0;i<N;++i) {
    if(i)std::cout<<',';
    if constexpr(std::is_arithmetic_v<T>)std::cout<<a[i];else array(a[i]);
  }
  std::cout<<']';
}
template<class R> void fields(const R& r) {
  std::cout<<"\"eta\":"<<r.eta<<",\"ne\":"<<r.electron_density<<",\"b\":"<<r.b_thermal
      <<",\"energy_scale\":"<<r.energy_scale<<",\"conductivity\":"<<r.conductivity
      <<",\"enthalpy\":";array(r.transport_enthalpy);
  std::cout<<",\"mobility\":";array(r.mobility);
}
template<std::size_t N> void response(const ember::BasicCollisionTransportDerivatives<N>& r) {
  std::cout<<'{';fields(r.value);std::cout<<",\"active\":";array(r.value.active_species);
  std::cout<<",\"backward_error\":"<<r.value.maximum_solve_backward_error<<",\"defined\":";array(r.defined);
  std::cout<<",\"partials\":[";
  for(std::size_t k=0;k<6;++k) {if(k)std::cout<<',';std::cout<<'{';fields(r.partials[k]);std::cout<<'}';}
  std::cout<<"]}";
}
}
int main(int argc,char** argv) {
  if(argc!=2)return 2;
  try {
    ember::ScreenedCollisionTransport collision(argv[1]);
    std::cout<<std::setprecision(17);std::string line;
    while(std::getline(std::cin,line)) {
      double T,rho,X,Y3,Z,length;std::istringstream input(line);
      if(!(input>>T>>rho>>X>>Y3>>Z>>length))return 2;
      try {
        const auto bulk=collision.bulk_metal_derivatives(T,rho,X,Y3,Z,length);
        const auto fixed=collision.derivatives(T,rho,X,Y3,Z,length);
        std::cout<<"{\"bulk\":";response(bulk);std::cout<<",\"fixed\":";response(fixed);std::cout<<"}\n";
      }catch(const std::exception& e){std::cout<<"{\"error\":\""<<e.what()<<"\"}\n";}
      std::cout.flush();
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

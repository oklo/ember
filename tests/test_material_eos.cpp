#include "ember/eos_material.hpp"
#include "ember/constants.hpp"
#include "ember/detail/taylor3.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
using namespace ember;
namespace {
// An independent, composition-dependent caloric potential. No table or phase
// solver is involved, so this checks the common structure/transport interface.
class AnalyticMaterial final:public MaterialEos {
 public:
  std::array<HelmholtzJet,10> material_jets(double T,double rho,const Composition&c,std::size_t n)const override{
    validate_composition_domain(T,rho,c);
    using J=detail::Taylor3<5>;
    const auto t=J::variable(std::log(T),0),r=J::variable(std::log(rho),1);
    const auto x=J::variable(c.X[0],2),y=J::variable(c.X[1],3),z=J::variable(c.Z(),4);
    const auto dt=t-std::log(1e5);
    const auto f=constants::R_gas*((1+x+.1*z+.03*x*y)*r-(1.5+.4*x+.5*y)*t+.02*z*dt*dt+2*x*z);
    constexpr std::array<std::array<unsigned,3>,10> powers{{{0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}}};
    std::array<HelmholtzJet,10> result{};
    for(unsigned k=0;k<n;++k){auto p=powers[k];unsigned order=p[0]+p[1]+p[2];
      for(unsigned i=0;i+order<=3;++i)for(unsigned j=0;i+j+order<=3;++j)
        result[k][i][j]=f.derivative({i,j,p[0],p[1],p[2]});}
    return result;
  }
  void validate_composition_domain(double T,double rho,const Composition&c)const override{
    if(!(T>0&&rho>0&&c.X[2]>0))throw std::domain_error("invalid analytic material state");
  }
  const char* name()const override{return "analytic material regression";}
};
double relative(double a,double b){return std::abs(a-b)/std::max({1.,std::abs(a),std::abs(b)});}
Composition shifted(Composition c,unsigned k,double h){c[k==2?Species::Zrest:static_cast<Species>(k)]+=h;c.X[2]-=h;return c;}
}
int main(){
 AnalyticMaterial eos;double force_error=0,heat_error=0,response_error=0;unsigned failures=0;
 Composition c;c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
 c.X[0]=.2;c.X[1]=.001;c.X[2]=.759;c[Species::Zrest]=.04;
 for(double T:{1e4,1e6})for(double rho:{1e-3,1e3}){
  const auto state=eos.eval(T,rho,c);const auto p=eos.composition_potential(T,rho,c);
  const auto h=eos.composition_heat(T,rho,c);const auto iso=eos.isobaric_composition_response(T,rho,c);
  for(unsigned k=0;k<3;++k){
   const double dx=1e-4*(k==2?c.Z():c.X[k]);auto a=shifted(c,k,dx),b=shifted(c,k,-dx);
   const auto pa=eos.composition_potential(T,rho,a),pb=eos.composition_potential(T,rho,b);
   force_error=std::max(force_error,relative((pa.phi-pb.phi)/(2*dx),p.gradient[k]));
   for(unsigned j=0;j<3;++j)force_error=std::max(force_error,relative((pa.gradient[j]-pb.gradient[j])/(2*dx),p.hessian[j][k]));
   const double ra=eos.rho_from_PT(T,state.P,a,rho),rb=eos.rho_from_PT(T,state.P,b,rho);
   const auto ea=eos.eval(T,ra,a),eb=eos.eval(T,rb,b);
   heat_error=std::max(heat_error,relative(((ea.E+ea.P/ra)-(eb.E+eb.P/rb))/(2*dx),h.exchange_enthalpy[k]+h.radiation_enthalpy[k]));
   response_error=std::max(response_error,relative(std::log(ra/rb)/(2*dx),iso.dlnRho[k]));
  }
  for(unsigned k=0;k<2;++k){
   const double d=1e-5,ta=k==0?T*std::exp(d):T,tb=k==0?T*std::exp(-d):T;
   const double ra=k==1?rho*std::exp(d):rho,rb=k==1?rho*std::exp(-d):rho;
   const auto a=eos.composition_heat(ta,ra,c),b=eos.composition_heat(tb,rb,c);
   for(unsigned j=0;j<3;++j){
    heat_error=std::max(heat_error,relative((a.exchange_enthalpy[j]-b.exchange_enthalpy[j])/(2*d),h.enthalpy_partials[j][k]));
    heat_error=std::max(heat_error,relative((a.radiation_enthalpy[j]-b.radiation_enthalpy[j])/(2*d),h.radiation_enthalpy_partials[j][k]));
   }
  }
 }
 if(force_error>1e-5 || heat_error>1e-5 || response_error>1e-5)++failures;
 std::cout<<"material EOS: force "<<force_error<<", heat "<<heat_error<<", density response "<<response_error<<", failures "<<failures<<'\n';
 return failures?1:0;
}

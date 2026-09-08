// Offline audit utility. stdin: ln(rho), ln(T); argv[1]: potential table.
// stdout: status, P,E,S,chiRho,chiT,cv,cp,grad_ad,delta,dE_dlnrho,
// dcp_dlnT,dcp_dlnrho,ddelta_dlnT,ddelta_dlnrho,dgrad_ad_dlnT,dgrad_ad_dlnrho.
#include "ember/eos_helmholtz.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
int main(int argc,char** argv) {
 using namespace ember;
 HelmholtzTableEos eos(argc>1?argv[1]:"data/eos/freeeos300_hhe_x070_potential.dat",HelmholtzTableEos::Mixture::allow_documented_proxy);
 const auto c=solar_scaled(.7,.02);
 double lr,lt;
 while(std::cin>>lr>>lt) {
  try {
   const auto a=eos.eval_with_derivatives(std::exp(lt),std::exp(lr),c);const auto&e=a.state;
   std::cout<<std::setprecision(17)<<0<<' '<<e.P<<' '<<e.E<<' '<<e.S<<' '<<e.chiRho<<' '<<e.chiT<<' '<<e.cv<<' '<<e.cp<<' '<<e.grad_ad<<' '<<e.delta<<' '<<a.dE_dlnRho<<' '<<a.dcp_dlnT<<' '<<a.dcp_dlnRho<<' '<<a.ddelta_dlnT<<' '<<a.ddelta_dlnRho<<' '<<a.dgrad_ad_dlnT<<' '<<a.dgrad_ad_dlnRho<<'\n';
  }catch(const std::exception&e) {std::cout<<"1 "<<e.what()<<'\n';}
 }
}

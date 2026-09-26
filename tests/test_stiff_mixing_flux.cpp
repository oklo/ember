#include "../src/flux_chain.hpp"
#include <iostream>
#include <iomanip>
using namespace ember::detail;
int main() {
  try {
    constexpr std::size_t n=7;
    std::size_t checks=0;double worst=0;
    auto check=[&](bool ok){++checks;if(!ok)throw std::runtime_error("stiff diffusion state or conservative flux differs from analytic solution");};
    for(double g:{0.,1e-8,1.,1e6,1e12,1e18,1e30}) {
      std::vector<FluxMatrix<n>> E(2),A(1),B(1);
      std::vector<FluxVector<n>> rhs(2),offset(1),flux;
      FluxVector<n> xleft{},xright{};
      for(std::size_t k=0;k<n;++k) {
        const double scale=k==6?1e-300:std::pow(.1,static_cast<double>(k));
        xleft[k]=.8*scale;xright[k]=.2*scale;
        E[0][k][k]=.25;E[1][k][k]=.75;A[0][k][k]=B[0][k][k]=g;
        rhs[0][k]=.25*xleft[k];rhs[1][k]=.75*xright[k];
      }
      const auto answer=solve_flux_chain<n>(E,A,B,rhs,offset,&flux);
      for(std::size_t k=0;k<n;++k) {
        const double scale=xleft[k],mean=.25*xleft[k]+.75*xright[k];
        // Exact two-cell backward-Euler exchange, formed without subtracting
        // the two almost equal final abundances.
        const double f=g==0?0:(xleft[k]-xright[k])/(1/g+1/.25+1/.75);
        const double a=xleft[k]-f/.25,b=xright[k]+f/.75;
        const double error=std::max({std::abs(answer[0][k]-a),std::abs(answer[1][k]-b),std::abs(flux[0][k]-f)})/scale;
        worst=std::max(worst,error);check(error<2e-14);
        if(g==0)check(flux[0][k]==0);
        else check(std::abs(flux[0][k]/f-1)<2e-13);
        check(std::abs((.25*answer[0][k]+.75*answer[1][k]-mean)/scale)<2e-14);
        check(std::abs((.25*(answer[0][k]-xleft[k])+flux[0][k])/scale)<2e-14);
        check(std::abs((.75*(answer[1][k]-xright[k])-flux[0][k])/scale)<2e-14);
      }
      if(g>=1e18)check(flux[0][0]>.1 && std::abs(answer[0][0]-answer[1][0])<1e-15);
    }
    std::cout<<std::setprecision(17)<<checks<<" analytic stiff-mixing checks; maximum scaled error "<<worst<<'\n';
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}

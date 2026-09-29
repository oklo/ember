#pragma once
#include "ember/eos_helmholtz.hpp"
#include "ember/constants.hpp"
#include <cmath>
namespace ember::detail {
// Linear mixing at a common electron density. Differentiate both the species
// weights and the electron density as H1, He3 or GS98 metals replace He4.
// The evaluator returns F/T and thermal derivatives for one pure ion species.
template<class Evaluate> std::array<HelmholtzJet,10> common_density_ion_jets(
    double T,double rho,const Composition& c,std::size_t channels,Evaluate evaluate) {
  const double Z=c.Z(),Ye=c.X[0]+(2./3)*c.X[1]+.5*c.X[2]+gs98_ion_moment(1)*Z;
  const double logT=std::log(T),logne=std::log(rho*Ye/constants::amu);
  std::array<HelmholtzJet,4> species{};
  species[0]=evaluate(logT,logne,1,1);
  species[1]=evaluate(logT,logne,3,2);
  species[2]=evaluate(logT,logne,4,2);
  if(Z>0 || channels>1)for(const auto& m:gs98_metals) {
    const auto f=evaluate(logT,logne,m.mass_number,m.charge);
    for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)species[3][i][j]+=m.fraction*f[i][j];
  }
  const std::array<double,4> fractions{c.X[0],c.X[1],c.X[2],Z};
  const std::array<double,3> beta{.5/Ye,(1./6)/Ye,(gs98_ion_moment(1)-.5)/Ye};
  std::array<HelmholtzJet,10> out{};
  for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)
    for(unsigned s=0;s<4;++s)out[0][i][j]+=fractions[s]*species[s][i][j];
  constexpr std::array<unsigned,3> which{0,1,3};
  if(channels>1)for(unsigned a=0;a<3;++a)
    for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j<=2;++j)
      out[1+a][i][j]=species[which[a]][i][j]-species[2][i][j]+beta[a]*out[0][i][j+1];
  if(channels>4) {
    unsigned ch=4;
    for(unsigned a=0;a<3;++a)for(unsigned b=a;b<3;++b,++ch)
      for(unsigned i=0;i<2;++i)for(unsigned j=0;i+j<=1;++j)
        out[ch][i][j]=beta[a]*(species[which[b]][i][j+1]-species[2][i][j+1])
          +beta[b]*(species[which[a]][i][j+1]-species[2][i][j+1])
          +beta[a]*beta[b]*(out[0][i][j+2]-out[0][i][j+1]);
  }
  return out;
}
} // namespace ember::detail

#include "ember/ion_quantum.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
using namespace ember;
int main() {
 try {
  int failed=0;
  auto check=[&](bool ok,const char* what,double value=0.) {
    std::printf("%s %s %.5g\n",ok?"PASS":"FAIL",what,value);failed+=!ok;
  };
  auto mixture=[](double X,double Y3,double Z) {
    auto c=solar_scaled(X,Z);c.basis=AbundanceBasis::baryon_mass;
    c.metal_inventory=MetalInventory::gs98;c.X[1]=Y3;c.X[2]-=Y3;return c;
  };
  // Independent values from the author's LIQUBC (eos22.f, June 2022),
  // SHA256 f514e1781f4f476c8e23ad03848488445aa539dd7ccbdda049dd9ed91d7f8ae9.
  // Compare F,U,P,Cv; derivatives are checked against the free energy below,
  // rather than relying on LIQUBC's known pressure-derivative transcription.
  struct Reference{double X,Y3,T,rho,A,f,u,p,cv;};
  const Reference refs[]{
    {0,0,1e6,3e4,4,0.019121148459931334,0.038146240041516677,0.019072805451553331,-0.037764178358202005},
    {1,0,2e6,1e3,1,.000638922166744306,.0012777253919715736,.0006388612038547122,-.001277249722407947},
    {0,1,2e6,1e4,3,.0028388528457831894,.005675576292507589,.0028377810734108253,-.005667065867136278},
    {0,0,3e6,4e4,4,.002838851119511321,.00567556939253544,.00283777710745971,-.005667045197860543},
    {0,0,1e7,1e4,4,.00006389762853371145,.0001277941842567909,.00006389708960625695,-.00012778989304018662}};
  double source_error=0;
  for(const auto& r:refs) {
    const auto f=ion_quantum_liquid_jets(r.T,r.rho,mixture(r.X,r.Y3,0),1)[0];
    const double scale=r.A/constants::R_gas;
    const double actual[]{f[0][0]*scale,-f[1][0]*scale,f[0][1]*scale,-(f[1][0]+f[2][0])*scale};
    const double expected[]{r.f,r.u,r.p,r.cv};
    for(int k=0;k<4;++k)source_error=std::max(source_error,std::abs(actual[k]/expected[k]-1));
  }
  check(source_error<2e-8,"pure-ion F/U/P/Cv agree with independent Fortran source",source_error);
  double thermal_error=0,composition_error=0;
  auto c=mixture(.10,.08,.13);const double T=3e6,rho=4e4,h=2e-5;
  const auto f=ion_quantum_liquid_jets(T,rho,c);
  // Differentiate every available channel with respect to both thermal
  // coordinates, including chemical-potential and Hessian responses.
  for(unsigned k=0;k<10;++k) {
    const unsigned order=k==0?0:(k<4?1:2);
    for(unsigned direction=0;direction<2;++direction) {
      const auto p=ion_quantum_liquid_jets(T*std::exp(direction==0?h:0),rho*std::exp(direction==1?h:0),c);
      const auto m=ion_quantum_liquid_jets(T*std::exp(direction==0?-h:0),rho*std::exp(direction==1?-h:0),c);
      for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j+order<3;++j) {
        const double fd=(p[k][i][j]-m[k][i][j])/(2*h);
        const double exact=f[k][i+(direction==0)][j+(direction==1)];
        thermal_error=std::max(thermal_error,std::abs(fd-exact)/std::max(std::abs(f[0][0][0]),std::abs(exact)));
      }
    }
  }
  auto perturb=[&](Composition in,unsigned k,double dx) {
    if(k<2)in.X[k]+=dx;
    else {const double z=in.Z();for(unsigned i=METAL_BEGIN;i<METAL_END;++i)in.X[i]*=(z+dx)/z;}
    in.X[2]-=dx;return in;
  };
  for(unsigned a=0;a<3;++a) {
    const auto p=ion_quantum_liquid_jets(T,rho,perturb(c,a,h));
    const auto m=ion_quantum_liquid_jets(T,rho,perturb(c,a,-h));
    for(unsigned i=0;i<3;++i)for(unsigned j=0;i+j<=2;++j)
      composition_error=std::max(composition_error,std::abs((p[0][i][j]-m[0][i][j])/(2*h)-f[1+a][i][j])/std::abs(f[0][0][0]));
    for(unsigned b=0;b<3;++b) {
      const auto lo=std::min(a,b),hi=std::max(a,b);const unsigned ch=4+lo*3-lo*(lo-1)/2+hi-lo;
      for(unsigned i=0;i<2;++i)for(unsigned j=0;i+j<=1;++j)
        composition_error=std::max(composition_error,std::abs((p[1+b][i][j]-m[1+b][i][j])/(2*h)-f[ch][i][j])/std::abs(f[0][0][0]));
    }
  }
  check(thermal_error<2e-7,"all thermal derivatives follow the same free energy",thermal_error);
  check(composition_error<2e-7,"mixture forces and enthalpies include electron-density response",composition_error);
  const auto hot=ion_quantum_liquid_jets(T*10,rho,c,1)[0];
  // Weak-quantum free energy varies as T^-2 and obeys E=2F, P=rho F.
  check(std::abs(-hot[1][0]/hot[0][0]-2)<1e-4 && std::abs(hot[0][1]/hot[0][0]-1)<1e-4,
        "Wigner-Kirkwood limit recovered for an isotope/metal mixture");
  bool refused=false;try{(void)ion_quantum_liquid_jets(1e4,1e5,c);}catch(const std::domain_error&){refused=true;}
  check(refused,"strongly quantum and freezing regime requires further physics");
  refused=false;try{(void)ion_quantum_liquid_jets(1e4,1,c);}catch(const std::domain_error&){refused=true;}
  check(refused,"significant partial-ionization correction refused");
  const auto dense_h=mixture(.99,.003,0);
  const auto cool=ion_quantum_liquid_jets(250000,90,dense_h,1)[0];
  const auto cool_p=ion_quantum_liquid_jets(250000*std::exp(h),90,dense_h,1)[0];
  const auto cool_m=ion_quantum_liquid_jets(250000*std::exp(-h),90,dense_h,1)[0];
  check(cool[0][0]>0 && std::abs((cool_p[0][0]-cool_m[0][0])/(2*h)/cool[1][0]-1)<2e-7,
        "pressure-ionized hydrogen keeps the same potential and thermal response below300kK");
  refused=false;try{(void)ion_quantum_liquid_jets(250000,90,mixture(.9,.003,0));}
  catch(const std::domain_error&){refused=true;}
  check(refused,"dense hydrogen extension does not admit unsupported helium-rich material");
  const auto envelope_h=mixture(.988,.0038,0);
  const auto envelope=ion_quantum_liquid_jets(150000,40,envelope_h,1)[0];
  const auto envelope_p=ion_quantum_liquid_jets(150000*std::exp(h),40,envelope_h,1)[0];
  const auto envelope_m=ion_quantum_liquid_jets(150000*std::exp(-h),40,envelope_h,1)[0];
  check(envelope[0][0]>0 && std::abs((envelope_p[0][0]-envelope_m[0][0])/(2*h)/envelope[1][0]-1)<2e-7,
        "cool hydrogen retains the same free energy and thermal derivative");
  refused=false;try{(void)ion_quantum_liquid_jets(100000,9,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen retains its assessed density boundary");
  refused=false;try{(void)ion_quantum_liquid_jets(90000,22,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen retains its assessed temperature boundary");
  refused=false;try{(void)ion_quantum_liquid_jets(100000,150,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen retains its smaller plasma-temperature bound");
  refused=false;try{(void)ion_quantum_liquid_jets(150000,40,mixture(.97,.0038,0));}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen approximation does not extend to helium-rich layers");
  refused=false;try{(void)ion_quantum_liquid_jets(1e6,3e4,mixture(.01,0,.98));}catch(const std::domain_error&){refused=true;}
  check(refused,"strong ionic coupling is not treated as an assessed liquid");
  const auto dilute=ion_quantum_liquid_jets(3000,1e-7,c,1)[0];
  check(dilute[0][0]>0 && dilute[0][0]<1,"dilute cool material has negligible continuous correction");
  return failed?1:0;
 }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}
}

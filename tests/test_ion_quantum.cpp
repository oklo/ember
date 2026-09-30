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
    {0,0,510000,100000,4,0.23815767721227404,0.46229481396279182,0.23108161861309917,-0.40989100858630345},
    {0,1,550000,100000,3,0.35837431682064352,0.68579776391843983,0.34270014033333823,-0.57383920637448582},
    {0,0,800000,50000,4,0.049594653597738458,0.098547211813387658,0.049271160919621339,-0.096016084818755676},
    {0,0,900000,50000,4,.039239300629936025,.07807529008828462,.03903611044704816,-.0764806265431045},
    {0,1,1050000,50000,3,.051167433198013895,.10164473765723248,.05081875065165928,-.09892582855332743},
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
  const double h=2e-5;
  for(double T:{3e6,1.05e6,8e5,5.1e5,4.5e5}) for(double rho:{8500.,4e4}) {
  if(T<5e5 && rho>1e4)continue;
  auto c=T<5e5?mixture(.005,.001,.0002):mixture(.02,.02,.13);
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
  }
  auto c=mixture(.10,.08,.13);const double T=3e6,rho=4e4;
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
  for(auto state:{std::array<double,4>{399999,8500,.005,.0002},
                  std::array<double,4>{450000,8500,.02,.0002},
                  std::array<double,4>{450000,8500,.005,.002},
                  std::array<double,4>{450000,20000,.005,.0002}}) {
    refused=false;
    try{(void)ion_quantum_liquid_jets(state[0],state[1],mixture(state[2],.001,state[3]));}
    catch(const std::domain_error&){refused=true;}
    check(refused,"cool helium extension retains temperature, density and composition limits");
  }
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
  refused=false;try{(void)ion_quantum_liquid_jets(49999,10,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen retains its assessed temperature boundary");
  const auto colder=ion_quantum_liquid_jets(86000,23.6,envelope_h,1)[0];
  const auto colder_p=ion_quantum_liquid_jets(86000*std::exp(h),23.6,envelope_h,1)[0];
  const auto colder_m=ion_quantum_liquid_jets(86000*std::exp(-h),23.6,envelope_h,1)[0];
  check(colder[0][0]>0 && std::abs((colder_p[0][0]-colder_m[0][0])/(2*h)/colder[1][0]-1)<2e-7,
        "cooling envelope keeps the quantum potential derivative");
  check(ion_quantum_liquid_jets(50000,10,dense_h,1)[0][0][0]>0,
        "assessed 50 kK hydrogen state is supported");
  const auto warm_dense=ion_quantum_liquid_jets(250000,350,dense_h,1)[0];
  const auto warm_dense_p=ion_quantum_liquid_jets(250000*std::exp(h),350,dense_h,1)[0];
  const auto warm_dense_m=ion_quantum_liquid_jets(250000*std::exp(-h),350,dense_h,1)[0];
  check(warm_dense[0][0]>0 && std::abs((warm_dense_p[0][0]-warm_dense_m[0][0])/(2*h)/warm_dense[1][0]-1)<2e-7,
        "denser pressure-ionized hydrogen retains the free-energy derivative");
  refused=false;try{(void)ion_quantum_liquid_jets(299000,501,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cooling hydrogen retains its upper density boundary");
  refused=false;try{(void)ion_quantum_liquid_jets(199999,201,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"denser hydrogen extension retains its warmer temperature boundary");
  refused=false;try{(void)ion_quantum_liquid_jets(100000,150,dense_h);}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen retains its smaller plasma-temperature bound");
  refused=false;try{(void)ion_quantum_liquid_jets(150000,40,mixture(.97,.0038,0));}
  catch(const std::domain_error&){refused=true;}
  check(refused,"cool hydrogen approximation does not extend to helium-rich layers");
  refused=false;try{(void)ion_quantum_liquid_jets(1e6,3e4,mixture(.01,0,.98));}catch(const std::domain_error&){refused=true;}
  check(refused,"strong ionic coupling is not treated as an assessed liquid");
  refused=false;try{(void)ion_quantum_liquid_jets(900000,50000,mixture(.1,0,0));}
  catch(const std::domain_error&){refused=true;}
  check(refused,"larger quantum ratio requires a helium-rich core");
  refused=false;try{(void)ion_quantum_liquid_jets(499999,50000,mixture(0,0,0));}
  catch(const std::domain_error&){refused=true;}
  check(refused,"helium core retains assessed temperature floor");
  refused=false;try{(void)ion_quantum_liquid_jets(900000,50000,mixture(0,0,.17));}
  catch(const std::domain_error&){refused=true;}
  check(refused,"helium core retains the assessed metal fraction");
  check(ion_quantum_liquid_jets(70000,39,mixture(.98,.01,0),1)[0][0][0]>0,
        "cool hydrogen supports the independently checked moderate quantum range");
  const auto dilute=ion_quantum_liquid_jets(3000,1e-7,c,1)[0];
  check(dilute[0][0]>0 && dilute[0][0]<1,"dilute cool material has negligible continuous correction");
  return failed?1:0;
 }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}
}

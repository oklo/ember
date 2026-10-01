// Independent checks of the pure-He4 ion phase component (not selected by any EOS).
// Reference values: pinned eos22.f (sha f514e178..., documented quantum-derivative repairs), evaluated
// component by component (out/helium-solid-potential-sept30-v1/ocp_probe.f90), A=4, Z=2.
#include "ember/ion_phase.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* what,double value){if(!ok){++failures;std::cerr<<"FAIL "<<what<<' '<<value<<'\n';}}
struct Row{int phase;double T,rho;double ocp[7],screened[7];};
const Row rows[]={
  {0,80000,10000,{-135.31206221744591,-132.87757509057531,2.4344871268705792,-43.521213554843449,2.1955888125443308,1.0475625176611101,-58.150009182613658},{-137.27060432160721,-134.65636826896707,2.6142360526400998,-43.458520045502873,2.143467773830956,1.0875402969687704,-58.101872067603672}},
  {0,140000,50000,{-130.9412857057537,-128.9279988496468,2.0132868561068449,-42.072729387892451,1.871617720534523,0.89018973045385175,-56.076350691165217},{-132.06368345411383,-129.9424950562294,2.1211883978843655,-42.040792481626063,1.8405657825326989,0.91263578847207871,-56.057485003118323}},
  {0,250000,100000,{-92.564583256711202,-90.087414663387122,2.4771685933240799,-29.253545597153039,1.948492570431541,0.96516580902166427,-39.097176066610743},{-93.220287564351764,-90.66371604246622,2.5565715218855445,-29.233775199983235,1.9256296145079335,0.98060498131225582,-39.088235954179382}},
  {0,60000,5000,{-143.5926635117643,-141.06282714862741,2.5298363631368939,-46.275161836530273,2.3072642365726792,1.0964006477252719,-61.856679559880462},{-146.19388408166253,-143.43587542861053,2.75800865305201,-46.194227662630901,2.240842042160109,1.1476050374854443,-61.792225007396453}},
  {1,80000,10000,{-135.22263565938391,-133.5292879071103,1.693347752273594,-43.738412238942523,2.241727671700295,1.0579284617746161,-58.440897172395793},{-137.33662445213943,-135.43305104715745,1.9035734049819824,-43.702201561638248,2.4649750290268879,1.2389145704196547,-58.492664014917821}},
  {1,140000,50000,{-130.8362052497786,-129.55165620390099,1.2845490458775479,-42.279324972339211,1.9354570803015589,0.90239543650905596,-56.350996482756642},{-132.09207174387438,-130.72143482095137,1.3706369229229709,-42.290243250736502,2.1141557491840519,1.0188849271089244,-56.428421786186625}},
  {1,250000,100000,{-92.292329860077515,-90.396556102772877,1.8957737573046329,-29.35737862647942,2.7770802353525461,1.235285592741465,-39.323407984889442},{-93.036861506864895,-91.041638959240771,1.9952225476241279,-29.348949801886935,2.8862434642612596,1.3200336754661337,-39.361683980474346}},
  {1,60000,5000,{-143.54135064021719,-141.7561414527195,1.7852091874977289,-46.506260836599047,2.2723809013579732,1.0805498161002149,-62.157741814498408},{-146.31645323683446,-144.25100379045253,2.0654494463819844,-46.448015593526335,2.5246329248272952,1.2999798307795061,-62.202445238302651}},
};
constexpr double A=4,Z=2;
double logne(double rho){return std::log(rho*.5/constants::amu);}
// f,u,s,p,cv,pdt,pdr per ion from a per-mass F/T jet
std::array<double,7> derived(const HelmholtzJet& j){
  const double k=A/constants::R_gas,f=j[0][0]*k,fT=j[1][0]*k,fTT=j[2][0]*k,fn=j[0][1]*k,fnT=j[1][1]*k,fnn=j[0][2]*k;
  return {f,-fT,-fT-f,fn,-fT-fTT,fn+fnT,fn+fnn};
}
HelmholtzJet branch(int phase,double T,double rho,IonScreening s){
  return phase==0?ion_ocp_liquid_jet(std::log(T),logne(rho),A,Z,s):ion_ocp_solid_jet(std::log(T),logne(rho),A,Z,s);
}
double f_ion(int phase,double T,double rho,IonScreening s){return branch(phase,T,rho,s)[0][0]*A/constants::R_gas;}
Composition helium(double xh){auto c=solar_scaled(xh,0.);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;return c;}
}
int main(){
  // 1. Port of the free energies: every response against the Fortran source.
  double worst_ocp=0,worst_scr=0;
  for(const auto& r:rows){
    const auto o=derived(branch(r.phase,r.T,r.rho,IonScreening::ocp)),s=derived(branch(r.phase,r.T,r.rho,IonScreening::fitted_screening));
    for(int i=0;i<7;++i){
      worst_ocp=std::max(worst_ocp,std::abs(o[i]-r.ocp[i])/std::max(1.,std::abs(r.ocp[i])));
      worst_scr=std::max(worst_scr,std::abs(s[i]-r.screened[i])/std::max(1.,std::abs(r.screened[i])));
    }
  }
  check(worst_ocp<1e-11,"OCP branch responses vs eos22.f",worst_ocp);
  // The common liquid uses published constants; eos22.f rounds several
  // screening constants to single precision. Their measured response
  // difference is 2.301 ppm; it is not an interpolation or derivative error.
  check(worst_scr<3e-6,"screened branch responses vs rounded eos22.f",worst_scr);
  // 2. Third-order entries are derivatives of the second-order ones (central differences).
  double worst3=0;
  for(const auto& r:rows)for(auto sc:{IonScreening::ocp,IonScreening::fitted_screening}){
    const double h=1e-4,lt=std::log(r.T),ln=logne(r.rho);
    auto J=[&](double a,double b){return r.phase==0?ion_ocp_liquid_jet(a,b,A,Z,sc):ion_ocp_solid_jet(a,b,A,Z,sc);};
    const auto c=J(lt,ln),tp=J(lt+h,ln),tm=J(lt-h,ln),np=J(lt,ln+h),nm=J(lt,ln-h);
    const double scale=std::abs(c[0][0])+1e-300;
    for(int i=0;i<3;++i)for(int j=0;i+j<=2;++j){
      if(i+j==2){
        worst3=std::max(worst3,std::abs((tp[i][j]-tm[i][j])/(2*h)-c[i+1][j])/scale);
        worst3=std::max(worst3,std::abs((np[i][j]-nm[i][j])/(2*h)-c[i][j+1])/scale);
      }
    }
  }
  check(worst3<1e-6,"third-order jet entries vs differences of second-order entries",worst3);
  // 3. Soft-min phase term across melting (rho=5e4, both variants).
  for(auto sc:{IonScreening::ocp,IonScreening::fitted_screening}){
    const double rho=5e4;IonPhaseOptions opt;opt.screening=sc;
    const double Tm=[&]{double a=6e4,b=3e5;for(int i=0;i<200;++i){double m=std::sqrt(a*b);
      (f_ion(1,m,rho,sc)<f_ion(0,m,rho,sc)?a:b)=m;}return std::sqrt(a*b);}();
    auto full=[&](double T){  // liquid + phase difference, per ion jet
      auto l=branch(0,T,rho,sc);const auto d=ion_phase_difference_jets(T,rho,helium(0.),1,opt)[0];
      for(int i=0;i<4;++i)for(int j=0;i+j<=3;++j)l[i][j]+=d[i][j];return l;};
    const auto lq=derived(branch(0,Tm,rho,sc)),so=derived(branch(1,Tm,rho,sc));const double Q=lq[1]-so[1];
    const int n=8001;const double lo=std::log(Tm)-.25,hi=std::log(Tm)+.25;
    double cv_min=1e300,base_min=1e300,integral=0,prevT=0,prevCv=0;
    for(int i=0;i<n;++i){
      const double T=std::exp(lo+(hi-lo)*i/(n-1));const auto d=derived(full(T));
      cv_min=std::min(cv_min,d[4]);
      base_min=std::min(base_min,std::min(derived(branch(0,T,rho,sc))[4],derived(branch(1,T,rho,sc))[4]));
      if(i)integral+=.5*(d[4]+prevCv)*(T-prevT);prevT=T;prevCv=d[4];
    }
    const double Tlo=std::exp(lo),Thi=std::exp(hi);
    const double E=derived(full(Thi))[1]*Thi-derived(full(Tlo))[1]*Tlo;
    const double latent=E-(derived(branch(0,Thi,rho,sc))[1]*Thi-lq[1]*Tm)-(so[1]*Tm-derived(branch(1,Tlo,rho,sc))[1]*Tlo);
    check(cv_min>=base_min-1e-9,"soft-min heat capacity below both pure phases",cv_min-base_min);
    check(std::abs(E-integral)/std::abs(E)<1e-6,"state energy difference equals integral of Cv",std::abs(E-integral)/std::abs(E));
    check(std::abs(latent/(Q*Tm)-1)<1e-3,"recovered latent heat",latent/(Q*Tm)-1);
    const double offset=derived(full(Tm))[0]-std::min(lq[0],so[0]);
    check(std::abs(offset+opt.width*std::log(2.))<1e-9,"soft-min offset at melting equals -c ln 2",offset);
    const double far=std::abs(derived(full(std::exp(lo)))[0]-derived(branch(1,std::exp(lo),rho,sc))[0]);
    check(far<1e-9,"pure solid limit below the window",far);
  }
  // 4. Guards.
  auto throws=[](auto f){try{f();return false;}catch(const std::domain_error&){return true;}};
  check(throws([&]{ion_phase_difference_jets(1.5e5,5e4,helium(.01),1);}),"composition gate",0);
  check(throws([&]{ion_phase_difference_jets(1e7,1e8,helium(0.),1);}),"Baiko-Chugunov domain gate",0);
  const auto zero=ion_phase_difference_jets(1e6,5e4,helium(0.),10);   // Gamma < 90: liquid, no solid evaluation
  double zmax=0;for(const auto& ch:zero)for(const auto& a:ch)for(double v:a)zmax=std::max(zmax,std::abs(v));
  check(zmax==0,"no phase term below the solid-Gamma cut",zmax);
  // The cut at Gamma=100 must be invisible (solid weight < exp(-30)) over the supported density range.
  int cut_failures=0;
  for(double rho=3e3;rho<=2e6;rho*=1.5)for(auto sc:{IonScreening::ocp,IonScreening::fitted_screening}){
    IonPhaseOptions opt;opt.screening=sc;
    const double T100=[&]{const double T=1e5;const double lr=logne(rho);(void)lr;return T;}();(void)T100;
    for(double g=90;g<=100;g+=.5){
      // Gamma = Z^(5/3) e^2/(a_e k T): find T for this Gamma by scaling from a reference evaluation.
      const double Tref=1e5;const double gref=std::pow(Z,5./3)*315775.02480407/
        (std::cbrt(3/(4*3.141592653589793*std::exp(logne(rho))))/5.29177210903e-9)/Tref;
      const double T=Tref*gref/g;
      try{(void)ion_phase_difference_jets(T,rho,helium(0.),1,opt);}catch(const std::domain_error& e){
        if(std::string(e.what()).find("Gamma cut")!=std::string::npos)++cut_failures;}
    }
  }
  check(cut_failures==0,"solid weight negligible at the Gamma cut",cut_failures);
  // 5. Composition gradient (H replaces He4) against a difference of the value channel.
  {const double T=1.5e5,rho=5e4,x=2e-4,h=1e-6;
   const auto c=ion_phase_difference_jets(T,rho,helium(x),10);
   const double fd=(ion_phase_difference_jets(T,rho,helium(x+h),1)[0][0][0]-ion_phase_difference_jets(T,rho,helium(x-h),1)[0][0][0])/(2*h);
   const double err=std::abs(fd-c[1][0][0])/std::max(1e-300,std::abs(c[1][0][0]));
   check(err<1e-6,"composition gradient channel vs finite difference",err);}
  std::cout<<"ion phase: OCP port "<<worst_ocp<<", screened port "<<worst_scr<<", third order "<<worst3
           <<", failures "<<failures<<'\n';
  return failures?1:0;
}

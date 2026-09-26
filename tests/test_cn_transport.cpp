#include "ember/cn_burning.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace ember;
namespace {
using Six=std::array<double,6>;
Six physical(const Composition& c) {
  const auto q=cn_physical_ledger(c,*c.cn_molality);const auto& y=*c.cn_molality;
  return {c.X[0],c.X[1],q.helium4,12*y[0],13*y[1],14*y[2]};
}
} // namespace
int main() {
  int checks=0,failures=0;
  auto check=[&](bool ok,const char* label,double error=0.) {
    ++checks;if(!ok){++failures;std::printf("FAIL %s %.8g\n",label,error);}
  };
  auto rejects=[](auto f){try{f();}catch(const std::exception&){return true;}return false;};
  PPCNNetwork nuclear;Model m;m.M=10;m.m={1,2,4,7,10};
  for(int i=0;i<5;++i) {
    auto c=solar_scaled(.15+.1*i,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=i==0?0:.001*i;c.X[2]-=c.X[1];const auto y=initial_gs98_cn(c);
    c.cn_molality=CNAbundances{y[0]*(1-.2*i),.03*i*y[0],y[2]+.17*i*y[0]};
    m.comp.push_back(c);m.y.push_back({std::log(1.+.1*i),0,std::log(1e4),1});
  }
  const auto w=nodal_mass_weights(m);MixingRegions single{{0,1},{1,2},{2,3},{3,4},{4,5}};
  double max_residual=0,max_balance=0,max_mass_error=0,max_stiff_error=0;
  const auto original=m.comp;
  for(bool gradient:{false,true})for(bool warm:{false,true}) {
    m.comp=original;
    if(gradient)for(double& y:*m.comp[0].cn_molality)y*=.7;
    for(std::size_t i=0;i<m.size();++i) {
      m.y[i].lnT=std::log(warm?1.2e7+1e6*static_cast<double>(i):1e4);m.y[i].lnrho=std::log(warm?300.:1.);
    }
    const double dt=warm?1e13:2.;
    for(int configuration=0;configuration<4;++configuration) {
      const auto regions=configuration==3?MixingRegions{{0,2},{2,3},{3,5}}:single;
      std::vector<double> D(4);
      for(std::size_t i=0;i<D.size();++i) {
        const double r=.5*(m.r(i)+m.r(i+1)),rho=.5*(m.rho(i)+m.rho(i+1));
        const double a=4*M_PI*r*r*rho;
        const double g=configuration==1?1e30:(configuration==2 && i==2?0:.01*static_cast<double>(i+1));
        D[i]=g*m.M/dt*(m.m[i+1]-m.m[i])/(a*a);
      }
      const auto after=burn_and_transport(m,m,nuclear,regions,D,dt,1e-14);
      Six global{};double heat=0,rest=0;
      std::vector<Six> values,old,sources;
      for(std::size_t i=0;i<m.size();++i) {
        values.push_back(physical(after[i]));old.push_back(physical(m.comp[i]));
        const auto p=nuclear.pp().eval(m.T(i),m.rho(i),after[i]);
        const auto c=nuclear.cn().response(m.T(i),m.rho(i),after[i],*after[i].cn_molality).physical.state;
        Six source{};for(int k=0;k<6;++k)source[k]=p.dXdt[k]+c.dXdt[k];sources.push_back(source);
        for(int k=0;k<6;++k) {
          check(values.back()[k]>=0,"physical isotope stays nonnegative");
          global[k]+=w[i]/m.M*(values.back()[k]-old.back()[k]-dt*source[k]);
          rest-=w[i]/dt*(values.back()[k]-old.back()[k])*(nuclides[k].A/mass_numbers[k]-1)*constants::c*constants::c;
        }
        heat+=w[i]*(p.eps+c.eps+p.eps_neutrino+c.eps_neutrino);
      }
      for(double b:global){max_balance=std::max(max_balance,std::abs(b));check(std::abs(b)<2e-14,"global physical-isotope conservation",b);}
      if(warm){const double e=std::abs(rest/heat-1);max_mass_error=std::max(max_mass_error,e);check(e<2e-6,"physical nuclear mass-energy balance",e);}
      if(configuration!=1)for(auto [begin,end]:regions) {
        Six residual{};
        for(std::size_t i=begin;i<end;++i)for(int k=0;k<6;++k)
          residual[k]+=w[i]/m.M*(values[i][k]-old[i][k]-dt*sources[i][k]);
        for(int side:{-1,1}) {
          if((side<0 && begin==0) || (side>0 && end==m.size()))continue;
          const std::size_t face=side<0?begin-1:end-1;
          const double r=.5*(m.r(face)+m.r(face+1)),rho=.5*(m.rho(face)+m.rho(face+1));
          const double a=4*M_PI*r*r*rho,g=dt/m.M*a*a*D[face]/(m.m[face+1]-m.m[face]);
          const auto inside=side<0?face+1:face,outside=side<0?face:face+1;
          for(int k=0;k<6;++k)residual[k]+=g*(values[inside][k]-values[outside][k]);
        }
        for(double e:residual){max_residual=std::max(max_residual,std::abs(e));check(std::abs(e)<2e-14,"independent region flux equations",e);}
        for(std::size_t i=begin+1;i<end;++i)check(after[i]==after[begin],"convective reservoir homogeneous");
      }
      if(configuration==1) {
        const auto fully=burn_and_mix(m,m,nuclear,{{0,m.size()}},dt,1e-14);
        for(std::size_t i=0;i<m.size();++i)for(int k=0;k<6;++k) {
          const double e=std::abs(values[i][k]-physical(fully[i])[k]);max_stiff_error=std::max(max_stiff_error,e);
          check(e<2e-13,"stiff transport matches independent burning/convection",e);
        }
      }
    }
    check(burn_and_transport(m,m,nuclear,single,{0,0,0,0},dt,1e-14)==burn_and_mix(m,m,nuclear,single,dt,1e-14),"zero diffusion exact local path");
  }
  auto bad=m;bad.comp[0].cn_molality.reset();
  check(rejects([&]{burn_and_transport(bad,bad,nuclear,single,{1,1,1,1},2);}),"missing CN rejected");
  check(rejects([&]{burn_and_transport(m,m,nuclear,{{1,5}},{1,1,1,1},2);}),"bad partition rejected");
  check(rejects([&]{burn_and_transport(m,m,nuclear,single,{1,-1,1,1},2);}),"negative diffusivity rejected");
  std::printf("%d checks, %d failures; max region %.4g global %.4g mass-energy %.4g stiff %.4g\n",
    checks,failures,max_residual,max_balance,max_mass_error,max_stiff_error);
  return failures?1:0;
}

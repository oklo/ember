#include "ember/eos_additive_volume.hpp"
#include "ember/constants.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>
using namespace ember;
void close(double a,double b,double scale,const char* what){
  if(std::abs(a-b)>2e-10*scale)throw std::runtime_error(what);
}
int main(int argc,char** argv){try{
  if(argc==2){
    AdditiveVolumePotential eos(argv[1]);std::cout<<std::setprecision(17);
    double T,rho,X,Y;
    while(std::cin>>T>>rho>>X>>Y){auto q=eos.residual_jets(T,rho,X,Y);for(auto& a:q)for(auto& b:a)for(double v:b)std::cout<<v<<' ';std::cout<<'\n';}
    return 0;
  }
  const auto path=std::filesystem::temp_directory_path()/
      ("ember-additive-volume-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".dat");
  struct Cleanup {std::filesystem::path p;~Cleanup(){std::error_code e;std::filesystem::remove(p,e);}}cleanup{path};
  const double R=constants::R_gas;
  {std::ofstream out(path);out<<std::setprecision(17)<<"EMBER_ADDITIVE_VOLUME_POTENTIAL 1\n";
   for(int a:{1,4}){out<<(a==1?"H":"He")<<" 2 2\n8 16\n-12 8\n";
    std::array<double,16> c{};c[0]=R/a*(-12-1.5*8+double(a));c[1]=R/a;c[4]=-1.5*R/a;
    for(auto v:c)out<<v<<' ';out<<'\n';}}
  AdditiveVolumePotential ideal(path);
  for(double T:{1e4,1e5,1e6})for(double rho:{.01,1.,100.})for(auto c:{std::array{.7,.01},std::array{.992,.003},std::array{1.,0.},std::array{0.,0.}}){
    const auto q=ideal.residual_jets(T,rho,c[0],c[1]);const auto scalar=ideal.residual_jets(T,rho,c[0],c[1],false);
    double he=(1-c[0]+c[1]/3)/4,n=c[0]+he,w=std::log(rho)-1.5*std::log(T);
    close(q[0][0][0],R*(n*w+c[0]+he*(4+std::log(4.))),R*100,"ideal potential");
    close(q[0][0][1],R*n,R,"ideal pressure");close(q[0][1][0],-1.5*R*n,R,"ideal energy");
    close(q[1][0][0],R*(.75*w+1-.25*(4+std::log(4.))),R*100,"hydrogen derivative");
    close(q[2][0][0],R/12*(w+4+std::log(4.)),R*100,"helium-3 derivative");
    close(q[1][0][1],.75*R,R,"pressure composition derivative");close(q[2][1][0],-R/8,R,"energy composition derivative");
    for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j){
      close(q[0][i][j],scalar[0][i][j],R*100,"thermal and composition modes disagree");
      if(i+j>=2)close(q[0][i][j],0.,R,"ideal higher thermal derivative");
    }
    for(int k:{3,4,5})for(auto& a:q[k])for(double v:a)close(v,0.,R,"ideal residual composition Hessian");
  }
  bool rejected=false;try{ideal.residual_jets(1e5,1.,.99,.02);}catch(const std::domain_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("unphysical composition accepted");
  rejected=false;try{ideal.residual_jets(1.,1.,.99,0.);}catch(const std::domain_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("out-of-domain temperature accepted");
  const auto bounded_file=[&](bool invalid){
    std::ofstream out(path);out<<std::setprecision(17)<<"EMBER_ADDITIVE_VOLUME_POTENTIAL 2\n";
    for(int a:{1,4}){
      out<<(a==1?"H":"He")<<" 3 2\n8 12 16\n-12 8\ndensity_bounds\n";
      out<<(invalid?2:-12)<<' '<<(a==1?0:8)<<"\n-12 8\n";
      for(double t:{8.,12.}){
        std::array<double,16> c{};c[0]=R/a*(-12-1.5*t+double(a));c[1]=R/a;c[4]=-1.5*R/a;
        for(auto v:c)out<<v<<' ';out<<'\n';
      }
    }
  };
  bounded_file(false);AdditiveVolumePotential bounded(path);
  for(auto point:{std::array{1e4,.5},std::array{1e6,100.}}){
    const auto a=ideal.residual_jets(point[0],point[1],.7,.01);
    const auto b=bounded.residual_jets(point[0],point[1],.7,.01);
    for(unsigned k=0;k<6;++k)for(unsigned i=0;i<4;++i)for(unsigned j=0;i+j<=3;++j)
      close(a[k][i][j],b[k][i][j],R*100,"density restriction changed supported physics");
  }
  const auto range=bounded.density_range(1e4,.7,0.);
  close(range.max,1/(.7+.3/4),1.,"cold mixture density bound");
  rejected=false;try{bounded.residual_jets(1e4,2.,.7,0.);}catch(const std::domain_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("excluded cold dense state accepted");
  bounded_file(true);rejected=false;
  try{AdditiveVolumePotential invalid(path);}catch(const std::runtime_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("invalid density bounds accepted");
  std::cout<<"additive-volume ideal, isotope, derivative and domain checks passed\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

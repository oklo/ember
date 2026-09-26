#include "ember/material_transport.hpp"
#include "ember/eos_smooth_mixture.hpp"
#include "material_test_eos.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace ember;
namespace {
void read(MaterialVector& a,std::size_t n){for(std::size_t i=0;i<n;++i)std::cin>>a[i];}
void read(MaterialMatrix& a,std::size_t n){for(std::size_t i=0;i<n;++i)read(a[i],n);}
void print(const MaterialVector& a,std::size_t n) {
  std::cout<<'[';for(std::size_t i=0;i<n;++i)std::cout<<(i?",":"")<<a[i];std::cout<<']';
}
void print(const MaterialMatrix& a,std::size_t n) {
  std::cout<<'[';for(std::size_t i=0;i<n;++i){if(i)std::cout<<',';print(a[i],n);}std::cout<<']';
}
void print(const std::vector<MaterialVector>& a,std::size_t n) {
  std::cout<<'[';for(std::size_t i=0;i<a.size();++i){if(i)std::cout<<',';print(a[i],n);}std::cout<<']';
}
Composition composition(double x,double y) {
  auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=y;c.X[2]-=y;return c;
}
}
int main(int argc,char** argv) {
  try {
    std::unique_ptr<SmoothMetalHelmholtzEos> native;
    if(argc==2) {
      const auto start=std::chrono::steady_clock::now();
      native=std::make_unique<SmoothMetalHelmholtzEos>(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
      std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    } else if(argc!=1)return 2;
    std::cout<<std::setprecision(17);std::string mode;
    while(std::cin>>mode) {
      try {
        if(mode=="chain") {
          std::size_t count,n;std::cin>>count>>n;if(count<1 || count>10000 || n<1 || n>3)return 2;
          std::vector<MaterialMatrix> a(count),k(count-1);std::vector<MaterialVector> b(count);
          for(auto& v:a)read(v,n);for(auto& v:k)read(v,n);for(auto& v:b)read(v,n);
          const auto answer=solve_material_chain(a,k,b,n);
          std::cout<<"{\"ok\":true,\"answer\":";print(answer,n);std::cout<<"}\n";
          continue;
        }
        if(mode=="point") {
          double x,y,T,rho,e0,s0;int derivatives;std::cin>>x>>y>>T>>rho>>e0>>s0>>derivatives;
          if(!native)return 2;
          const auto p=native_material_point(*native,T,rho,composition(x,y),e0,s0,derivatives!=0);
          std::cout<<"{\"ok\":true,\"conserved\":";print(p.conserved,3);
          std::cout<<",\"potential\":";print(p.potential,3);std::cout<<",\"capacity\":";print(p.capacity,3);
          std::cout<<",\"primitive_from_conserved\":";print(p.primitive_from_conserved,3);
          std::cout<<",\"entropy\":"<<p.entropy<<"}\n";continue;
        }
        if(mode!="analytic" && mode!="native")return 2;
        std::size_t count,n;double dt,tolerance,e0,s0;int changing;
        std::cin>>count>>n>>dt>>tolerance>>e0>>s0>>changing;
        if(count<1 || count>10000 || n<1 || n>3)return 2;
        AnalyticMaterial ideal;ideal.n=n;ideal.energy_scale=e0;
        if(mode=="analytic") {read(ideal.masses,n);read(ideal.charges,n);read(ideal.binding,n);std::cin>>ideal.degeneracy;}
        ideal.density.resize(count);std::vector<double> mass(count);
        for(auto& v:ideal.density)std::cin>>v;for(auto& v:mass)std::cin>>v;
        std::vector<MaterialVector> old(count),guess(count);std::vector<MaterialMatrix> k(count-1);
        for(auto& v:old)read(v,n);for(auto& v:guess)read(v,n);for(auto& v:k)read(v,n);
        if(!std::cin)return 2;
        MaterialThermodynamics eos=ideal;
        MaterialTransportOptions options;options.components=n;options.tolerance=tolerance;
        if(mode=="native") {
          if(!native || n!=3)return 2;options.composition_sum_limit=.98;
          eos=[&](std::size_t i,const MaterialVector& q,bool derivatives) {
            return native_material_point(*native,std::exp(q[2]),ideal.density[i],composition(q[0],q[1]),e0,s0,derivatives);
          };
        }
        auto faces=[&](std::span<const MaterialVector> q) {
          auto out=k;
          if(changing) {
            if(n!=2)throw std::invalid_argument("Fick limit requires two components");
            for(std::size_t i=0;i+1<count;++i) {
              const double a=q[i][0],b=q[i+1][0],mid=.5*(a+b);
              const double mobility=std::abs(a-b)<1e-7*std::min(mid,1-mid)?mid*(1-mid):
                  (b-a)/(std::log(b)-std::log1p(-b)-std::log(a)+std::log1p(-a));
              out[i][0][0]*=mobility;
            }
          }
          return out;
        };
        const auto start=std::chrono::steady_clock::now();
        const auto a=transport_material(old,guess,mass,dt,eos,faces,options);
        std::cout<<"{\"ok\":true,\"primitives\":";print(a.primitives,n);
        std::cout<<",\"conserved\":";print(a.conserved,n);std::cout<<",\"face_flux\":";print(a.face_flux,n);
        std::cout<<",\"integrated_conservation_error\":";print(a.integrated_conservation_error,n);
        std::cout<<",\"residual\":"<<a.residual<<",\"iterations\":"<<a.iterations
          <<",\"entropy_change\":"<<a.entropy_change<<",\"backward_euler_entropy_bound\":"<<a.backward_euler_entropy_bound
          <<",\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
      } catch(const std::exception& e) {
        std::cout<<"{\"ok\":false,\"error\":"<<std::quoted(e.what())<<"}\n";
      }
    }
    return std::cin.eof()?0:2;
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

#include "ember/species_transport.hpp"
#include "ember/eos_smooth_mixture.hpp"
#include "ember/material_flux.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>

using namespace ember;
namespace {
template<class T> void read(T& a){for(auto& value:a)std::cin>>value;}
void print(const SpeciesVector& a){std::cout<<'['<<a[0]<<','<<a[1]<<']';}
void print(const std::vector<SpeciesVector>& a) {
  std::cout<<'[';for(std::size_t i=0;i<a.size();++i){if(i)std::cout<<',';print(a[i]);}std::cout<<']';
}
Composition composition(double x,double y) {
  auto c=solar_scaled(x,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
  c.X[1]=y;c.X[2]-=y;return c;
}
struct ZeroNuclear final:Nuclear {
  NuclearState eval(double,double,const Composition&)const override{return {};}
  NuclearResponse composition_response(double,double,const Composition&)const override{return {};}
  const char* name()const override{return "no burning";}
};
}
int main(int argc,char** argv) {
  try {
    std::unique_ptr<SmoothMetalHelmholtzEos> eos;
    if(argc==2) {
      const auto start=std::chrono::steady_clock::now();
      eos=std::make_unique<SmoothMetalHelmholtzEos>(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
      std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    }else if(argc!=1)return 2;
    std::cout<<std::setprecision(17);std::string mode;
    const PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);const ZeroNuclear zero;
    while(std::cin>>mode) {
      try {
        if(mode=="chain") {
          std::size_t n;std::cin>>n;if(n<1 || n>10000)return 2;
          std::vector<SpeciesMatrix> e(n),a(n-1),b(n-1);std::vector<SpeciesVector> rhs(n),d(n-1);
          for(auto& block:e)for(auto& row:block)read(row);
          for(auto& block:a)for(auto& row:block)read(row);
          for(auto& block:b)for(auto& row:block)read(row);
          for(auto& v:rhs)read(v);
          for(auto& v:d)read(v);
          if(!std::cin)return 2;
          const auto result=solve_species_flux_chain(e,a,b,rhs,d);
          std::cout<<"{\"ok\":true,\"answer\":";print(result);std::cout<<"}\n"<<std::flush;continue;
        }
        if(mode!="native" || !eos)return 2;
        std::size_t n,nregions;double dt,tolerance,e0,s0;int burning;
        std::cin>>n>>nregions>>dt>>tolerance>>e0>>s0>>burning;
        if(n<2 || n>10000 || nregions<1 || nregions>n || (burning!=0 && burning!=1))return 2;
        MixingRegions regions;std::size_t begin=0;
        for(std::size_t i=0;i<nregions;++i){std::size_t end;std::cin>>end;regions.emplace_back(begin,end);begin=end;}
        Model model;
        for(std::size_t i=0;i<n;++i) {
          double mass,rho,T,x,y;std::cin>>mass>>rho>>T>>x>>y;
          model.m.push_back(mass);model.y.push_back({0,std::log(rho),std::log(T),0});model.comp.push_back(composition(x,y));
        }
        model.M=model.m.back();std::vector<MaterialMatrix> k(n-1);
        for(auto& matrix:k)for(auto& row:matrix)read(row);
        if(!std::cin)return 2;
        struct Cache {Composition c;CompositionPotentialResponse p;};
        std::vector<std::optional<Cache>> cache(n);std::size_t queries=0;
        auto potential=[&](std::size_t i,const Composition& c)->const CompositionPotentialResponse& {
          if(!cache[i] || cache[i]->c.X!=c.X) {
            cache[i]=Cache{c,eos->composition_potential(model.T(i),model.rho(i),c)};++queries;
          }
          return cache[i]->p;
        };
        SpeciesFlux flux=[&](std::size_t i,const Composition& left,const Composition& right,bool derivatives) {
          bool active=false;for(const auto& row:k[i])for(double x:row)active|=x!=0;
          if(!active)return SpeciesFaceResponse{};
          const auto& p=potential(i,left);const auto& q=potential(i+1,right);
          SpeciesFaceResponse f;
          f.rate=material_face_flux(k[i],p.gradient,q.gradient,model.T(i),model.T(i+1),1,1,1,e0,s0).species_rate;
          if(derivatives)for(std::size_t j=0;j<2;++j)for(std::size_t col=0;col<2;++col)
            for(std::size_t v=0;v<2;++v) {
              f.dleft[j][col]+=k[i][j][v]*p.hessian[v][col]/s0;
              f.dright[j][col]-=k[i][j][v]*q.hessian[v][col]/s0;
            }
          return f;
        };
        SpeciesTransportOptions options;options.abundance_tolerance=tolerance;
        const Nuclear& nuclear=burning?static_cast<const Nuclear&>(pp):static_cast<const Nuclear&>(zero);
        const auto start=std::chrono::steady_clock::now();
        const auto result=burn_and_diffuse(model,model,nuclear,regions,flux,dt,options);
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::vector<SpeciesVector> final,rate,potential_values;
        std::vector<SpeciesMatrix> hessian,reaction_jacobian;
        for(std::size_t i=0;i<n;++i) {
          const auto& c=result.composition[i];final.push_back({c.X[0],c.X[1]});
          const auto reaction=nuclear.composition_response(model.T(i),model.rho(i),c);
          rate.push_back({reaction.state.dXdt[0],reaction.state.dXdt[1]});
          const auto& p=potential(i,c);potential_values.push_back(p.gradient);hessian.push_back(p.hessian);
          SpeciesMatrix jac{};
          for(std::size_t j=0;j<2;++j)for(std::size_t col=0;col<2;++col)
            jac[j][col]=reaction.d_dXdt_dX[j][col]-reaction.d_dXdt_dX[j][2];
          reaction_jacobian.push_back(jac);
        }
        std::cout<<"{\"ok\":true,\"composition\":";print(final);
        std::cout<<",\"nuclear_source\":";print(rate);
        std::cout<<",\"chemical_potential\":";print(potential_values);
        for(const auto& [name,array]:{std::pair{"chemical_hessian",&hessian},std::pair{"nuclear_jacobian",&reaction_jacobian}}) {
          std::cout<<",\""<<name<<"\":[";
          for(std::size_t i=0;i<array->size();++i) {
            if(i)std::cout<<',';
            std::cout<<'[';print((*array)[i][0]);std::cout<<',';print((*array)[i][1]);std::cout<<']';
          }
          std::cout<<']';
        }
        std::cout<<",\"balance\":";print(result.integrated_balance);
        std::cout<<",\"boundary_flux\":[";
        for(std::size_t i=0;i<result.boundary_fluxes.size();++i){if(i)std::cout<<',';print(result.boundary_fluxes[i].rate);}
        std::cout<<"],\"residual_history\":[";
        for(std::size_t i=0;i<result.residual_history.size();++i)std::cout<<(i?",":"")<<result.residual_history[i];
        std::cout<<"],\"residual\":"<<result.residual<<",\"iterations\":"<<result.iterations
                 <<",\"abundance_correction\":"<<result.abundance_correction
                 <<",\"native_queries\":"<<queries<<",\"elapsed_seconds\":"<<elapsed<<"}\n"<<std::flush;
      }catch(const std::exception& e){std::cout<<"{\"ok\":false,\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
    }
    return std::cin.eof()?0:2;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

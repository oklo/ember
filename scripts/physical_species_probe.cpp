#include "ember/screened_microscopic_transport.hpp"
#include "ember/species_transport.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <optional>

using namespace ember;
namespace {
void vec(const SpeciesVector& a){std::cout<<'['<<a[0]<<','<<a[1]<<']';}
void mat(const SpeciesMatrix& a){std::cout<<'[';vec(a[0]);std::cout<<',';vec(a[1]);std::cout<<']';}
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
  if(argc!=3)return 2;
  try {
    const auto start=std::chrono::steady_clock::now();
    SmoothMetalHelmholtzEos eos(argv[1],HelmholtzTableEos::Mixture::allow_documented_proxy);
    ScreenedCollisionTransport collisions(argv[2]);
    const PPChains pp(PPRates::solar_fusion_ii,PPScreening::salpeter_van_horn);const ZeroNuclear zero;
    std::cerr<<"load_seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::string line;std::cout<<std::setprecision(17);
    while(std::getline(std::cin,line)) {
      try {
        std::istringstream in(line);std::size_t n,nr;double dt,blend;int burning,ions,transport;
        if(!(in>>n>>nr>>dt>>burning>>ions>>blend>>transport) || n<2 || n>1000 || nr<1 || nr>n)return 2;
        MixingRegions regions;std::size_t begin=0;
        for(std::size_t r=0;r<nr;++r){std::size_t end;in>>end;regions.emplace_back(begin,end);begin=end;}
        Model model;
        for(std::size_t i=0;i<n;++i) {
          double mass,radius,rho,T,x,y;in>>mass>>radius>>rho>>T>>x>>y;
          model.m.push_back(mass);model.y.push_back({std::log(radius),std::log(rho),std::log(T),0});model.comp.push_back(composition(x,y));
        }
        if(!in)return 2;model.M=model.m.back();
        if((blend<0 && blend!=-1) || blend>1)throw std::domain_error("invalid guess blend");
        const bool default_guess=blend==-1;const double fraction=default_guess?0.:blend;
        const auto weights=nodal_mass_weights(model);
        const Nuclear& nuclear=burning?static_cast<const Nuclear&>(pp):static_cast<const Nuclear&>(zero);
        const auto begin_time=std::chrono::steady_clock::now();
        auto seed=burn_and_mix(model,model,nuclear,regions,dt,1e-14);
        std::array<long double,2> average{};
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<2;++k)average[k]+=static_cast<long double>(weights[i]/model.M)*seed[i].X[k];
        std::array<bool,2> active{average[0]>0,average[1]>0};
        SpeciesVector guess_balance{};
        for(std::size_t i=0;i<n;++i) {
          for(std::size_t k=0;k<2;++k) {
            const double old=seed[i].X[k];seed[i].X[k]=(1-fraction)*old+fraction*static_cast<double>(average[k]);
            guess_balance[k]+=weights[i]/model.M*(seed[i].X[k]-old);
          }
          for(std::size_t k=3;k<NSPEC;++k)seed[i].X[k]=model.comp.front().X[k];
          seed[i].X[2]=1-seed[i].Z()-seed[i].X[0]-seed[i].X[1];
        }
        ScreenedMicroscopicTransport physics(eos,collisions,ions,2e6,active);
        std::size_t evaluations=0;
        SpeciesFlux flux=[&](std::size_t i,const Composition& a,const Composition& b,bool derivatives) {
          if(!transport)return SpeciesFaceResponse{};
          ++evaluations;
          return microscopic_face(physics,i,model.m[i],model.m[i+1],model.y[i],a,model.y[i+1],b,derivatives).species;
        };
        SpeciesTransportOptions options;if(!default_guess)options.initial_guess=seed;options.seed_present_species=true;
        const auto result=burn_and_diffuse(model,model,nuclear,regions,flux,dt,options);
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin_time).count();
        std::vector<SpeciesVector> rates;std::vector<SpeciesMatrix> jac;
        double old_phi=0,new_phi=0;
        for(std::size_t i=0;i<n;++i) {
          const auto r=nuclear.composition_response(model.T(i),model.rho(i),result.composition[i]);rates.push_back({r.state.dXdt[0],r.state.dXdt[1]});
          SpeciesMatrix j{};for(std::size_t k=0;k<2;++k)for(std::size_t l=0;l<2;++l)j[k][l]=r.d_dXdt_dX[k][l]-r.d_dXdt_dX[k][2];jac.push_back(j);
          old_phi+=weights[i]/model.M*eos.regular_composition_potential(model.T(i),model.rho(i),model.comp[i]).phi;
          new_phi+=weights[i]/model.M*eos.regular_composition_potential(model.T(i),model.rho(i),result.composition[i]).phi;
        }
        std::vector<MicroscopicFaceResponse> final_faces;
        for(const auto& f:result.boundary_fluxes) {
          const auto i=f.face;
          final_faces.push_back(transport?microscopic_face(physics,i,model.m[i],model.m[i+1],model.y[i],result.composition[i],model.y[i+1],result.composition[i+1],true):MicroscopicFaceResponse{});
        }
        // All operations which may reject a state precede output.
        std::cout<<"{\"automatic_seed\":"<<result.interior_guess_adjusted<<",\"active\":["<<active[0]<<','<<active[1]<<"],\"composition\":[";
        for(std::size_t i=0;i<n;++i){if(i)std::cout<<',';vec({result.composition[i].X[0],result.composition[i].X[1]});}
        std::cout<<"],\"nuclear\":[";for(std::size_t i=0;i<n;++i){if(i)std::cout<<',';vec(rates[i]);}
        std::cout<<"],\"nuclear_jacobian\":[";for(std::size_t i=0;i<n;++i){if(i)std::cout<<',';mat(jac[i]);}
        std::cout<<"],\"faces\":[";
        for(std::size_t k=0;k<final_faces.size();++k) {
          if(k)std::cout<<',';const auto& f=final_faces[k];
          std::cout<<"{\"index\":"<<result.boundary_fluxes[k].face<<",\"rate\":";vec(f.species.rate);
          std::cout<<",\"dleft\":";mat(f.species.dleft);std::cout<<",\"dright\":";mat(f.species.dright);
          std::cout<<",\"carried\":"<<f.carried_luminosity<<",\"conductivity\":"<<f.conductivity<<'}';
        }
        std::cout<<"],\"balance\":";vec(result.integrated_balance);std::cout<<",\"guess_balance\":";vec(guess_balance);
        std::cout<<",\"old_phi\":"<<old_phi<<",\"new_phi\":"<<new_phi<<",\"seconds\":"<<elapsed<<",\"face_evaluations\":"<<evaluations
          <<",\"iterations\":"<<result.iterations<<",\"abundance_correction\":"<<result.abundance_correction<<",\"residual\":"<<result.residual
          <<",\"residual_history\":[";
        for(std::size_t i=0;i<result.residual_history.size();++i)std::cout<<(i?",":"")<<result.residual_history[i];
        std::cout<<"]}\n"<<std::flush;
      }catch(const std::exception& e){std::cout<<"{\"error\":"<<std::quoted(e.what())<<"}\n"<<std::flush;}
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

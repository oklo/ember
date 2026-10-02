#include "ember/fingering_transport.hpp"
#include "ember/opacity_radiative.hpp"
#include "ember/convection.hpp"
#include "ember/opacity_conductive_interior.hpp"
#include "ember/conduction_table.hpp"
#include "ember/convective_material_heat.hpp"
#include "ember/envelope_transport.hpp"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
using namespace ember;

int main(int argc,char**argv)try {
 const char* configuration=argc>1?argv[1]:std::getenv("EMBER_LIFETIME_CONFIG");
 if(!configuration){std::cout<<"requires a lifetime configuration with runtime tables\n";return 77;}
 std::map<std::string,std::string> cfg;std::ifstream input(configuration);std::string key,value;
 while(input>>key>>std::quoted(value))cfg[key]=value;
 ColdHeliumOptions o;o.mixture_phase=true;o.liquid_continuation_gamma=200;o.dense_transition=true;o.deep_solid=true;
 VariableMetalHelmholtzEos eos(cfg.at("eos"),HelmholtzTableEos::Mixture::allow_documented_proxy,
  VariableMetalHelmholtzEos::LowMetalInterpolation::quadratic,{},true,o,
  VariableMetalHelmholtzEos::IsotopeInterpolation::number_density);
 RadiativeOpacity raw({cfg.at("opacity_low"),cfg.at("opacity_warm"),cfg.at("opacity_bridge"),cfg.at("opacity_hot"),cfg.at("opacity_cold_dense")},
  {cfg.at("opacity_hydrogen_response"),0,.16,1,3.5,6.6});
 TabulatedConduction conduction(cfg.at("conduction"));
 auto domain=ConductiveInteriorOpacity::Domain{};domain.minimum_T=std::pow(10.,5.7);
 ConductiveInteriorOpacity interior(raw,conduction,.001,1,domain);
 ConductiveInteriorOpacity radiation(interior,conduction,.003,1,ConductiveInteriorOpacity::ionized_hydrogen_envelope_domain());
 ScreenedCollisionTransport collisions(cfg.at("collisions"));
 ScreenedMetalMicroscopicTransport base(eos,collisions,true,5e4,{true,true,true},true);
 base.use_phase_mobility(o,.01);
 driver::ConvectiveMaterialHeat material_heat(eos,conduction);
 EnvelopeTransport envelope_heat(material_heat,base,2e6,3e6);
 BrownFingeringTransport transport(eos,radiation,envelope_heat,collisions,OscillatoryMixing::growth_squared);
 set_composition_buoyancy_reuse(0);
 std::ifstream fixture(std::string(EMBER_TEST_DATA_DIR)+"/cold_fingering_faces.dat");
 unsigned faces=0,entries=0,failures=0;double worst=0;
 std::size_t index;double ml,mh;
 while(fixture>>index>>ml>>mh) {
  Point lo,hi;Composition a,b;
  for(int side=0;side<2;++side) {
   auto& y=side?hi:lo;auto& c=side?b:a;
   fixture>>y.lnr>>y.lnrho>>y.lnT>>y.L;
   for(auto& x:c.X)fixture>>x;
   CNAbundances cn;for(auto& x:cn)fixture>>x;
   c.cn_molality=cn;c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
   c.cn_mass_convention=CNMassConvention::explicit_metal_mass;
  }
  if(!fixture)throw std::runtime_error("incomplete fingering fixture");
  ++faces;
  MetalCNFaceResponse derivative;
  transport.add_mixing_flux(derivative,index,ml,mh,lo,a,hi,b,true);
  const auto f=transport.face(index,ml,mh,lo,a,hi,b);
  if((index==408)!=(f.mass_conductance==0))throw std::runtime_error("fixture mixing classification changed");
  if(index==408)continue;
  // These physical faces reproduce a stalled species solve at trace He4.
  // A larger independent perturbation checks the derivative against changes
  // in the actual flux, including the composition-dependent mixing response.
  // Both scales remain small compared with the available He4 fraction.
  for(int side=0;side<2;++side)for(std::size_t col=0;col<2;++col) {
   const auto& c=side?b:a;const auto x=metal_cn_abundances(c);
   for(double fraction:{.003,.01}) {
    const double h=std::min(3e-6*x[col],fraction*c.X[2]);
    auto xp=x,xm=x;xp[col]+=h;xm[col]-=h;
    const auto cp=metal_cn_composition(c,xp),cm=metal_cn_composition(c,xm);
    MetalCNFaceResponse plus,minus;
    transport.add_mixing_flux(plus,index,ml,mh,lo,side?a:cp,hi,side?cp:b,false);
    transport.add_mixing_flux(minus,index,ml,mh,lo,side?a:cm,hi,side?cm:b,false);
    for(std::size_t row=0;row<2;++row) {
     const double actual=(side?derivative.dright:derivative.dleft)[row][col];
     const double expected=(plus.rate[row]-minus.rate[row])/(xp[col]-xm[col]);
     const double error=std::abs(actual-expected)/std::max(std::abs(expected),1.);
     ++entries;worst=std::max(worst,error);
     if(error>.05 || actual*expected<0) {
      ++failures;std::cerr<<"face "<<index<<" side "<<side<<" row "<<row<<" col "<<col
        <<" step fraction "<<fraction<<" derivative "<<actual<<" expected "<<expected<<" relative error "<<error<<'\n';
     }
    }
   }
  }
 }
 if(faces!=4 || entries!=48)throw std::runtime_error("missing fingering regression faces");
 std::cout<<entries<<" flux derivative checks, maximum relative difference "<<worst<<", "<<failures<<" failures\n";
 return failures?1:0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}

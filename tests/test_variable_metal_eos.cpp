#include "ember/eos_variable_metal.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <optional>

using namespace ember;
namespace {
int failures=0;
void check(bool ok,const char* name,double value=0) {
  failures+=!ok;std::printf("[%s] %s (%.10g)\n",ok?"PASS":"FAIL",name,value);
}
template<class F>bool throws(F f){try{f();}catch(const std::exception&){return true;}return false;}
Composition composition(double H,double Y3,double Z) {
  auto c=solar_scaled(H,Z);c.basis=AbundanceBasis::baryon_mass;
  c.metal_inventory=MetalInventory::gs98;c.X[1]=Y3;c.X[2]-=Y3;return c;
}
double coefficient(double u,double v,double Z,bool cubic) {
  return constants::R_gas*(1+u+.3*v+(cubic?.7*Z+.4*Z*Z+.2*Z*Z*Z:std::exp(2*Z)));
}
// A thermodynamically regular manufactured potential with known pressure
// and energy. Non-polynomial Z dependence exercises joins independently of
// the cubic-reproduction check. The final source plane can be fully masked.
void fixture(const std::filesystem::path& dir,bool cubic,bool masked,std::size_t first_masked=4) {
  std::filesystem::create_directories(dir);
  const std::array<double,6> zs{0,.005,.02,.04,.16,.3};
  const std::array<double,4> us{0,.25,.7,1};
  const std::array<double,3> vs{0,.2,1};
  for(std::size_t iz=0;iz<zs.size();++iz)for(std::size_t iu=0;iu<us.size();++iu)
    for(std::size_t iv=0;iv<vs.size();++iv) {
      std::ofstream f(dir/(std::to_string(iz)+"-"+std::to_string(iu)+"-"+std::to_string(iv)+".dat"));
      f<<std::setprecision(17)<<"EMBER_HELMHOLTZ 2\nsource \"analytic test\"\n"
       <<"composition_proxy \"none\"\nbasis baryon_mass\nmetal_inventory gs98\ncomposition";
      auto c=composition((1-zs[iz])*us[iu],(1-zs[iz])*(1-us[iu])*vs[iv],zs[iz]);
      c.X[2]=(1-zs[iz])*(1-us[iu])*(1-vs[iv]);
      for(std::size_t k=0;k<TABLE_NSPEC;++k)f<<' '<<c.X[k];
      f<<"\nlog_t 4 4 5 6 7\nlog_q 4 -6 -2 2 6\ndata\n";
      const double a=coefficient(us[iu],vs[iv],zs[iz],cubic);
      for(int it=0;it<4;++it)for(double lq:{-6.,-2.,2.,6.}) {
        const bool valid=!(masked && iz>=first_masked);
        f<<(valid?1:0)<<' '<<a*lq*std::log(10.)<<' '<<a<<" 0 0 0 0 0 0 0\n";
      }
    }
  for(int n:{4,5,6}) {
    std::ofstream f(dir/("family"+std::to_string(n)+".dat"));
    f<<std::setprecision(17)<<"EMBER_VARIABLE_METAL_HELMHOLTZ "<<(n==4?1:2)<<"\nmetals "<<n;
    for(int i=0;i<n;++i)f<<' '<<zs[i];
    f<<"\nhydrogen_share 4 0 0.25 0.7 1\nhelium3_share 3 0 0.2 1\n";
    for(int iz=0;iz<n;++iz)for(int iu=0;iu<4;++iu)for(int iv=0;iv<3;++iv)
      f<<'"'<<iz<<'-'<<iu<<'-'<<iv<<".dat\"\n";
  }
}
double response_difference(const VariableMetalHelmholtzEos& a,const VariableMetalHelmholtzEos& b,
                           double T,double rho,const Composition& c) {
  const auto ex=a.eval_with_derivatives(T,rho,c),ey=b.eval_with_derivatives(T,rho,c);
  const auto& x=ex.state;const auto& y=ey.state;
  const auto p=a.composition_potential(T,rho,c),q=b.composition_potential(T,rho,c);
  const auto h=a.composition_heat(T,rho,c),j=b.composition_heat(T,rho,c);
  double error=0;
  auto add=[&](double r,double s){error=std::max(error,std::abs(r-s)/std::max(1.,std::abs(r)));};
  add(x.P,y.P);add(x.E,y.E);add(x.S,y.S);add(x.cp,y.cp);add(x.grad_ad,y.grad_ad);
  add(x.cv,y.cv);add(x.chiT,y.chiT);add(x.chiRho,y.chiRho);add(x.delta,y.delta);add(x.Gamma1,y.Gamma1);
  add(ex.dE_dlnRho,ey.dE_dlnRho);add(ex.dcp_dlnT,ey.dcp_dlnT);add(ex.dcp_dlnRho,ey.dcp_dlnRho);
  add(ex.ddelta_dlnT,ey.ddelta_dlnT);add(ex.ddelta_dlnRho,ey.ddelta_dlnRho);
  add(ex.dgrad_ad_dlnT,ey.dgrad_ad_dlnT);add(ex.dgrad_ad_dlnRho,ey.dgrad_ad_dlnRho);
  add(p.phi,q.phi);
  for(std::size_t k=0;k<3;++k) {
    add(p.gradient[k],q.gradient[k]);add(h.exchange_enthalpy[k],j.exchange_enthalpy[k]);
    for(std::size_t l=0;l<3;++l)add(p.hessian[k][l],q.hessian[k][l]);
    for(std::size_t l=0;l<5;++l)add(h.enthalpy_partials[k][l],j.enthalpy_partials[k][l]);
  }
  return error;
}
}
int main() {
  const auto dir=std::filesystem::temp_directory_path()/
    ("ember-metal-extension-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  try {
    fixture(dir/"cubic",true,false);fixture(dir/"smooth",false,false);fixture(dir/"masked",false,true);
    using M=HelmholtzTableEos::Mixture;
    VariableMetalHelmholtzEos cubic(dir/"cubic/family6.dat",M::allow_documented_proxy);
    VariableMetalHelmholtzEos old(dir/"smooth/family4.dat",M::allow_documented_proxy);
    VariableMetalHelmholtzEos one(dir/"smooth/family5.dat",M::allow_documented_proxy);
    VariableMetalHelmholtzEos two(dir/"smooth/family6.dat",M::allow_documented_proxy);
    VariableMetalHelmholtzEos masked(dir/"masked/family6.dat",M::allow_documented_proxy);
    VariableMetalHelmholtzEos::pack_binary(dir/"smooth/family6.dat",dir/"smooth.bin");
    VariableMetalHelmholtzEos::pack_binary(dir/"smooth/family4.dat",dir/"old.bin");
    VariableMetalHelmholtzEos::pack_binary(dir/"masked/family6.dat",dir/"masked.bin");
    VariableMetalHelmholtzEos cached(dir/"smooth.bin",M::allow_documented_proxy);
    VariableMetalHelmholtzEos cached_old(dir/"old.bin",M::allow_documented_proxy);
    VariableMetalHelmholtzEos cached_masked(dir/"masked.bin",M::allow_documented_proxy);
    const double T=3e5,rho=.3;
    double exact=0,preserved=0,extended_preserved=0,cache_error=0;
    bool thermal_channels_identical=true,unrequested_channels_nan=true;
    for(double z:{1e-6,.004,.015,.02,std::nextafter(.02,1.),.03,.04,.040001,.07,.12,.16,.18,.24,.3})
      for(double u:{.53,.73})for(double v:{.07,.15}) {
        const auto c=composition((1-z)*u,(1-z)*(1-u)*v,z);
        const auto s=cubic.eval(T,rho,c);const double a=coefficient(u,v,z,true);
        const double pr=constants::a_rad*std::pow(T,4)/3;
        exact=std::max({exact,std::abs(s.P/(rho*T*a+pr)-1),std::abs(s.E/(1.5*T*a+3*pr/rho)-1)});
        if(z<=.04) {
          preserved=std::max(preserved,response_difference(old,two,T,rho,c));
          preserved=std::max(preserved,response_difference(old,masked,T,rho,c));
          cache_error=std::max({cache_error,response_difference(old,cached_old,T,rho,c),
                               response_difference(masked,cached_masked,T,rho,c)});
        }
        cache_error=std::max(cache_error,response_difference(two,cached,T,rho,c));
        if(z<=.16)extended_preserved=std::max(extended_preserved,response_difference(one,two,T,rho,c));
        const auto full_p=two.composition_potential(T,rho,c);
        const auto thermal_p=two.composition_potential(T,rho,c,{true,true,true},false);
        const auto full_h=two.composition_heat(T,rho,c);
        const auto thermal_h=two.composition_heat(T,rho,c,{true,true,true},true,false);
        thermal_channels_identical &= full_p.phi==thermal_p.phi
          && full_p.gradient==thermal_p.gradient && full_p.dgradient_dlnT==thermal_p.dgradient_dlnT
          && full_p.dgradient_dlnRho==thermal_p.dgradient_dlnRho
          && full_h.material_delta==thermal_h.material_delta
          && full_h.exchange_enthalpy==thermal_h.exchange_enthalpy
          && full_h.radiation_enthalpy==thermal_h.radiation_enthalpy;
        for(std::size_t k=0;k<3;++k)for(std::size_t l=0;l<2;++l)
          thermal_channels_identical &= full_h.enthalpy_partials[k][l]==thermal_h.enthalpy_partials[k][l]
            && full_h.radiation_enthalpy_partials[k][l]==thermal_h.radiation_enthalpy_partials[k][l]
            && full_h.delta_partials[l]==thermal_h.delta_partials[l];
        for(std::size_t k=0;k<3;++k)for(std::size_t l=0;l<3;++l)
          unrequested_channels_nan &= std::isnan(thermal_p.hessian[k][l])
            && std::isnan(thermal_h.enthalpy_partials[k][2+l])
            && std::isnan(thermal_h.radiation_enthalpy_partials[k][2+l])
            && std::isnan(thermal_h.delta_partials[2+l]);
      }
    check(exact<2e-11,"arbitrary cubic source pressure/energy reproduced across both extensions",exact);
    check(preserved==0,"all low-Z responses unchanged, including ULP perturbation and masked distant planes",preserved);
    check(extended_preserved==0,"a further source plane preserves every previously covered response",extended_preserved);
    check(cache_error==0,"binary families preserve thermal and composition responses exactly",cache_error);
    check(thermal_channels_identical,"omitting unused composition Hessians leaves all requested channels identical");
    check(unrequested_channels_nan,"unrequested composition derivatives remain explicitly unavailable");
    {
      const auto c=composition(.3,.03,.1);
      std::optional<VariableMetalHelmholtzEos> replacement;
      replacement.emplace(dir/"smooth.bin",M::allow_documented_proxy);
      const void* address=&*replacement;
      const auto before=replacement->eval(T,rho,c);
      // Interleave derivative requests, then revisit the same state. A value
      // query must not supply zero-filled unrequested derivatives to a later
      // composition-force or heat query.
      const auto potential=replacement->composition_potential(T,rho,c);
      const auto heat=replacement->composition_heat(T,rho,c);
      replacement->composition_potential(T,rho,c,{true,true,true},false);
      const auto repeated=replacement->composition_potential(T,rho,c);
      const auto repeated_heat=replacement->composition_heat(T,rho,c);
      check(potential.gradient==repeated.gradient && potential.hessian==repeated.hessian
          && heat.enthalpy_partials==repeated_heat.enthalpy_partials,
            "repeated state retains the requested force and heat derivatives");
      check(throws([&]{replacement->eval(0,rho,c);})
          && throws([&]{replacement->eval(T,0,c);}),
            "a supported cached state does not admit invalid temperature or density");
      auto wrong_basis=c;wrong_basis.basis=static_cast<AbundanceBasis>(255);
      check(throws([&]{replacement->eval(T,rho,wrong_basis);}),
            "cached state still requires the correct abundance basis");
      replacement.reset();
      replacement.emplace(dir/"masked.bin",M::allow_documented_proxy);
      check(address==&*replacement,"replacement EOS occupies the same address");
      check(throws([&]{replacement->eval(T,rho,c);}),
            "replacement EOS retains its own source masks after a cached query");
      replacement.reset();
      replacement.emplace(dir/"cubic/family6.dat",M::allow_documented_proxy);
      const auto after=replacement->eval(T,rho,c),expected=cubic.eval(T,rho,c);
      check(after.P==expected.P && after.E==expected.E && after.P!=before.P,
            "replacement EOS uses the new material at the same composition");
    }
    check(!throws([&]{cached.validate_composition_domain(T,rho,composition(.3,.03,.1));}),
          "response reuse allows a supported composition state");
    check(throws([&]{cached_masked.validate_composition_domain(T,rho,composition(.3,.03,.1));}),
          "response reuse rejects masked composition derivative planes");
    check(throws([&]{cached.validate_composition_domain(T,rho,composition(.3,.03,.301));})
        && throws([&]{cached.validate_composition_domain(1.,rho,composition(.3,.03,.1));}),
          "response reuse retains composition and temperature source limits");
    check(throws([&]{cached_masked.eval(T,rho,composition(.3,.03,.1));}),"binary family preserves source masks");
    check(throws([&]{cached.eval(T,rho,composition(.3,.03,.301));}),"binary family still rejects extrapolation");
    check(throws([&]{VariableMetalHelmholtzEos::pack_binary(dir/"smooth/family6.dat",dir/"smooth.bin");}),
          "packing refuses to overwrite an existing input");
    std::filesystem::copy_file(dir/"smooth.bin",dir/"truncated.bin");
    std::filesystem::resize_file(dir/"truncated.bin",std::filesystem::file_size(dir/"truncated.bin")-1);
    check(throws([&]{VariableMetalHelmholtzEos bad(dir/"truncated.bin",M::allow_documented_proxy);}),
          "truncated binary input is rejected");
    std::filesystem::copy_file(dir/"smooth.bin",dir/"corrupt.bin");
    {
      std::fstream file(dir/"corrupt.bin",std::ios::in|std::ios::out|std::ios::binary);
      file.seekp(-80,std::ios::end);auto value=std::bit_cast<std::uint64_t>(2.);
      if constexpr(std::endian::native!=std::endian::little)value=std::byteswap(value);
      file.write(reinterpret_cast<const char*>(&value),8);
    }
    check(throws([&]{VariableMetalHelmholtzEos bad(dir/"corrupt.bin",M::allow_documented_proxy);}),
          "invalid binary mask is rejected");
    check(throws([&]{masked.eval(T,rho,composition(.3,.03,.1));}),"needed masked high-Z source is still rejected");
    check(throws([&]{two.eval(T,rho,composition(.3,.03,.301));}),"metal extrapolation remains rejected");
    double continuity=0,derivative=0,firstlaw=0;
    for(double z:{.04,.16}) {
      const auto a=two.composition_potential(T,rho,composition(.3,.03,z-1e-10));
      const auto b=two.composition_potential(T,rho,composition(.3,.03,z+1e-10));
      for(std::size_t k=0;k<3;++k) {
        continuity=std::max(continuity,std::abs(a.gradient[k]-b.gradient[k])/constants::R_gas);
        for(std::size_t l=0;l<3;++l)
          continuity=std::max(continuity,std::abs(a.hessian[k][l]-b.hessian[k][l])/constants::R_gas);
      }
    }
    check(continuity<2e-6,"potential gradient and Hessian continuous at both new joins",continuity);
    for(double z:{.039,.041,.1,.159,.161,.23}) {
      const auto c=composition(.3,.03,z);const auto s=two.eval(T,rho,c);
      const auto p=two.composition_potential(T,rho,c);const auto heat=two.composition_heat(T,rho,c);
      const double step=2e-6;
      for(std::size_t k=0;k<3;++k) {
        std::array<double,3> coords{.3,.03,z};coords[k]+=step;
        const auto plus=composition(coords[0],coords[1],coords[2]);coords[k]-=2*step;
        const auto minus=composition(coords[0],coords[1],coords[2]);
        const auto a=two.composition_potential(T,rho,plus),b=two.composition_potential(T,rho,minus);
        derivative=std::max(derivative,std::abs((a.phi-b.phi)/(2*step)-p.gradient[k])/constants::R_gas);
        for(std::size_t l=0;l<3;++l)
          derivative=std::max(derivative,std::abs((a.gradient[l]-b.gradient[l])/(2*step)-p.hessian[l][k])/constants::R_gas);
        const auto ha=two.composition_heat(T,rho,plus),hb=two.composition_heat(T,rho,minus);
        for(std::size_t l=0;l<3;++l)
          derivative=std::max(derivative,std::abs((ha.exchange_enthalpy[l]-hb.exchange_enthalpy[l])/(2*step)
                                                -heat.enthalpy_partials[l][2+k])/(T*constants::R_gas));
      }
      const auto a=two.eval(T,rho*std::exp(step),c),b=two.eval(T,rho*std::exp(-step),c);
      firstlaw=std::max(firstlaw,std::abs((a.E-b.E)/(2*step)-s.P/rho*(1-s.chiT))/(T*s.cv));
    }
    check(derivative<5e-6,"composition forces, Hessian and material-heat derivatives match finite differences",derivative);
    check(firstlaw<2e-7,"enriched states retain thermodynamic first-law identity",firstlaw);
    {
      using L=VariableMetalHelmholtzEos::LowMetalInterpolation;
      fixture(dir/"trace_masked",false,true,3);
      VariableMetalHelmholtzEos local(dir/"smooth/family6.dat",M::allow_documented_proxy,L::quadratic);
      VariableMetalHelmholtzEos trace_masked(dir/"trace_masked/family6.dat",M::allow_documented_proxy,L::quadratic);
      VariableMetalHelmholtzEos original_masked(dir/"trace_masked/family6.dat",M::allow_documented_proxy);
      const auto trace=composition(.993,.003,1e-30);
      check(!throws([&]{trace_masked.eval(T,rho,trace);trace_masked.composition_potential(T,rho,trace);
            trace_masked.composition_heat(T,rho,trace);trace_masked.validate_composition_domain(T,rho,trace);}),
            "trace-metal values, forces, heat and reuse need only the three low-metal planes");
      check(throws([&]{original_masked.eval(T,rho,trace);}),
            "default cubic retains the original source requirements");
      check(throws([&]{trace_masked.eval(T,rho,composition(.3,.03,.01));}),
            "blended interval still requires its fourth source plane");
      double joins=0,fd=0,preserved_high=0;
      for(double z:{.005,.02}) {
        const auto a=local.composition_potential(T,rho,composition(.3,.03,z-1e-10));
        const auto b=local.composition_potential(T,rho,composition(.3,.03,z+1e-10));
        for(std::size_t k=0;k<3;++k) {
          joins=std::max(joins,std::abs(a.gradient[k]-b.gradient[k])/constants::R_gas);
          for(std::size_t l=0;l<3;++l)
            joins=std::max(joins,std::abs(a.hessian[k][l]-b.hessian[k][l])/constants::R_gas);
        }
      }
      for(double z:{.001,.004,.005,.0075,.013,.019,.02,.025}) {
        const auto c=composition(.3,.03,z);
        const auto p=local.composition_potential(T,rho,c);
        const auto h=local.composition_heat(T,rho,c);
        const double step=5e-7;
        for(std::size_t k=0;k<3;++k) {
          std::array<double,3> x{.3,.03,z};x[k]+=step;const auto plus=composition(x[0],x[1],x[2]);
          x[k]-=2*step;const auto minus=composition(x[0],x[1],x[2]);
          const auto a=local.composition_potential(T,rho,plus),b=local.composition_potential(T,rho,minus);
          const auto ha=local.composition_heat(T,rho,plus),hb=local.composition_heat(T,rho,minus);
          fd=std::max(fd,std::abs((a.phi-b.phi)/(2*step)-p.gradient[k])/constants::R_gas);
          for(std::size_t l=0;l<3;++l) {
            fd=std::max(fd,std::abs((a.gradient[l]-b.gradient[l])/(2*step)-p.hessian[l][k])/constants::R_gas);
            fd=std::max(fd,std::abs((ha.exchange_enthalpy[l]-hb.exchange_enthalpy[l])/(2*step)
                                     -h.enthalpy_partials[l][2+k])/(T*constants::R_gas));
          }
        }
        if(c.Z()>=.02)preserved_high=std::max(preserved_high,response_difference(two,local,T,rho,c));
      }
      check(joins<2e-6,"low-metal joins preserve continuous force and curvature",joins);
      check(fd<5e-6,"low-metal forces, curvature and heat differentiate the same potential",fd);
      check(preserved_high==0,"source interpolation at and above the second nonzero metal node is unchanged",preserved_high);
      const auto a=local.composition_potential(T,rho,composition(.993,.003,1e-30));
      const auto b=local.composition_potential(T,rho,composition(.993,.003,1e-29));
      double ions=0;for(const auto& m:gs98_metals)if(m.charge!=19)ions+=m.fraction/m.mass_number;
      check(std::abs((b.gradient[2]-a.gradient[2])/(constants::R_gas*ions*std::log(10.))-1)<1e-10,
            "trace metal chemical potential retains exact logarithmic ideal mixing");
    }
  }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());++failures;}
  std::filesystem::remove_all(dir);
  return failures?1:0;
}

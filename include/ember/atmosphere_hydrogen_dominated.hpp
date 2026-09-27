#pragma once
#include "ember/atmosphere.hpp"
#include "ember/constants.hpp"
#include "ember/interp.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
namespace ember {
// Explicit warm-atmosphere approximation. Only the matching gas T and P
// neglect trace helium. Every density inversion uses the actual composition.
class HydrogenDominatedAtmosphereGrid final: public Atmosphere {
public:
 enum class Approximation { neglect_trace_atmospheric_helium };
 HydrogenDominatedAtmosphereGrid(const Eos& eos,const std::filesystem::path& path,
     Approximation approximation,double maximum_helium):eos_(eos) {
  std::ifstream in(path);if(!in)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: cannot open source");
  read(in,approximation,maximum_helium);
 }
 HydrogenDominatedAtmosphereGrid(const Eos& eos,std::istream& in,
     Approximation approximation,double maximum_helium):eos_(eos){read(in,approximation,maximum_helium);}
 bool covers(double t,double g,const Composition& c)const{
  if(!positive(t)||!positive(g)||c.basis!=AbundanceBasis::baryon_mass||
     c.metal_inventory!=MetalInventory::gs98||!std::isfinite(c.sum())||std::abs(c.sum()-1)>1e-10)return false;
  for(double x:c.X)if(!std::isfinite(x)||x<0||x>1)return false;
  if(c[Species::H2]!=0 || c.X[1]+c.X[2]>maximum_helium_+8*std::numeric_limits<double>::epsilon())return false;
  const double z=c.Z();
  for(std::size_t k=0;k<pattern_.size();++k)
   if(std::abs(c.X[k+METAL_BEGIN]-z*pattern_[k])>1e-12*std::max(z,1e-30))return false;
  const double eps=8*std::numeric_limits<double>::epsilon()*axes_[0].back();
  if(z<axes_[0].front()-eps||z>axes_[0].back()+eps)return false;
  auto q=coordinates(t,g,z);
  for(std::size_t k=0;k<3;++k)if(q[k]<axes_[k].front()||q[k]>axes_[k].back())return false;
  return true;
 }
 AtmosphereState eval(double t,double g,const Composition& c)const override {
  if(!covers(t,g,c))throw std::domain_error("HydrogenDominatedAtmosphereGrid: outside explicit composition or source domain");
  auto q=coordinates(t,g,c.Z());auto T=sample(logT_,q),pg=sample(logPg_,q);
  const double pr=constants::a_rad*std::pow(T[0],4)/3;
  AtmosphereState s{};s.T=T[0];s.Pgas=pg[0];s.P=s.Pgas+pr;s.tau=tau_;
  if(!positive(s.P)||s.P==pr)throw std::domain_error("HydrogenDominatedAtmosphereGrid: invalid pressure");
  s.rho=eos_.rho_from_PT(s.T,s.P,c);
  s.dlnT_dlnTeff=T[2];s.dlnT_dlng=T[3];
  s.dlnP_dlnTeff=(s.Pgas*pg[2]+4*pr*T[2])/s.P;
  s.dlnP_dlng=(s.Pgas*pg[3]+4*pr*T[3])/s.P;return s;
 }
 double tau_match()const{return tau_;}
 const char* name()const override{return source_.c_str();}
private:
 static bool positive(double x){return std::isfinite(x)&&x>0;}
 std::array<double,3> coordinates(double t,double g,double z)const {
  // Only an endpoint sum's floating-point roundoff can reach this clamp:
  // covers() first rejects any actual metallicity outside the source grid.
  return {std::clamp(z,axes_[0].front(),axes_[0].back()),std::log10(t),std::log10(g)};
 }
 void read(std::istream& in,Approximation a,double maximum){
  auto label=[&](const char* expected){std::string s;if(!(in>>s)||s!=expected)throw std::runtime_error(std::string("HydrogenDominatedAtmosphereGrid: expected ")+expected);};
  label("EMBER_HYDROGEN_DOMINATED_ATMOSPHERE");int version{};in>>version;if(version!=1)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: version");
  label("source");in>>std::quoted(source_);label("approximation");in>>std::quoted(approximation_);
  if(!in||source_.empty()||approximation_.empty()||a!=Approximation::neglect_trace_atmospheric_helium||
     !positive(maximum)||maximum>=1)throw std::invalid_argument("HydrogenDominatedAtmosphereGrid: explicit approximation required");
  label("basis");label("baryon_mass");label("tau");in>>tau_;
  label("maximum_helium");in>>maximum_helium_;label("source_helium");in>>source_helium_;
  if(!in||!positive(tau_)||!positive(maximum_helium_)||maximum_helium_>maximum||
     (!std::isfinite(source_helium_)||source_helium_<0)||source_helium_>maximum_helium_)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: approximation bounds");
  label("metal_pattern");for(double& x:pattern_){in>>x;if(!in||!std::isfinite(x)||x<0||x>1)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: metal pattern");}
  if(std::abs(std::accumulate(pattern_.begin(),pattern_.end(),0.)-1)>1e-12)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: metal normalization");
  std::size_t count=1;const char* labels[]={"metallicity","log_teff","log_g"};
  for(std::size_t k=0;k<3;++k){
   label(labels[k]);std::size_t n{};in>>n;if(!in||n<2||n>10000||count>1000000/n)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: size");
   count*=n;axes_[k].resize(n);
   for(std::size_t j=0;j<n;++j){double& v=axes_[k][j];in>>v;if(!in||!std::isfinite(v)||(j&&v<=axes_[k][j-1])||
      (k==0?(v<0||v>=1):!positive(std::pow(10.,v))))throw std::runtime_error("HydrogenDominatedAtmosphereGrid: axis");}
  }
  label("data");logT_.resize(count);logPg_.resize(count);
  for(std::size_t j=0;j<count;++j){in>>logT_[j]>>logPg_[j];if(!in||!positive(std::pow(10.,logT_[j]))||!positive(std::pow(10.,logPg_[j])))throw std::runtime_error("HydrogenDominatedAtmosphereGrid: incomplete source");}
  std::string extra;if(in>>extra)throw std::runtime_error("HydrogenDominatedAtmosphereGrid: trailing data");
 }
 std::array<double,4> sample(const std::vector<double>& f,const std::array<double,3>& q)const {
  std::array<std::size_t,3> lo{};std::array<double,3> u{},width{};
  for(std::size_t k=0;k<3;++k){lo[k]=interp::locate(axes_[k],q[k]);width[k]=axes_[k][lo[k]+1]-axes_[k][lo[k]];u[k]=(q[k]-axes_[k][lo[k]])/width[k];}
  std::array<double,4> out{};
  for(unsigned corner=0;corner<8;++corner){
   std::size_t idx=0;std::array<double,3>w{};
   for(unsigned k=0;k<3;++k){bool high=corner&(1U<<k);idx=idx*axes_[k].size()+lo[k]+high;w[k]=high?u[k]:1-u[k];}
   out[0]+=f[idx]*w[0]*w[1]*w[2];
   for(unsigned k=0;k<3;++k){double d=(corner&(1U<<k)?1.:-1.)/width[k];for(unsigned j=0;j<3;++j)if(j!=k)d*=w[j];out[k+1]+=f[idx]*d;}
  }
  out[0]=std::pow(10.,out[0]);out[1]*=std::log(10.);return out;
 }
 const Eos& eos_;std::string source_,approximation_;double tau_{},maximum_helium_{},source_helium_{};
 std::array<double,NMETALS> pattern_{};
 std::array<std::vector<double>,3> axes_;
 std::vector<double> logT_,logPg_;
};

} // namespace ember

#include "ember/envelope_gibbs.hpp"
#include "ember/constants.hpp"
#include "jet2.hpp"
#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <iomanip>

namespace ember {
namespace {
using J=detail::Jet2;
void require(bool ok,const char* why) {
  if(!ok)throw std::domain_error(std::string("Gibbs envelope: ")+why);
}
class Potential {
  std::vector<double> t_,p_,coeff_,source_t_,source_p_;
  std::vector<int> support_;
  std::string source_;
  int degree_{};
  std::size_t nt_{},np_{};
  static std::pair<std::size_t,std::vector<J>> basis(const std::vector<double>& axis,int k,const J& x) {
    const std::size_t n=axis.size()-k-1;
    require(x.value>=axis[k] && x.value<=axis[n],"outside potential axes");
    const auto span=std::min(static_cast<std::size_t>(std::upper_bound(axis.begin(),axis.end(),x.value)-axis.begin()-1),n-1);
    std::vector<J> values(k+1),left(k+1),right(k+1);values[0]=1.;
    for(int j=1;j<=k;++j) {
      left[j]=x-axis[span+1-j];right[j]=axis[span+j]-x;J saved;
      for(int r=0;r<j;++r) {
        const J term=values[r]/(right[r+1]+left[j-r]);
        values[r]=saved+right[r+1]*term;saved=left[j-r]*term;
      }
      values[j]=saved;
    }
    return {span-k,values};
  }
public:
  explicit Potential(const std::string& path):source_(path) {
    std::ifstream in(path);std::string tag;in>>tag>>degree_;
    require(tag=="EMBER_GIBBS_GAS_V1" && degree_>=3 && degree_<=5,"invalid potential format");
    auto read=[&](std::vector<double>& a) {
      std::size_t n{};in>>n;require(in.good() && n>=2 && n<=1000000,"invalid array size");
      a.resize(n);for(auto& v:a){in>>v;require(in.good() && std::isfinite(v),"invalid table value");}
    };
    read(t_);read(p_);read(coeff_);read(source_t_);read(source_p_);
    for(const auto* axis:{&t_,&p_}) {
      require(axis->size()>=2*static_cast<std::size_t>(degree_+1),"short knot axis");
      const auto n=axis->size()-degree_-1;
      require((*axis)[degree_]<(*axis)[n] && std::is_sorted(axis->begin(),axis->end()),"unordered knot axis");
      for(int k=0;k<=degree_;++k)
        require((*axis)[k]==(*axis)[degree_] && (*axis)[n+k]==(*axis)[n],"unclamped knot axis");
      for(std::size_t k=degree_+1;k<=n;++k)
        require((*axis)[k]>(*axis)[k-1],"repeated interior knot");
    }
    nt_=t_.size()-degree_-1;np_=p_.size()-degree_-1;
    require(coeff_.size()==nt_*np_,"coefficient count does not match axes");
    for(const auto* axis:{&source_t_,&source_p_})
      for(std::size_t i=1;i<axis->size();++i)require((*axis)[i]>(*axis)[i-1],"unordered source axis");
    require(source_t_.size()<=10000000/source_p_.size(),"support mask too large");
    support_.resize(source_t_.size()*source_p_.size());
    for(auto& v:support_){in>>v;require(in.good() && (v==0 || v==1),"invalid source mask");}
    for(std::size_t i=0;i<source_t_.size();++i)for(std::size_t j=0;j<source_p_.size();++j)
      if(support_[i*source_p_.size()+j])
        require(source_t_[i]>=t_[degree_] && source_t_[i]<=t_[nt_]
             && source_p_[j]>=p_[degree_] && source_p_[j]<=p_[np_],"supported source node exceeds potential");
    std::string extra;require(!(in>>extra),"unexpected trailing data");
  }
  double minimum_temperature() const { return std::max(source_t_.front(),t_[degree_]); }
  bool supported(double t,double p) const {
    if(!std::isfinite(t+p) || t<source_t_.front() || t>source_t_.back()
       || p<source_p_.front() || p>source_p_.back())return false;
    auto cell=[](const std::vector<double>& a,double x) {
      return std::min(static_cast<std::size_t>(std::upper_bound(a.begin(),a.end(),x)-a.begin()-1),a.size()-2);
    };
    const auto i=cell(source_t_,t),j=cell(source_p_,p);
    for(int a=0;a<2;++a)for(int b=0;b<2;++b)
      if(!support_[(i+a)*source_p_.size()+j+b])return false;
    return true;
  }
  J value(const J& t,const J& p) const {
    if(!supported(t.value,p.value)) {
      std::ostringstream message;message<<std::setprecision(4)
        <<"Gibbs envelope: state outside supported source; T="<<std::exp(t.value)
        <<" K, gas P="<<std::exp(p.value)<<" dyn/cm2, source="<<source_;
      throw std::domain_error(message.str());
    }
    const auto [it,bt]=basis(t_,degree_,t);const auto [ip,bp]=basis(p_,degree_,p);J f;
    for(int a=0;a<=degree_;++a)for(int b=0;b<=degree_;++b)
      f=f+coeff_[(it+a)*np_+ip+b]*bt[a]*bp[b];
    return f;
  }
  void check_join(const Potential& cold) const {
    const double t=minimum_temperature();std::size_t checked=0;
    for(std::size_t i=1;i<source_p_.size();++i) {
      const double p=.5*(source_p_[i-1]+source_p_[i]);
      if(!supported(t,p) || !cold.supported(t,p))continue;
      const J a=value(J::variable(t,0),J::variable(p,1));
      const J b=cold.value(J::variable(t,0),J::variable(p,1));
      auto same=[](double x,double y){return std::abs(x-y)<=1e-7*(1+std::max(std::abs(x),std::abs(y)));};
      require(same(a.value,b.value),"discontinuous potential join");
      for(int j=0;j<2;++j) {
        require(same(a.d[j],b.d[j]),"discontinuous first derivative at join");
        for(int k=0;k<2;++k)require(same(a.h[j][k],b.h[j][k]),"discontinuous response at join");
      }
      ++checked;
    }
    require(checked>0,"potential join has no common supported states");
  }
};

struct MetalGas {
  double ions{},electrons{},ion_reference{},number_entropy{},electron_reference{};
  MetalGas() {
    using namespace constants;
    const auto reference=[](double mass,double spin) {
      return 1.5*std::log(h*h/(2*M_PI*mass*kB))-std::log(spin*kB);
    };
    for(const auto& e:gs98_metals) {
      const double n=e.fraction/e.mass_number;
      ions+=n;electrons+=n*e.charge;
      ion_reference+=n*reference(e.atomic_weight/NA,1.);
      number_entropy+=n*std::log(n);
    }
    electron_reference=reference(me,2.);
  }
};
} // namespace

struct GibbsEnvelope::Impl {
  Potential h,he,hw,hew;
  EnvelopeMetals metals;
  explicit Impl(const Files& f,EnvelopeMetals m):h(f.hydrogen),he(f.helium),hw(f.hydrogen_warm),hew(f.helium_warm),metals(m) {
    hw.check_join(h);hew.check_join(he);
  }
};
GibbsEnvelope::GibbsEnvelope(const Files& f,EnvelopeMetals m):impl_(std::make_unique<Impl>(f,m)) {}
GibbsEnvelope::~GibbsEnvelope()=default;

GibbsEnvelope::Thermodynamics GibbsEnvelope::evaluate(double T,double P,const Composition& c) const {
  require(std::isfinite(T+P) && T>0 && P>0,"invalid temperature or pressure");
  require(c.basis==AbundanceBasis::baryon_mass,"baryon mass fractions required");
  require(std::isfinite(c.sum()) && std::abs(c.sum()-1)<1e-12,"composition is not normalized");
  for(double x:c.X)require(std::isfinite(x) && x>=0,"invalid mass fraction");
  const double X=c.h1(),D=c[Species::H2],Y3=c[Species::He3],Y4=c[Species::He4],Z=c.Z();
  require(D<=1e-4 && Z<=.04,"composition outside trace-D and metal approximation");
  require(impl_->metals!=EnvelopeMetals::reject || Z<=1e-12,"metal approximation must be selected explicitly");
  require(impl_->metals==EnvelopeMetals::reject || Z==0 || c.metal_inventory==MetalInventory::gs98,"GS98 metal distribution required");
  const J t=J::variable(std::log(T),0),p=J::variable(std::log(P),1);
  const J Pg=exp(p)-constants::a_rad/3.*exp(4.*t);
  require(Pg.value>0,"nonpositive gas pressure");
  const J pg=log(Pg);
  J f;
  if(X+D>0) {
    const auto& h=t.value>impl_->hw.minimum_temperature()?impl_->hw:impl_->h;
    f=(X+D/2.)*(constants::R_gas*1.00782503)*h.value(t,pg);
  }
  const double heweight=Y4+4./3.*Y3;
  if(heweight>0) {
    const auto& he=t.value>impl_->hew.minimum_temperature()?impl_->hew:impl_->he;
    f=f+heweight*(constants::R_gas*4.00260325/4.)*he.value(t,pg);
  }
  // The isotope terms alter the classical entropy reference but not heat
  // capacity. Trace D omits molecular rotational and zero-point isotope shifts.
  f=f+constants::R_gas*(.5*Y3*std::log(4.00260325/3.01602932)
                        +.75*D*std::log(1.00782503/2.01410177812));
  const double n[]={X,D/2.,Y3/3.,Y4/4.};double total=0,number_entropy=0;
  for(double v:n)if(v>0){total+=v;number_entropy+=v*std::log(v);}
  if(impl_->metals!=EnvelopeMetals::reject && Z>0) {
    static const MetalGas metal;
    f=f+constants::R_gas*Z*(metal.ions*(pg-2.5*t)+metal.ion_reference);
    total+=Z*metal.ions;
    number_entropy+=Z*(metal.ions*std::log(Z)+metal.number_entropy);
    if(impl_->metals==EnvelopeMetals::ionized) {
      const double electrons=Z*metal.electrons;
      f=f+constants::R_gas*electrons*(pg-2.5*t+metal.electron_reference);
      total+=electrons;number_entropy+=electrons*std::log(electrons);
    }
  }
  f=f-constants::R_gas*(total*std::log(total)-number_entropy);
  const double gp=f.d[1],cp=-(f.d[0]+f.h[0][0]),delta=1+f.h[0][1]/gp;
  const double chi=1/(1-f.h[1][1]/gp),cv=cp-gp*delta*delta*chi;
  const double rho=P/(T*gp),grad=gp*delta/cp;
  require(std::isfinite(rho+cp+delta+chi+cv) && rho>0 && cp>0 && chi>0 && cv>0,"nonphysical mixture response");
  const double entropy=-(f.value+f.d[0]),energy=-T*(f.d[0]+gp);
  require(std::isfinite(entropy+energy),"nonfinite entropy or energy");
  return {{rho,cp,delta,grad,chi},entropy,energy};
}
EnvelopeThermodynamics::State GibbsEnvelope::at_pressure(double t,double p,const Composition& c,double& guess) const {
  const auto result=evaluate(std::exp(t),std::exp(p),c).pressure;guess=std::log(result.rho);return result;
}
double GibbsEnvelope::rho_from_PT(double T,double P,const Composition& c,double) const {
  return evaluate(T,P,c).pressure.rho;
}
} // namespace ember

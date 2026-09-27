#include "ember/collision_transport.hpp"
#include "ember/gs98_mixture.hpp"
#include "ember/detail/differential.hpp"
#include <type_traits>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace ember {
namespace {
constexpr double kb=1.380649e-16, mu=1.66053906660e-24;
constexpr double me=9.1093837015e-28, hbar=1.054571817e-27;
constexpr double e2=2.3070775523417355e-19, pi=std::numbers::pi;
void require(bool ok,const char* why) {
  if(!ok) throw std::domain_error(std::string("collision transport: ")+why);
}
using detail::Differential;
using std::sqrt;using std::exp;using std::log;using std::log1p;using std::pow;using std::cbrt;
double value(double x) {return x;}
template<std::size_t N> double value(const Differential<N>& x) {return x.value;}
bool finite(double x) {return std::isfinite(x);}
template<std::size_t N> bool finite(const Differential<N>& x) {
  if(!std::isfinite(x.value)) return false;
  for(double d:x.d) if(!std::isfinite(d)) return false;
  return true;
}
template<class S> void positive(const S& x) {require(finite(x) && value(x)>0,"positive finite scale required");}
template<class S> S magnitude(const S& x) {return value(x)>=0 ? x : -x;}
template<class S> S positive_part(const S& x) {return value(x)>0 ? x : S(0);}
bool exactly_zero(double x) {return x==0;}
template<std::size_t N> bool exactly_zero(const Differential<N>& x) {
  if(x.value!=0)return false;
  return std::all_of(x.d.begin(),x.d.end(),[](double d) {return d==0;});
}

template<class S> struct Matrix {
  std::size_t rows{},cols{};
  std::vector<S> a;
  Matrix(std::size_t r,std::size_t c):rows(r),cols(c),a(r*c) {}
  S& operator()(std::size_t i,std::size_t j) {return a[i*cols+j];}
  S operator()(std::size_t i,std::size_t j) const {return a[i*cols+j];}
};
template<std::size_t N> Matrix<Differential<N>> compose_binary_response(
    const Matrix<Differential<2>>& local,
    const Differential<N>& x,const Differential<N>& y) {
  Matrix<Differential<N>> out(local.rows,local.cols);
  for(std::size_t i=0;i<local.a.size();++i) {
    out.a[i].value=local.a[i].value;
    for(std::size_t k=0;k<N;++k)
      out.a[i].d[k]=local.a[i].d[0]*x.d[k]+local.a[i].d[1]*y.d[k];
  }
  return out;
}
template<class S> Matrix<S> cholesky(const Matrix<S>& a) {
  Matrix<S> l(a.rows,a.cols);
  for(std::size_t i=0;i<a.rows;++i) for(std::size_t j=0;j<=i;++j) {
    S v=a(i,j);
    for(std::size_t k=0;k<j;++k) v-=l(i,k)*l(j,k);
    if(i==j) {positive(v);l(i,j)=sqrt(v);}
    else l(i,j)=v/l(j,j);
  }
  return l;
}
template<class S> Matrix<S> solve(const Matrix<S>& a,const Matrix<S>& rhs,double& error) {
  require(a.rows==a.cols && rhs.rows==a.rows,"matrix dimensions differ");
  const std::size_t n=a.rows;
  std::vector<S> s(n);
  Matrix<S> scaled(n,n),r(rhs.rows,rhs.cols),x(rhs.rows,rhs.cols);
  for(std::size_t i=0;i<n;++i) {positive(a(i,i));s[i]=sqrt(a(i,i));}
  for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<n;++j)
    scaled(i,j)=a(i,j)/s[i]/s[j];
  const auto l=cholesky(scaled);
  for(std::size_t k=0;k<rhs.cols;++k) {
    for(std::size_t i=0;i<n;++i) {
      S v=rhs(i,k)/s[i];
      for(std::size_t j=0;j<i;++j) v-=l(i,j)*r(j,k);
      r(i,k)=v/l(i,i);
    }
    for(std::size_t i=n;i-->0;) {
      S v=r(i,k);
      for(std::size_t j=i+1;j<n;++j) v-=l(j,i)*x(j,k);
      x(i,k)=v/l(i,i);
    }
  }
  double an=0,xn=0,bn=0,rn=0;
  for(std::size_t i=0;i<n;++i) {
    double ar=0,xr=0,br=0,rr=0;
    for(std::size_t j=0;j<n;++j) ar+=std::abs(value(scaled(i,j)));
    for(std::size_t k=0;k<rhs.cols;++k) {
      const double b=value(rhs(i,k)/s[i]);double v=-b;
      for(std::size_t j=0;j<n;++j) v+=value(scaled(i,j))*value(x(j,k));
      xr+=std::abs(value(x(i,k)));br+=std::abs(b);rr+=std::abs(v);
    }
    an=std::max(an,ar);xn=std::max(xn,xr);bn=std::max(bn,br);rn=std::max(rn,rr);
  }
  const double residual=(an*xn+bn)>0 ? rn/(an*xn+bn) : 0;
  require(std::isfinite(residual) && residual<1e-12,"linear solve residual failed");
  error=std::max(error,residual);
  for(std::size_t i=0;i<n;++i) for(std::size_t k=0;k<rhs.cols;++k) {
    x(i,k)/=s[i];require(finite(x(i,k)),"nonfinite linear response");
  }
  return x;
}
template<class S> struct Mixture {
  std::vector<S> x,n;
  std::vector<double> mass,z;
  std::vector<std::size_t> independent,output;
  std::size_t reference{};
  S ne{},nz2{};
  Mixture(S rho,S X,S Y3,S Z) {
    positive(rho);
    for(S v:{X,Y3,Z}) require(finite(v) && value(v)>=0,"negative mass fraction");
    const S Y4=1-X-Y3-Z;positive(Y4);
    auto add=[&](S fraction,double mass_number,double charge) {
      if(value(fraction)==0) return;
      x.push_back(fraction);mass.push_back(mass_number);z.push_back(charge);
      n.push_back(rho/mu*fraction/mass_number);
      ne+=n.back()*charge;nz2+=n.back()*charge*charge;
    };
    if(value(X)>0) {independent.push_back(x.size());output.push_back(0);add(X,1,1);}
    if(value(Y3)>0) {independent.push_back(x.size());output.push_back(1);add(Y3,3,2);}
    reference=x.size();add(Y4,4,2);
#if defined(EMBER_TWO_METAL_COLLISIONS) && EMBER_TWO_METAL_COLLISIONS
    // Collision carriers, not physical isotopes. Match number/charge moments
    // through cubic order and metal mass and mass-weighted charge. This also
    // approximates the ion contribution to static screening. The EOS and
    // nuclear inventory retain the physical mixture.
    struct Representative {double fraction,mass_number,charge;};
    constexpr std::array<Representative,2> representative{{
      {0.82531077145474385,15.160180402274264,7.5942225823740097},
      {0.17468922854525681,45.470377851069671,21.681166798218783},
    }};
    for(const auto& metal:representative) add(Z*metal.fraction,metal.mass_number,metal.charge);
#else
    for(const auto& metal:gs98_metals) add(Z*metal.fraction,metal.mass_number,metal.charge);
#endif
    positive(ne);positive(nz2);
  }
};
template<class S> S occupation(S eta,S x) {
  const S z=exp(-magnitude(x-eta));
  return value(x)>=value(eta) ? z/(1+z) : 1/(1+z);
}
template<class S> S occupation_derivative(S eta,S x) {
  const S z=exp(-magnitude(x-eta));return z/((1+z)*(1+z));
}
template<class S> S bracket(S b) {
  if(value(b)<1e-3) return b*b*(.5+b*(-2./3+b*(.75+b*(-.8+b*(5./6+b*(-6./7+b*7./8))))));
  return log1p(b)-b/(1+b);
}
template<class S> std::vector<S> spline_basis(const std::vector<double>& t,S x) {
  require(finite(x) && value(x)>=t.front() && value(x)<=t.back(),"electron pair query outside table");
  std::vector<S> b(t.size()-1);
  for(std::size_t i=0;i<b.size();++i)
    b[i]=(value(x)>=t[i] && value(x)<t[i+1]) || (value(x)==t.back() && t[i]<value(x) && t[i+1]==value(x)) ? 1 : 0;
  for(std::size_t d=1;d<=3;++d) for(std::size_t i=0;i<t.size()-d-1;++i) {
    const double low=t[i+d]-t[i],high=t[i+d+1]-t[i+1];
    b[i]=(low>0 ? (x-t[i])*b[i]/low : S(0))+(high>0 ? (t[i+d+1]-x)*b[i+1]/high : S(0));
  }
  b.resize(t.size()-4);return b;
}
struct Polynomial {
  std::vector<double> x,c;
  std::size_t order{};
  template<class S> std::array<S,4> eval(S v) const {
    require(value(v)>=x.front() && value(v)<=x.back(),"ion collision query outside table");
    const std::size_t count=x.size()-1;
    const auto upper=std::upper_bound(x.begin(),x.end(),value(v));
    const auto i=std::min(count-1,static_cast<std::size_t>(upper-x.begin()-1));
    std::array<S,4> out{};
    for(std::size_t k=0;k<order;++k) for(std::size_t m=0;m<4;++m)
      out[m]=out[m]*(v-x[i])+c[(k*count+i)*4+m];
    for(S& y:out) {y=exp(y);positive(y);}return out;
  }
};
template<class S,std::size_t Species=2> struct Response {
  std::array<std::array<S,Species+1>,Species+1> mobility{};
  std::array<S,Species> transport_enthalpy{};
  std::array<bool,Species> active_species{};
  S conductivity{},energy_scale{},eta{},electron_density{},b_thermal{};
  double maximum_solve_backward_error{};
};
template<class S> S screened_length(S T,S rho,S X,S Y3,S Z,S stiffness,bool include_ions) {
  positive(T);positive(stiffness);const Mixture<S> mixture(rho,X,Y3,Z);
  S inverse=4*pi*e2*mixture.ne/stiffness;
  if(include_ions) for(std::size_t i=0;i<mixture.x.size();++i) {
    const double z=mixture.z[i];const S radius=cbrt(3*z/(4*pi*mixture.ne));
    const S gamma=z*z*e2/(radius*kb*T);
    inverse+=4*pi*e2*mixture.n[i]*z*z/(kb*T*(1+3*gamma));
  }
  positive(inverse);return 1/sqrt(inverse);
}
} // namespace

struct ScreenedCollisionTransport::Data {
  std::vector<double> tx,ty,coeff,jx,jw,gx,gw;
  Polynomial base,extension;
  double density_min{},density_max{};
  explicit Data(const std::string& path) {
    std::ifstream f(path);std::string version;f>>version;
    require(version=="EMBER_COLLISION_TRANSPORT_V1","unsupported table format");
    auto count=[&]() {std::size_t n=0;f>>n;require(f.good() && n>=4 && n<=4096,"invalid table dimension");return n;};
    auto read=[&](std::vector<double>& a,std::size_t n) {
      a.resize(n);for(double& x:a) {f>>x;require(f.good() && std::isfinite(x),"invalid table value");}
    };
    auto increasing=[](const std::vector<double>& a,bool repeated) {
      for(std::size_t i=1;i<a.size();++i)
        require(repeated ? a[i]>=a[i-1] : a[i]>a[i-1],"unordered table grid");
      require(a.back()>a.front(),"empty table interval");
    };
    const auto nx=count(),ny=count();
    require(nx<=64 && ny<=64,"excessive pair spline size");
    read(tx,nx);read(ty,ny);increasing(tx,true);increasing(ty,true);
    require(tx[3]==tx.front() && tx[nx-4]==tx.back() && ty[3]==ty.front() && ty[ny-4]==ty.back(),"nonclamped spline");
    read(coeff,45*(nx-4)*(ny-4));
    for(auto* p:{&base,&extension}) {
      const auto n=count();p->order=count();
      require(p->order<=6,"unsupported ion polynomial order");
      read(p->x,n);read(p->c,(n-1)*p->order*4);increasing(p->x,false);
    }
    require(base.x.back()==extension.x.front(),"ion interpolation join differs");
    auto n=count();read(jx,n);read(jw,n);n=count();read(gx,n);read(gw,n);
    for(const auto* weights:{&jw,&gw}) for(double w:*weights) positive(w);
    for(const auto* nodes:{&jx,&gx}) {increasing(*nodes,false);require(nodes->front()>-1 && nodes->back()<1,"invalid quadrature node");}
    f>>std::ws;require(f.eof(),"trailing table content");
    density_min=fhalf(tx.front());density_max=fhalf(tx.back());
  }
  template<class S> Matrix<S> pair(S eta,S bthermal) const {
    if constexpr(std::is_same_v<S,Differential<6>>) {
      // This matrix depends on two variables. Differentiate those once,
      // then apply the chain rule to the six external state coordinates.
      using Local=Differential<2>;
      return compose_binary_response(pair(Local::variable(eta.value,0),
          Local::variable(bthermal.value,1)),eta,bthermal);
    }
    positive(bthermal);const auto bx=spline_basis(tx,eta),by=spline_basis(ty,log(bthermal));
    // Cubic splines have local support. Retain every nonzero value or
    // derivative, in the original summation order, and reuse the products
    // across the 45 coefficient surfaces. No small term is discarded.
    struct Weight {std::size_t index;S value;};
    std::vector<Weight> weights;weights.reserve(16);
    for(std::size_t a=0;a<bx.size();++a) if(!exactly_zero(bx[a]))
      for(std::size_t b=0;b<by.size();++b) if(!exactly_zero(by[b]))
        weights.push_back({a*by.size()+b,bx[a]*by[b]});
    Matrix<S> lower(9,9),ee(10,10);std::size_t k=0;
    for(std::size_t i=0;i<9;++i) {
      for(std::size_t j=0;j<=i;++j,++k) {
        S v=0;
        for(const auto& weight:weights)
          v+=weight.value*coeff[k*bx.size()*by.size()+weight.index];
        lower(i,j)=v;
      }
      const S diagonal=exp(lower(i,i));positive(diagonal);
      for(std::size_t j=0;j<i;++j) lower(i,j)*=diagonal;
      lower(i,i)=diagonal;
    }
    for(std::size_t i=0;i<9;++i) for(std::size_t j=0;j<=i;++j) {
      S v=0;for(std::size_t a=0;a<=j;++a) v+=lower(i,a)*lower(j,a);
      ee(i+1,j+1)=ee(j+1,i+1)=v;
    }
    return ee;
  }
  double fhalf(double eta) const {
    const double upper=std::sqrt(std::max(eta,0.)+80.);double value=0;
    for(std::size_t i=0;i<gx.size();++i) {
      const double t=(gx[i]+1)*upper/2;
      value+=gw[i]*upper*t*t*occupation(eta,t*t);
    }
    return value;
  }
  double eta(double ne,double T) const {
    const double target=ne/(std::pow(2*me*kb*T,1.5)/(2*pi*pi*hbar*hbar*hbar));
    double lo=tx.front(),hi=tx.back();
    require(target>=density_min && target<=density_max,"electron density outside pair-table degeneracy range");
    if(target==density_min)return lo;
    if(target==density_max)return hi;
    // Boltzmann and degenerate limits seed a bracketed Newton solve. Every
    // iterate uses the original quadrature; the table-domain check is exact.
    const double classical=std::log(target/(std::sqrt(pi)/2));
    double x=std::clamp(classical<0?classical:std::pow(1.5*target,2./3),lo,hi);
    constexpr double tolerance=8*std::numeric_limits<double>::epsilon();
    for(int iteration=0;iteration<64;++iteration) {
      const double upper=std::sqrt(std::max(x,0.)+80.);double f=0,df=0;
      for(std::size_t i=0;i<gx.size();++i) {
        const double t=(gx[i]+1)*upper/2,w=gw[i]*upper*t*t;
        f+=w*occupation(x,t*t);
        df+=w*occupation_derivative(x,t*t);
      }
      const double residual=f-target;
      if(std::abs(residual)<=tolerance*target)return x;
      if(residual>0)hi=x;else lo=x;
      if(hi-lo<=tolerance*std::max(1.,std::abs(x)))return (lo+hi)/2;
      const double trial=x-residual/df;
      x=std::isfinite(trial) && trial>lo && trial<hi?trial:(lo+hi)/2;
    }
    throw std::domain_error("collision transport: electron-density inversion did not converge");
  }

  template<class S> S chemical_eta(S ne,S T) const {
    const double eta_value=eta(value(ne),value(T));
    if constexpr(std::is_same_v<S,double>) return eta_value;
    else {
      const double upper=std::sqrt(std::max(eta_value,0.)+80.);double derivative=0;
      for(std::size_t i=0;i<gx.size();++i) {
        const double t=(gx[i]+1)*upper/2;
        derivative+=gw[i]*upper*t*t*occupation_derivative(eta_value,t*t);
      }
      positive(derivative);const double ratio=fhalf(eta_value)/derivative;
      S out(eta_value);
      for(std::size_t k=0;k<out.d.size();++k)
        out.d[k]=ratio*(ne.d[k]/ne.value-1.5*T.d[k]/T.value);
      return out;
    }
  }
  template<class S> Matrix<S> electron_ion(S eta,S bthermal) const {
    if constexpr(std::is_same_v<S,Differential<6>>) {
      using Local=Differential<2>;
      return compose_binary_response(electron_ion(Local::variable(eta.value,0),
          Local::variable(bthermal.value,1)),eta,bthermal);
    }
    const S length=positive_part(eta)+80.;
    std::vector<S> x(jx.size()),w(jx.size());
    S norm=0,mean=0,variance=0;
    for(std::size_t i=0;i<jx.size();++i) {
      x[i]=(jx[i]+1)*length/2;
      w[i]=jw[i]*pow(length/2,2.5)*occupation_derivative(eta,x[i]);
      norm+=w[i];mean+=w[i]*(x[i]-eta);
    }
    positive(norm);mean/=norm;
    for(std::size_t i=0;i<jx.size();++i) {w[i]/=norm;variance+=w[i]*pow(x[i]-eta-mean,2);}
    positive(variance);const S sd=sqrt(variance);
    Matrix<S> q(10,jx.size());std::array<S,10> alpha{},beta{};beta[1]=1;
    for(std::size_t i=0;i<jx.size();++i) {x[i]=(x[i]-eta-mean)/sd;q(0,i)=1;q(1,i)=x[i];}
    for(std::size_t n=1;n<9;++n) {
      for(std::size_t i=0;i<jx.size();++i) alpha[n]+=w[i]*x[i]*q(n,i)*q(n,i);
      S v=0;
      for(std::size_t i=0;i<jx.size();++i) {
        q(n+1,i)=(x[i]-alpha[n])*q(n,i)-beta[n]*q(n-1,i);
        v+=w[i]*q(n+1,i)*q(n+1,i);
      }
      positive(v);beta[n+1]=sqrt(v);
      for(std::size_t i=0;i<jx.size();++i) q(n+1,i)/=beta[n+1];
    }
    Matrix<S> out(10,10);const S upper=sqrt(length);
    for(std::size_t i=0;i<gx.size();++i) {
      const S t=(gx[i]+1)*upper/2,energy=t*t,u=(energy-eta-mean)/sd;
      const S weight=gw[i]*upper*t*occupation_derivative(eta,energy)*bracket(bthermal*energy);
      std::array<S,10> p{1,u};
      for(std::size_t n=1;n<9;++n) p[n+1]=((u-alpha[n])*p[n]-beta[n]*p[n-1])/beta[n+1];
      p[1]/=sd;
      for(std::size_t a=0;a<10;++a) for(std::size_t b=0;b<=a;++b) out(a,b)+=weight*p[a]*p[b];
    }
    for(std::size_t a=0;a<10;++a) for(std::size_t b=0;b<a;++b) out(b,a)=out(a,b);
    return out;
  }
  template<bool Metals=false,class S> Response<S,Metals?3:2> evaluate(S T,S rho,S X,S Y3,S Z,S length) const;
};

template<bool Metals,class S> Response<S,Metals?3:2> ScreenedCollisionTransport::Data::evaluate(S T,S rho,S X,S Y3,S Z,S length) const {
  positive(T);positive(length);const Mixture<S> mix(rho,X,Y3,Z);
#if defined(EMBER_TWO_METAL_COLLISIONS) && EMBER_TWO_METAL_COLLISIONS
  // Retain the physical mixture's table domain even when its representative
  // charges occupy a smaller interval.
  const double minimum_charge=value(X)>0 ? 1. : 2.;
  const double maximum_charge=value(Z)>0 ? 28. : 2.;
  const double denominator=kb*value(T)*value(length);
  const double minimum_strength=std::log(minimum_charge*minimum_charge*e2/denominator)/std::log(10.);
  const double maximum_strength=std::log(maximum_charge*maximum_charge*e2/denominator)/std::log(10.);
  require(minimum_strength>=base.x.front() && maximum_strength<=extension.x.back(),
      "full GS98 ion collision query outside table");
#endif
  constexpr std::size_t species=Metals?3:2;
  Response<S,species> result;result.electron_density=mix.ne;result.energy_scale=kb/mu*T;
  result.eta=chemical_eta(mix.ne,T);result.b_thermal=8*me*kb*T*length*length/(hbar*hbar);
  auto ee=pair(result.eta,result.b_thermal);
  const auto ei=electron_ion(result.eta,result.b_thermal);
  const S ratio=mix.ne/mix.nz2;
  for(std::size_t i=0;i<10;++i) for(std::size_t j=0;j<10;++j) ee(i,j)=ei(i,j)+ratio*ee(i,j);
  Matrix<S> identity(10,2);identity(0,0)=identity(1,1)=1;
  const auto response=solve(ee,identity,result.maximum_solve_backward_error);
  Matrix<S> reduced(2,2),id2(2,2);id2(0,0)=id2(1,1)=1;
  for(std::size_t i=0;i<2;++i) for(std::size_t j=0;j<2;++j) reduced(i,j)=(response(i,j)+response(j,i))/2;
  const auto electron_collision=solve(reduced,id2,result.maximum_solve_backward_error);
  const auto count=mix.x.size(),light=mix.independent.size();
  const bool moving_metals=Metals && value(Z)>0;
  const auto nc=light+std::size_t(moving_metals),size=nc+count+1;
  Matrix<S> velocity(count+1,nc),m(size,size);
  for(std::size_t c=0;c<light;++c) {
    const auto i=mix.independent[c];result.active_species[mix.output[c]]=true;
    velocity(i,c)=1;velocity(mix.reference,c)=-mix.x[i]/mix.x[mix.reference];
    velocity(count,c)=rho/mu*mix.x[i]*(mix.z[i]/mix.mass[i]-mix.z[mix.reference]/mix.mass[mix.reference])/mix.ne;
  }
  if(moving_metals) {
    result.active_species[2]=true;
    // One group velocity is the generalized coordinate. Its conjugate flux
    // is rho*Z*v_Z. Reference helium balances the total group mass; the
    // electron velocity cancels the resulting charge current exactly.
    S charge=0;
    for(std::size_t i=mix.reference+1;i<count;++i) {
      velocity(i,light)=1;
      charge+=mix.x[i]*(mix.z[i]/mix.mass[i]-mix.z[mix.reference]/mix.mass[mix.reference]);
    }
    velocity(mix.reference,light)=-Z/mix.x[mix.reference];
    velocity(count,light)=rho/mu*charge/mix.ne;
  }
  for(std::size_t i=0;i<count;++i) for(std::size_t j=i;j<count;++j) {
    const double mi=mix.mass[i],mj=mix.mass[j],total=mi+mj,redmass=mu*(mi*mj/total);
    const double interaction=mix.z[i]*mix.z[j]*e2;
    const S logstrength=log(interaction/(kb*T*length))/std::log(10.);
    const auto moments=(value(logstrength)<=base.x.back()?base:extension).eval(logstrength);
    const S omega=std::sqrt(2*pi/redmass)*interaction*interaction/pow(kb*T,1.5)*moments[0];
    const S drag=16./3*mix.n[i]*mix.n[j]*redmass*omega;
    const S z=1-.4*moments[1]/moments[0];
    const S zp=2.5-2*moments[1]/moments[0]+.4*moments[2]/moments[0],zpp=moments[3]/moments[0];
    require(value(drag)>=0 && value(zp)>=0 && value(zpp)>=0,"negative ion collision coefficient");
    if(i==j) {m(nc+i,nc+i)+=.16*drag*zpp;continue;}
    for(std::size_t a=0;a<nc;++a) {
      const S diff=velocity(i,a)-velocity(j,a);
      for(std::size_t b=0;b<nc;++b) m(a,b)+=drag*diff*(velocity(i,b)-velocity(j,b));
      m(a,nc+i)-=drag*z*mj/total*diff;m(a,nc+j)+=drag*z*mi/total*diff;
    }
    m(nc+i,nc+i)+=.4*drag*(3*mi*mi+mj*mj*zp+.8*mi*mj*zpp)/(total*total);
    m(nc+j,nc+j)+=.4*drag*(3*mj*mj+mi*mi*zp+.8*mi*mj*zpp)/(total*total);
    const S cross=-.4*drag*mi*mj*(3+zp-.8*zpp)/(total*total);
    m(nc+i,nc+j)+=cross;m(nc+j,nc+i)+=cross;
  }
  for(std::size_t a=0;a<nc;++a) for(std::size_t i=0;i<count;++i) m(nc+i,a)=m(a,nc+i);
  const S prefactor=2*me*me*e2*e2*mix.nz2/(3*pi*hbar*hbar*hbar);
  std::vector<S> mean(nc),relative(nc);
  for(std::size_t a=0;a<nc;++a) {
    for(std::size_t i=0;i<count;++i) mean[a]+=mix.n[i]*mix.z[i]*mix.z[i]/mix.nz2*velocity(i,a);
    relative[a]=velocity(count,a)-mean[a];
  }
  for(std::size_t i=0;i<count;++i) {
    const S weight=mix.n[i]*mix.z[i]*mix.z[i]/mix.nz2;
    for(std::size_t a=0;a<nc;++a) for(std::size_t b=0;b<nc;++b)
      m(a,b)+=prefactor*ei(0,0)*weight*(velocity(i,a)-mean[a])*(velocity(i,b)-mean[b]);
    m(nc+i,nc+i)+=1.2*prefactor*ei(0,0)*weight;
  }
  for(std::size_t a=0;a<nc;++a) {
    for(std::size_t b=0;b<nc;++b) m(a,b)+=prefactor*electron_collision(0,0)*relative[a]*relative[b];
    m(a,size-1)=m(size-1,a)=prefactor*(electron_collision(0,1)+electron_collision(1,0))/2*relative[a];
  }
  m(size-1,size-1)=prefactor*electron_collision(1,1);
  Matrix<S> g(nc+1,size),rhs(size,nc+1);
  for(std::size_t a=0;a<light;++a) g(a,a)=rho*mix.x[mix.independent[a]];
  if(moving_metals)g(light,light)=rho*Z;
  for(std::size_t i=0;i<count;++i) g(nc,nc+i)=mix.n[i]*kb*T/result.energy_scale;
  g(nc,size-1)=mix.ne*kb*T/result.energy_scale;
  for(std::size_t i=0;i<size;++i) for(std::size_t a=0;a<=nc;++a) rhs(i,a)=g(a,i);
  const auto variables=solve(m,rhs,result.maximum_solve_backward_error);
  Matrix<S> mobility(nc+1,nc+1);
  for(std::size_t a=0;a<=nc;++a) for(std::size_t b=0;b<=a;++b) {
    S v=0;for(std::size_t i=0;i<size;++i) v+=(g(a,i)*variables(i,b)+g(b,i)*variables(i,a))*T/2;
    mobility(a,b)=mobility(b,a)=v;
    const auto ia=a==nc?species:(a==light?2:mix.output[a]);
    const auto ib=b==nc?species:(b==light?2:mix.output[b]);
    result.mobility[ia][ib]=result.mobility[ib][ia]=v;
  }
  Matrix<S> normalized(nc+1,nc+1);std::vector<S> scale(nc+1);
  for(std::size_t a=0;a<=nc;++a) {positive(mobility(a,a));scale[a]=sqrt(mobility(a,a));}
  for(std::size_t a=0;a<=nc;++a) for(std::size_t b=0;b<=nc;++b) normalized(a,b)=mobility(a,b)/scale[a]/scale[b];
  const auto chol=cholesky(normalized);
  result.conductivity=pow(result.energy_scale/T*scale[nc]*chol(nc,nc),2);
  if(nc) {
    Matrix<S> a(nc,nc),b(nc,1);
    for(std::size_t i=0;i<nc;++i) {b(i,0)=normalized(i,nc);for(std::size_t j=0;j<nc;++j) a(i,j)=normalized(i,j);}
    const auto h=solve(a,b,result.maximum_solve_backward_error);
    for(std::size_t i=0;i<nc;++i) result.transport_enthalpy[i==light?2:mix.output[i]]=result.energy_scale*scale[nc]/scale[i]*h(i,0);
  }
  positive(result.conductivity);return result;
}
namespace {
template<class S,std::size_t Species> BasicCollisionTransportResponse<Species> values(const Response<S,Species>& in) {
  BasicCollisionTransportResponse<Species> out;
  out.active_species=in.active_species;
  out.maximum_solve_backward_error=in.maximum_solve_backward_error;
  for(std::size_t i=0;i<=Species;++i) for(std::size_t j=0;j<=Species;++j) out.mobility[i][j]=value(in.mobility[i][j]);
  for(std::size_t i=0;i<Species;++i) out.transport_enthalpy[i]=value(in.transport_enthalpy[i]);
  out.conductivity=value(in.conductivity);out.energy_scale=value(in.energy_scale);
  out.eta=value(in.eta);out.electron_density=value(in.electron_density);out.b_thermal=value(in.b_thermal);
  return out;
}
std::array<Differential<6>,6> coordinates(double T,double rho,double X,double Y3,double Z,double last) {
  positive(T);positive(rho);positive(last);
  std::array<Differential<6>,6> out{T,rho,X,Y3,Z,last};
  for(std::size_t i=0;i<6;++i) out[i].d[i]=(i==0 || i==1 || i==5)?out[i].value:1.;
  // The positive reference population remains fixed by normalization. At an
  // absent species, only the derivatives within the current active set exist.
  for(std::size_t i=2;i<=4;++i) if(out[i].value==0) out[i].d[i]=0;
  return out;
}
using CollisionState=std::array<std::uint64_t,6>;
CollisionState state_bits(const std::array<double,6>& state) {
  return std::bit_cast<CollisionState>(state);
}
struct StateHash {
  std::size_t operator()(const CollisionState& state) const {
    std::size_t result=0;
    for(auto x:state)result^=std::hash<std::uint64_t>{}(x)+0x9e3779b9+(result<<6)+(result>>2);
    return result;
  }
};
template<class Result> struct CollisionCache {
  // Keeping the table alive prevents its address from being reused for a
  // different table. Each thread retains only its most recently used table.
  std::shared_ptr<const void> owner;
  std::unordered_map<CollisionState,Result,StateHash> states;
  void use(const std::shared_ptr<const void>& table) {
    if(owner!=table){states.clear();owner=table;}
  }
  void save(const CollisionState& state,const Result& result) {
    if(states.size()>=8192)states.clear();
    states.emplace(state,result);
  }
};
}
ScreenedCollisionTransport::ScreenedCollisionTransport(const std::string& path):data_(std::make_shared<Data>(path)) {}
void ScreenedCollisionTransport::check_domain(double T,double rho,double X,double Y3,double Z,double length) const {
  positive(T);positive(length);const Mixture<double> mix(rho,X,Y3,Z);
  const double target=mix.ne/(std::pow(2*me*kb*T,1.5)/(2*pi*pi*hbar*hbar*hbar));
  require(target>=data_->density_min && target<=data_->density_max,
      "electron density outside pair-table degeneracy range");
  const double log_b=std::log(8*me*kb*T*length*length/(hbar*hbar));
  require(log_b>=data_->ty.front() && log_b<=data_->ty.back(),"electron pair query outside table");
  const double minimum_charge=X>0?1.:2.,maximum_charge=Z>0?28.:2.;
  const double low=std::log(minimum_charge*minimum_charge*e2/(kb*T*length))/std::log(10.);
  const double high=std::log(maximum_charge*maximum_charge*e2/(kb*T*length))/std::log(10.);
  require(low>=data_->base.x.front() && high<=data_->extension.x.back(),
      "full GS98 ion collision query outside table");
}
CollisionTransportResponse ScreenedCollisionTransport::eval(double T,double rho,double X,double Y3,double Z,double length) const {
  // Exact-state reuse only: no rounding, interpolation, or reuse across
  // tables. Keep scalar and differentiated evaluations separate so their
  // original floating-point results are preserved independently.
  thread_local CollisionCache<CollisionTransportResponse> cache;
  cache.use(data_);const auto state=state_bits({T,rho,X,Y3,Z,length});
  if(const auto found=cache.states.find(state);found!=cache.states.end())return found->second;
  const auto result=values(data_->evaluate(T,rho,X,Y3,Z,length));
  cache.save(state,result);return result;
}
CollisionTransportDerivatives ScreenedCollisionTransport::derivatives(double T,double rho,double X,double Y3,double Z,double length) const {
  thread_local CollisionCache<CollisionTransportDerivatives> cache;
  cache.use(data_);const auto state=state_bits({T,rho,X,Y3,Z,length});
  if(const auto found=cache.states.find(state);found!=cache.states.end())return found->second;
  const auto c=coordinates(T,rho,X,Y3,Z,length);
  const auto in=data_->evaluate(c[0],c[1],c[2],c[3],c[4],c[5]);
  CollisionTransportDerivatives out;out.value=values(in);out.defined={true,true,X>0,Y3>0,Z>0,true};
  for(std::size_t k=0;k<6;++k) {
    auto& d=out.partials[k];
    for(std::size_t i=0;i<3;++i) for(std::size_t j=0;j<3;++j) d.mobility[i][j]=in.mobility[i][j].d[k];
    for(std::size_t i=0;i<2;++i) d.transport_enthalpy[i]=in.transport_enthalpy[i].d[k];
    d.conductivity=in.conductivity.d[k];d.energy_scale=in.energy_scale.d[k];d.eta=in.eta.d[k];
    d.electron_density=in.electron_density.d[k];d.b_thermal=in.b_thermal.d[k];
  }
  cache.save(state,out);return out;
}
BulkMetalCollisionResponse ScreenedCollisionTransport::bulk_metal_eval(double T,double rho,
    double X,double Y3,double Z,double length) const {
  thread_local CollisionCache<BulkMetalCollisionResponse> cache;
  cache.use(data_);const auto state=state_bits({T,rho,X,Y3,Z,length});
  if(const auto found=cache.states.find(state);found!=cache.states.end())return found->second;
  const auto result=values(data_->evaluate<true>(T,rho,X,Y3,Z,length));
  cache.save(state,result);return result;
}
BulkMetalCollisionDerivatives ScreenedCollisionTransport::bulk_metal_derivatives(double T,double rho,
    double X,double Y3,double Z,double length) const {
  thread_local CollisionCache<BulkMetalCollisionDerivatives> cache;
  cache.use(data_);const auto state=state_bits({T,rho,X,Y3,Z,length});
  if(const auto found=cache.states.find(state);found!=cache.states.end())return found->second;
  const auto c=coordinates(T,rho,X,Y3,Z,length);
  const auto in=data_->evaluate<true>(c[0],c[1],c[2],c[3],c[4],c[5]);
  BulkMetalCollisionDerivatives out;out.value=values(in);out.defined={true,true,X>0,Y3>0,Z>0,true};
  for(std::size_t k=0;k<6;++k) {
    auto& d=out.partials[k];
    for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)d.mobility[i][j]=in.mobility[i][j].d[k];
    for(std::size_t i=0;i<3;++i)d.transport_enthalpy[i]=in.transport_enthalpy[i].d[k];
    d.conductivity=in.conductivity.d[k];d.energy_scale=in.energy_scale.d[k];d.eta=in.eta.d[k];
    d.electron_density=in.electron_density.d[k];d.b_thermal=in.b_thermal.d[k];
  }
  cache.save(state,out);return out;
}
double ScreenedCollisionTransport::screening_length(double T,double rho,double X,double Y3,double Z,double stiffness,bool include_ions) {
  return screened_length(T,rho,X,Y3,Z,stiffness,include_ions);
}
CollisionScreeningDerivatives ScreenedCollisionTransport::screening_derivatives(double T,double rho,double X,double Y3,double Z,double stiffness,bool include_ions) {
  const auto c=coordinates(T,rho,X,Y3,Z,stiffness);
  const auto in=screened_length(c[0],c[1],c[2],c[3],c[4],c[5],include_ions);
  return {in.value,in.d,{true,true,X>0,Y3>0,Z>0,true}};
}

struct CollisionTaylorCache::Slot {
  std::mutex mutex;
  bool valid{};
  std::shared_ptr<const ScreenedCollisionTransport::Data> owner;
  std::array<double,6> anchor{};
  BulkMetalCollisionDerivatives exact;
};
CollisionTaylorCache::CollisionTaylorCache(double radius,std::size_t faces,bool verify)
    :radius_(radius),faces_(faces),verify_(verify),slots_(std::make_unique<Slot[]>(faces)) {
  require(std::isfinite(radius) && radius>0 && radius<=.01,"collision Taylor radius must lie in (0,0.01]");
  require(faces>0,"collision reuse requires at least one face");
}
CollisionTaylorCache::~CollisionTaylorCache()=default;
CollisionTaylorCache::Statistics CollisionTaylorCache::statistics() const {
  std::lock_guard lock(worst_mutex_);
  return {hits_.load(),misses_.load(),verified_.load(),worst_};
}
namespace {
bool positive_response(const BulkMetalCollisionResponse& v) {
  for(double x:{v.conductivity,v.energy_scale,v.electron_density,v.b_thermal})
    if(!(std::isfinite(x) && x>0))return false;
  if(!std::isfinite(v.eta))return false;
  for(double x:v.transport_enthalpy)if(!std::isfinite(x))return false;
  // Positive definiteness on the active mass/heat subspace preserves
  // nonnegative dissipation. Normalize to avoid mixing coefficient scales.
  std::array<std::size_t,4> index{};std::size_t n=0;
  for(std::size_t i=0;i<3;++i)if(v.active_species[i])index[n++]=i;
  index[n++]=3;
  std::array<double,4> scale{};std::array<std::array<double,4>,4> lower{};
  for(std::size_t i=0;i<n;++i) {
    const double diagonal=v.mobility[index[i]][index[i]];
    if(!(std::isfinite(diagonal) && diagonal>0))return false;
    scale[i]=std::sqrt(diagonal);
    for(std::size_t j=0;j<=i;++j) {
      double x=v.mobility[index[i]][index[j]]/scale[i]/scale[j];
      for(std::size_t k=0;k<j;++k)x-=lower[i][k]*lower[j][k];
      if(!std::isfinite(x) || (i==j && x<=0))return false;
      lower[i][j]=i==j?std::sqrt(x):x/lower[j][j];
    }
  }
  return true;
}
template<std::size_t S> void shift(BasicCollisionTransportResponse<S>& v,const BasicCollisionTransportPartial<S>& p,double d) {
  for(std::size_t i=0;i<=S;++i)for(std::size_t j=0;j<=S;++j)v.mobility[i][j]+=p.mobility[i][j]*d;
  for(std::size_t i=0;i<S;++i)v.transport_enthalpy[i]+=p.transport_enthalpy[i]*d;
  v.conductivity+=p.conductivity*d;v.energy_scale+=p.energy_scale*d;v.eta+=p.eta*d;
  v.electron_density+=p.electron_density*d;v.b_thermal+=p.b_thermal*d;
}
// Largest error relative to the scale of each coefficient group.
double relative_difference(const BulkMetalCollisionResponse& a,const BulkMetalCollisionResponse& b) {
  double m=0,dm=0;
  for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j){m=std::max(m,std::abs(b.mobility[i][j]));dm=std::max(dm,std::abs(a.mobility[i][j]-b.mobility[i][j]));}
  double e=m>0?dm/m:0,h=0,dh=0;
  for(std::size_t i=0;i<3;++i){h=std::max(h,std::abs(b.transport_enthalpy[i]));dh=std::max(dh,std::abs(a.transport_enthalpy[i]-b.transport_enthalpy[i]));}
  if(h>0)e=std::max(e,dh/h);
  if(b.conductivity!=0)e=std::max(e,std::abs(a.conductivity/b.conductivity-1));
  return e;
}
} // namespace
BulkMetalCollisionDerivatives CollisionTaylorCache::lookup(const ScreenedCollisionTransport& collision,
    std::size_t face,const std::array<double,6>& x,double T,double rho,double X,double Y3,double Z,double length) const {
  collision.check_domain(T,rho,X,Y3,Z,length);
  if(face>=faces_) {misses_.fetch_add(1,std::memory_order_relaxed);return collision.bulk_metal_derivatives(T,rho,X,Y3,Z,length);}
  auto& slot=slots_[face];
  std::unique_lock lock(slot.mutex);
  if(slot.valid && slot.owner==collision.data_) {
    bool inside=true;std::array<double,6> d{};
    for(std::size_t k=0;k<6 && inside;++k) {
      d[k]=x[k]-slot.anchor[k];
      inside=std::abs(d[k])<=radius_ && (d[k]==0 || slot.exact.defined[k]);
      if(k>=2 && k<=4)
        inside=inside && (x[k]>0)==(slot.anchor[k]>0) && std::abs(d[k])<=.01*slot.anchor[k];
    }
    const double helium=1-slot.anchor[2]-slot.anchor[3]-slot.anchor[4];
    inside=inside && std::abs(d[2]+d[3]+d[4])<=.01*helium;
    if(inside) {
      auto out=slot.exact;
      for(std::size_t k=0;k<6;++k)if(d[k]!=0)shift(out.value,slot.exact.partials[k],d[k]);
      if(positive_response(out.value)) {
        lock.unlock();
        hits_.fetch_add(1,std::memory_order_relaxed);
        if(verify_) {
          const auto exact=collision.bulk_metal_derivatives(T,rho,X,Y3,Z,length);
          const double e=relative_difference(out.value,exact.value);
          verified_.fetch_add(1,std::memory_order_relaxed);
          std::lock_guard w(worst_mutex_);worst_=std::max(worst_,e);
        }
        return out;
      }
    }
  }
  slot.exact=collision.bulk_metal_derivatives(T,rho,X,Y3,Z,length);
  slot.anchor=x;slot.owner=collision.data_;slot.valid=true;
  misses_.fetch_add(1,std::memory_order_relaxed);
  return slot.exact;
}
BulkMetalCollisionDerivatives CollisionTaylorCache::derivatives(const ScreenedCollisionTransport& collision,
    std::size_t face,double T,double rho,double X,double Y3,double Z,double length) const {
  return lookup(collision,face,{std::log(T),std::log(rho),X,Y3,Z,std::log(length)},T,rho,X,Y3,Z,length);
}
BulkMetalCollisionResponse CollisionTaylorCache::value(const ScreenedCollisionTransport& collision,
    std::size_t face,double T,double rho,double X,double Y3,double Z,double length) const {
  return lookup(collision,face,{std::log(T),std::log(rho),X,Y3,Z,std::log(length)},T,rho,X,Y3,Z,length).value;
}
} // namespace ember

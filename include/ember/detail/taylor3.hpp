#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace ember::detail {
// Taylor coefficients through total degree three, divided by multi-index
// factorials. Used for an implicit thermodynamic potential, whose pressure,
// heat and composition responses must all differentiate the same function.
template<std::size_t N> struct Taylor3 {
  static constexpr std::size_t count=(N+1)*(N+2)*(N+3)/6;
  static constexpr std::size_t codes=std::size_t(1)<<(2*N);
  using Powers=std::array<unsigned,N>;
  struct Layout {
    std::array<Powers,count> powers{};
    std::array<int,codes> lookup{};
    std::array<unsigned,count> degree{},code{};
    constexpr Layout() {
      lookup.fill(-1);std::size_t index=0;
      for(unsigned order=0;order<=3;++order)for(unsigned encoded=0;encoded<codes;++encoded) {
        Powers p{};unsigned d=0;
        for(std::size_t i=0;i<N;++i){p[i]=(encoded>>(2*i))&3;d+=p[i];}
        if(d==order){powers[index]=p;degree[index]=d;code[index]=encoded;lookup[encoded]=static_cast<int>(index++);}
      }
    }
  };
  inline static constexpr Layout layout{};
  struct Products {
    struct Term {unsigned a,b,out;};
    // Number of ordered pairs of multi-indices whose total degree is <=3.
    std::array<Term,(2*N+1)*(2*N+2)*(2*N+3)/6> terms{};
    constexpr Products(){
      std::size_t next=0;
      for(unsigned i=0;i<count;++i)for(unsigned j=0;j<count;++j)
        if(layout.degree[i]+layout.degree[j]<=3)
          terms[next++]={i,j,static_cast<unsigned>(layout.lookup[layout.code[i]+layout.code[j]])};
    }
  };
  inline static constexpr Products products{};
  std::array<double,count> c{};
  Taylor3()=default;
  Taylor3(double x){c[0]=x;}
  double value()const{return c[0];}
  static Taylor3 variable(double value,std::size_t index) {
    Taylor3 x(value);x.c[layout.lookup[std::size_t(1)<<(2*index)]]=1;return x;
  }
  double derivative(Powers p)const {
    unsigned degree=0,code=0,factor=1;
    for(std::size_t i=0;i<N;++i){degree+=p[i];code+=p[i]<<(2*i);if(p[i]==2)factor*=2;if(p[i]==3)factor*=6;}
    if(degree>3)throw std::invalid_argument("Taylor derivative exceeds total degree three");
    return c[layout.lookup[code]]*factor;
  }
  friend Taylor3 operator+(const Taylor3& a,const Taylor3& b){Taylor3 v;for(std::size_t i=0;i<count;++i)v.c[i]=a.c[i]+b.c[i];return v;}
  friend Taylor3 operator+(const Taylor3& a,double b){auto v=a;v.c[0]+=b;return v;}
  friend Taylor3 operator+(double a,const Taylor3& b){return b+a;}
  friend Taylor3 operator-(const Taylor3& a,const Taylor3& b){Taylor3 v;for(std::size_t i=0;i<count;++i)v.c[i]=a.c[i]-b.c[i];return v;}
  friend Taylor3 operator-(const Taylor3& a,double b){auto v=a;v.c[0]-=b;return v;}
  friend Taylor3 operator-(double a,const Taylor3& b){return -b+a;}
  friend Taylor3 operator-(const Taylor3& a){return Taylor3(0)-a;}
  friend Taylor3 operator*(const Taylor3& a,const Taylor3& b){
    Taylor3 v;
    for(const auto term:products.terms)v.c[term.out]+=a.c[term.a]*b.c[term.b];
    return v;
  }
  friend Taylor3 operator*(const Taylor3& a,double b){Taylor3 v;for(std::size_t i=0;i<count;++i)v.c[i]=a.c[i]*b;return v;}
  friend Taylor3 operator*(double a,const Taylor3& b){return b*a;}
  friend Taylor3 operator/(const Taylor3& a,const Taylor3& b){
    if(b.value()==0)throw std::domain_error("Taylor division by zero");
    Taylor3 h=b/b.value();h.c[0]=0;const Taylor3 h2=h*h;
    return (a*(Taylor3(1)-h+h2-h2*h))/b.value();
  }
  friend Taylor3 operator/(const Taylor3& a,double b){Taylor3 v;for(std::size_t i=0;i<count;++i)v.c[i]=a.c[i]/b;return v;}
};
template<std::size_t N>Taylor3<N> exp(const Taylor3<N>& a){auto h=a;h.c[0]=0;const auto h2=h*h;return std::exp(a.value())*(Taylor3<N>(1)+h+h2/2+h2*h/6);}
template<std::size_t N>Taylor3<N> log(const Taylor3<N>& a){
  if(a.value()<=0)throw std::domain_error("Taylor logarithm of nonpositive value");
  auto h=a/a.value();h.c[0]=0;const auto h2=h*h;return Taylor3<N>(std::log(a.value()))+h-h2/2+h2*h/3;
}
} // namespace ember::detail

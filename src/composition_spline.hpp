#pragma once
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

namespace ember::detail {
// Eliminate the two endpoint second derivatives with the not-a-knot
// conditions, then factor the remaining tridiagonal system once per axis.
class SplineFactor {
public:
  explicit SplineFactor(std::span<const double> x):h_(x.size()-1),d_(x.size()-2),
      upper_(x.size()-2),lower_factor_(x.size()-2) {
    if(x.size()<4)throw std::invalid_argument("composition spline needs four nodes");
    for(std::size_t i=0;i<h_.size();++i)h_[i]=x[i+1]-x[i];
    for(std::size_t k=0;k<d_.size();++k) {
      d_[k]=2*(h_[k]+h_[k+1]);upper_[k]=h_[k+1];
    }
    d_.front()+=h_[0]*(h_[0]+h_[1])/h_[1];
    upper_.front()-=h_[0]*h_[0]/h_[1];
    const auto last=h_.size()-1;
    d_.back()+=h_[last]*(h_[last-1]+h_[last])/h_[last-1];
    for(std::size_t k=1;k<d_.size();++k) {
      double lower=h_[k];
      if(k+1==d_.size())lower-=h_[last]*h_[last]/h_[last-1];
      lower_factor_[k]=lower/d_[k-1];d_[k]-=lower_factor_[k]*upper_[k-1];
    }
    for(double d:d_)if(!std::isfinite(d) || d==0)
      throw std::runtime_error("singular composition spline");
  }
  // Reuse caller-owned scratch. No heap allocation during source-node solves.
  void slopes(std::span<const double> v,std::span<double> m,std::span<double> s) const {
    const auto n=v.size();
    for(std::size_t i=1;i+1<n;++i)
      s[i]=6*((v[i+1]-v[i])/h_[i]-(v[i]-v[i-1])/h_[i-1]);
    for(std::size_t k=1;k<d_.size();++k)s[k+1]-=lower_factor_[k]*s[k];
    s[n-2]/=d_.back();
    for(std::size_t k=d_.size()-1;k>0;--k)s[k]=(s[k]-upper_[k-1]*s[k+1])/d_[k-1];
    s[0]=((h_[0]+h_[1])*s[1]-h_[0]*s[2])/h_[1];
    s[n-1]=((h_[n-3]+h_[n-2])*s[n-2]-h_[n-2]*s[n-3])/h_[n-3];
    for(std::size_t i=0;i+1<n;++i)m[i]=(v[i+1]-v[i])/h_[i]-h_[i]*(2*s[i]+s[i+1])/6;
    m[n-1]=(v[n-1]-v[n-2])/h_[n-2]+h_[n-2]*(s[n-2]+2*s[n-1])/6;
  }
private:
  std::vector<double> h_,d_,upper_,lower_factor_;
};

} // namespace ember::detail

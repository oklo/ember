#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

namespace ember::detail {
template<std::size_t N> using FluxVector=std::array<double,N>;
template<std::size_t N> using FluxMatrix=std::array<FluxVector<N>,N>;
template<std::size_t N> void flux_finite(const FluxVector<N>& v) {
  for(double x:v)if(!std::isfinite(x))throw std::domain_error("flux chain: nonfinite vector");
}
template<std::size_t N> void flux_finite(const FluxMatrix<N>& a) {for(const auto& row:a)flux_finite(row);}
template<std::size_t N> FluxMatrix<N> flux_add(FluxMatrix<N> a,const FluxMatrix<N>& b) {
  for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<N;++j)a[i][j]+=b[i][j];return a;
}
template<std::size_t N> FluxVector<N> flux_add(FluxVector<N> a,const FluxVector<N>& b) {
  for(std::size_t i=0;i<N;++i)a[i]+=b[i];return a;
}
template<std::size_t N> FluxMatrix<N> flux_product(const FluxMatrix<N>& a,const FluxMatrix<N>& b) {
  FluxMatrix<N> c{};
  for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<N;++j)for(std::size_t k=0;k<N;++k)c[i][j]+=a[i][k]*b[k][j];
  return c;
}
template<std::size_t N> FluxVector<N> flux_product(const FluxMatrix<N>& a,const FluxVector<N>& b) {
  FluxVector<N> c{};for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<N;++j)c[i]+=a[i][j]*b[j];return c;
}
template<std::size_t N> FluxMatrix<N> flux_inverse(FluxMatrix<N> a) {
  flux_finite(a);FluxMatrix<N> rhs{};
  for(std::size_t i=0;i<N;++i) {
    double s=0;for(double x:a[i])s=std::max(s,std::abs(x));
    if(!(s>0))throw std::domain_error("flux chain: singular block");
    rhs[i][i]=1/s;for(double& x:a[i])x/=s;
  }
  for(std::size_t k=0;k<N;++k) {
    std::size_t pivot=k;for(std::size_t i=k+1;i<N;++i)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;
    std::swap(a[k],a[pivot]);std::swap(rhs[k],rhs[pivot]);
    if(a[k][k]==0 || !std::isfinite(a[k][k]))throw std::domain_error("flux chain: singular pivot");
    for(std::size_t i=k+1;i<N;++i) {
      const double f=a[i][k]/a[k][k];
      for(std::size_t j=k+1;j<N;++j)a[i][j]-=f*a[k][j];
      for(std::size_t j=0;j<N;++j)rhs[i][j]-=f*rhs[k][j];
    }
  }
  for(std::size_t i=N;i-->0;)for(std::size_t j=0;j<N;++j) {
    for(std::size_t k=i+1;k<N;++k)rhs[i][j]-=a[i][k]*rhs[k][j];
    rhs[i][j]/=a[i][i];
  }
  flux_finite(rhs);return rhs;
}
// E_i x_i + f_i-f_(i-1)=b_i; f_i=A_i x_i-B_i x_(i+1)+d_i.
// Keep the finite mass/reaction block in the elimination instead of
// subtracting nearly equal diffusion diagonals. Shared by H/He and CN.
template<std::size_t N> std::vector<FluxVector<N>> solve_flux_chain(
    std::span<const FluxMatrix<N>> E,std::span<const FluxMatrix<N>> A,
    std::span<const FluxMatrix<N>> B,std::span<const FluxVector<N>> b,
    std::span<const FluxVector<N>> d,
    std::vector<FluxVector<N>>* face_fluxes=nullptr) {
  const auto n=E.size();
  if(!n || b.size()!=n || A.size()!=n-1 || B.size()!=n-1 || d.size()!=n-1)
    throw std::invalid_argument("flux chain: invalid dimensions");
  for(const auto& a:E)flux_finite(a);for(const auto& a:A)flux_finite(a);for(const auto& a:B)flux_finite(a);
  for(const auto& v:b)flux_finite(v);for(const auto& v:d)flux_finite(v);
  std::vector<FluxMatrix<N>> effective(E.begin(),E.end()),transfer(n-1);
  std::vector<FluxVector<N>> right(b.begin(),b.end()),shift(n-1),answer(n);
  for(std::size_t i=0;i+1<n;++i) {
    const auto inv=flux_inverse(flux_add(effective[i],A[i]));
    const auto mass_transfer=flux_product(effective[i],inv),flux_transfer=flux_product(A[i],inv);
    effective[i+1]=flux_add(effective[i+1],flux_product(mass_transfer,B[i]));
    right[i+1]=flux_add(right[i+1],flux_add(flux_product(flux_transfer,right[i]),flux_product(mass_transfer,d[i])));
    transfer[i]=flux_product(inv,B[i]);auto q=right[i];for(std::size_t k=0;k<N;++k)q[k]-=d[i][k];
    shift[i]=flux_product(inv,q);
  }
  answer.back()=flux_product(flux_inverse(effective.back()),right.back());
  for(std::size_t i=n-1;i-->0;)answer[i]=flux_add(shift[i],flux_product(transfer[i],answer[i+1]));
  for(const auto& x:answer)flux_finite(x);
  if(face_fluxes) {
    face_fluxes->resize(n-1);
    for(std::size_t i=0;i+1<n;++i) {
      // The eliminated equation is E_effective x_i + f_i = b_effective.
      // These blocks retain the finite storage/reaction terms. Computing
      // A*x_i-B*x_(i+1) instead loses a finite flux when stiff diffusion
      // makes the represented abundances equal to machine precision.
      const auto storage=flux_product(effective[i],answer[i]);
      for(std::size_t k=0;k<N;++k) {
        double direct=d[i][k],direct_scale=std::abs(d[i][k]);
        double storage_scale=std::abs(right[i][k]);
        for(std::size_t j=0;j<N;++j) {
          const double left=A[i][k][j]*answer[i][j],next=B[i][k][j]*answer[i+1][j];
          direct+=left-next;direct_scale+=std::abs(left)+std::abs(next);
          storage_scale+=std::abs(effective[i][k][j]*answer[i][j]);
        }
        // Weak exchange is best evaluated from its constitutive equation;
        // strong exchange is best recovered from finite conservation terms.
        // Choose the smaller sum of magnitudes of the arithmetic terms.
        // This changes no convergence or physical acceptance tolerance.
        (*face_fluxes)[i][k]=direct_scale<storage_scale?direct:right[i][k]-storage[k];
      }
      flux_finite((*face_fluxes)[i]);
    }
  }
  return answer;
}
} // namespace ember::detail

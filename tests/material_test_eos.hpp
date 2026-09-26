#pragma once
#include "ember/material_transport.hpp"
#include <cmath>
#include <stdexcept>

// Analytic test material only: classical ideal ions, prescribed binding terms,
// and optional cold degenerate electrons, in units Rgas=1. The audit checks
// derivatives by finite differences and diffusion against independent solves.
struct AnalyticMaterial {
  std::size_t n{3};
  ember::MaterialVector masses{},charges{},binding{};
  std::vector<double> density;
  double degeneracy{},energy_scale{1};
  ember::MaterialPoint operator()(std::size_t cell,const ember::MaterialVector& q,bool derivatives) const {
    using V=ember::MaterialVector;using M=ember::MaterialMatrix;
    V x{};x[n-1]=1;for(std::size_t i=0;i+1<n;++i){x[i]=q[i];x[n-1]-=q[i];}
    const double T=std::exp(q[n-1]),rho=density[cell];
    if(!(T>0) || !std::isfinite(T))throw std::domain_error("analytic material: invalid temperature");
    double cv=0,ye=0,cold=0,entropy=0;
    for(std::size_t i=0;i<n;++i) {
      if(x[i]<0 || (derivatives && x[i]==0))throw std::domain_error("analytic material: invalid composition");
      cv+=1.5*x[i]/masses[i];ye+=x[i]*charges[i]/masses[i];cold+=x[i]*binding[i];
      if(x[i]>0)entropy-=x[i]/masses[i]*std::log(rho*x[i]/masses[i]);
    }
    const double factor=degeneracy*std::pow(rho,2.0/3);
    cold+=factor*std::pow(ye,5.0/3);entropy+=cv*std::log(T);
    ember::MaterialPoint p;p.conserved=x;p.conserved[n-1]=(cv*T+cold)/energy_scale;p.entropy=entropy;
    if(!derivatives)return p;
    V dy{},ec{};M h{};
    for(std::size_t i=0;i+1<n;++i) {
      dy[i]=charges[i]/masses[i]-charges[n-1]/masses[n-1];
      const double dcv=1.5/masses[i]-1.5/masses[n-1];
      const double dcold=binding[i]-binding[n-1]+5.0/3*factor*std::pow(ye,2.0/3)*dy[i];
      ec[i]=dcv*T+dcold;
      p.potential[i]=dcv*(1-std::log(T))+(std::log(rho*x[i]/masses[i])+1)/masses[i]
          -(std::log(rho*x[n-1]/masses[n-1])+1)/masses[n-1]+dcold/T;
    }
    for(std::size_t i=0;i+1<n;++i) {
      for(std::size_t j=0;j+1<n;++j)h[i][j]=(i==j?1/(masses[i]*x[i]):0)+1/(masses[n-1]*x[n-1])
          +10.0/9*factor*std::pow(ye,-1.0/3)*dy[i]*dy[j]/T+ec[i]*ec[j]/(cv*T*T);
      h[i][n-1]=h[n-1][i]=-energy_scale*ec[i]/(cv*T*T);
      p.primitive_from_conserved[i][i]=1;p.primitive_from_conserved[n-1][i]=-ec[i]/(cv*T);
    }
    h[n-1][n-1]=energy_scale*energy_scale/(cv*T*T);
    p.capacity=ember::material_positive_inverse(h,n);p.potential[n-1]=-energy_scale/T;
    p.primitive_from_conserved[n-1][n-1]=energy_scale/(cv*T);return p;
  }
};

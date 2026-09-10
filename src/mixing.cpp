#include "ember/evolution.hpp"
#include "ember/convection.hpp"
#include "ember/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ember {
std::vector<double> secular_mixing_diffusivities(const Model& m,const Physics& p) {
  nodal_mass_weights(m);
  if(!std::isfinite(p.alpha_semiconvection+p.alpha_thermohaline) || p.alpha_semiconvection<0 || p.alpha_thermohaline<0)
    throw std::invalid_argument("secular mixing: invalid efficiency");
  std::vector<double> D(m.size()-1);
  if(p.alpha_semiconvection==0 && p.alpha_thermohaline==0)return D;
  if(p.criterion!=ConvectiveCriterion::ledoux || !p.eos || !p.opacity || m.comp.size()!=m.size())
    throw std::invalid_argument("secular mixing requires Ledoux criterion and material physics");
  std::vector<EosState> state;std::vector<double> opacity;
  for(std::size_t i=0;i<m.size();++i) {
    state.push_back(p.eos->eval(m.T(i),m.rho(i),m.comp[i]));
    opacity.push_back(p.opacity->eval(m.T(i),m.rho(i),m.comp[i]).kappa);
  }
  for(std::size_t i=0;i+1<m.size();++i) {
    const double T=.5*(m.T(i)+m.T(i+1)),rho=.5*(m.rho(i)+m.rho(i+1));
    const double P=.5*(state[i].P+state[i+1].P),cp=.5*(state[i].cp+state[i+1].cp);
    const double delta=.5*(state[i].delta+state[i+1].delta),ad=.5*(state[i].grad_ad+state[i+1].grad_ad);
    const double k=.5*(opacity[i]+opacity[i+1]),mass=.5*(m.m[i]+m.m[i+1]),L=.5*(m.y[i].L+m.y[i+1].L);
    const double rad=3*k*L*P/(16*M_PI*constants::a_rad*constants::c*constants::G*mass*std::pow(T,4));
    const double B=composition_buoyancy(*p.eos,T,P,delta,std::log(state[i+1].P)-std::log(state[i].P),m.comp[i],m.comp[i+1],rho).B;
    const double thermal=4*constants::a_rad*constants::c*std::pow(T,3)/(3*k*rho*rho*cp);
    // MLT handles Ledoux-unstable faces. These prescriptions apply only
    // where heat diffusion stabilizes one of the two buoyancy contributions.
    if(rad>ad+B)continue;
    if(B>0 && rad>ad && p.alpha_semiconvection>0)
      D[i]=p.alpha_semiconvection*thermal/6*(rad-ad)/(ad+B-rad);
    else if(B<0 && rad<ad && p.alpha_thermohaline>0)
      D[i]=-1.5*p.alpha_thermohaline*thermal*B/(ad-rad);
    if(!std::isfinite(D[i]) || D[i]<0)
      throw std::domain_error("secular mixing: unresolved neutral boundary or invalid diffusivity");
  }
  return D;
}

namespace {
using Vector=std::array<double,2>;
using Matrix=std::array<double,4>;
Vector multiply(const Matrix& a,const Vector& b) {
  return {a[0]*b[0]+a[1]*b[1],a[2]*b[0]+a[3]*b[1]};
}
Matrix product(const Matrix& a,const Matrix& b) {
  return {a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]};
}
Matrix inverse(Matrix a) {
  const double det=a[0]*a[3]-a[1]*a[2];
  if(!std::isfinite(det) || det<=0)throw std::runtime_error("burn_and_transport: singular chemical Jacobian");
  return {a[3]/det,-a[1]/det,-a[2]/det,a[0]/det};
}
// g*(E+g I)^-1, scaled before inversion. This lets the Schur complement
// retain the mass/reaction matrix even when diffusion is arbitrarily stiff:
// E_next += g E/(E+g I), instead of subtracting two large diffusion terms.
Matrix transfer(const Matrix& E,double g) {
  const double scale=std::max({g,std::abs(E[0]),std::abs(E[1]),std::abs(E[2]),std::abs(E[3])});
  const double ratio=g/scale;
  auto a=inverse({E[0]/scale+ratio,E[1]/scale,E[2]/scale,E[3]/scale+ratio});
  for(auto& v:a)v*=ratio;return a;
}
}

std::vector<Composition> burn_and_transport(const Model& thermal,const Model& previous,
    const Nuclear& nuclear,const MixingRegions& regions,const std::vector<double>& D,double dt,double tolerance) {
  const auto weights=nodal_mass_weights(thermal);
  if(D.size()+1!=thermal.size() || thermal.m!=previous.m || thermal.size()!=previous.size()
      || previous.comp.size()!=thermal.size() || thermal.comp.size()!=thermal.size()
      || !std::isfinite(dt+tolerance) || dt<=0 || tolerance<=0)
    throw std::invalid_argument("burn_and_transport: invalid mesh, step or diffusivity dimensions");
  for(double d:D)if(!std::isfinite(d) || d<0)throw std::invalid_argument("burn_and_transport: invalid diffusivity");
  if(std::all_of(D.begin(),D.end(),[](double d){return d==0;}))
    return burn_and_mix(thermal,previous,nuclear,regions,dt,tolerance);
  std::vector<Composition> old;std::vector<double> mass,g;std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())throw std::invalid_argument("burn_and_transport: invalid partition");
    next=end;Composition c{};c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=previous.comp[begin].metal_inventory;
    double w=0;
    for(std::size_t i=begin;i<end;++i) {
      const auto& v=previous.comp[i];
      if(v.basis!=c.basis || v.metal_inventory!=previous.comp.front().metal_inventory
          || std::abs(v.sum()-1)>1e-10)
        throw std::invalid_argument("burn_and_transport: inconsistent composition metadata");
      for(std::size_t j=0;j<NSPEC;++j) {
        if(!std::isfinite(v.X[j]) || v.X[j]<0)throw std::invalid_argument("burn_and_transport: invalid abundance");
        // All five metal carriers must be constant: only the active H/He
        // abundances are transported in this reduced network.
        if(j>=3 && std::abs(v.X[j]-previous.comp.front().X[j])>1e-12)
          throw std::invalid_argument("burn_and_transport: metal gradients require a larger network");
        c.X[j]+=weights[i]*v.X[j];
      }
      w+=weights[i];
    }
    for(auto& v:c.X)v/=w;
    mass.push_back(w/thermal.M);old.push_back(c);
    if(end<thermal.size()) {
      const double r=.5*(thermal.r(end-1)+thermal.r(end)),rho=.5*(thermal.rho(end-1)+thermal.rho(end));
      const double area_mass=4*M_PI*r*r*rho;
      const double coupling=dt/thermal.M*area_mass*area_mass*D[end-1]/(thermal.m[end]-thermal.m[end-1]);
      if(!std::isfinite(coupling))throw std::domain_error("burn_and_transport: unrepresentable diffusion coupling");
      g.push_back(coupling);
    }
  }
  if(next!=thermal.size())throw std::invalid_argument("burn_and_transport: incomplete partition");
  auto current=old;const auto count=regions.size();
  for(int iteration=0;iteration<60;++iteration) {
    std::vector<Matrix> E(count);std::vector<Vector> rhs(count),answer(count);
    for(std::size_t i=0;i<count;++i) {
      E[i]={mass[i],0,0,mass[i]};rhs[i]={mass[i]*old[i].X[0],mass[i]*old[i].X[1]};
      for(std::size_t j=regions[i].first;j<regions[i].second;++j) {
        const auto rates=nuclear.composition_response(thermal.T(j),thermal.rho(j),current[i]);
        const double factor=dt*weights[j]/thermal.M;
        for(int row=0;row<2;++row) {
          rhs[i][row]+=factor*rates.state.dXdt[row];
          for(int col=0;col<2;++col) {
            const double J=rates.d_dXdt_dX[row][col]-rates.d_dXdt_dX[row][2];
            E[i][2*row+col]-=factor*J;rhs[i][row]-=factor*J*current[i].X[col];
          }
        }
      }
      if(i && g[i-1]>0) {
        const auto W=transfer(E[i-1],g[i-1]),addition=product(W,E[i-1]);
        const auto b=multiply(W,rhs[i-1]);
        for(int j=0;j<4;++j)E[i][j]+=addition[j];
        for(int j=0;j<2;++j)rhs[i][j]+=b[j];
      }
    }
    for(std::size_t i=count;i-->0;) {
      const double coupling=i+1<count?g[i]:0;
      const double scale=std::max({coupling,std::abs(E[i][0]),std::abs(E[i][1]),std::abs(E[i][2]),std::abs(E[i][3])});
      const double q=coupling/scale;
      Vector b{rhs[i][0]/scale,rhs[i][1]/scale};
      if(i+1<count)for(int j=0;j<2;++j)b[j]+=q*answer[i+1][j];
      answer[i]=multiply(inverse({E[i][0]/scale+q,E[i][1]/scale,E[i][2]/scale,E[i][3]/scale+q}),b);
    }
    double damping=1;
    for(std::size_t i=0;i<count;++i) {
      const std::array<double,3> change{answer[i][0]-current[i].X[0],answer[i][1]-current[i].X[1],
          current[i].X[0]+current[i].X[1]-answer[i][0]-answer[i][1]};
      for(int j=0;j<3;++j)if(change[j]<-current[i].X[j])damping=std::min(damping,-.9*current[i].X[j]/change[j]);
    }
    if(!std::isfinite(damping) || damping<=0)throw std::runtime_error("burn_and_transport: no positive Newton step");
    double change=0;
    for(std::size_t i=0;i<count;++i) {
      for(int j=0;j<2;++j) {
        const double delta=damping*(answer[i][j]-current[i].X[j]);
        current[i].X[j]+=delta;change=std::max(change,std::abs(delta));
      }
      current[i].X[2]=1-current[i].Z()-current[i].X[0]-current[i].X[1];
      if(!std::isfinite(current[i].sum()) || *std::min_element(current[i].X.begin(),current[i].X.end())<0)
        throw std::runtime_error("burn_and_transport: invalid chemical iterate");
    }
    if(change<tolerance) {
      Vector balance{};
      for(std::size_t i=0;i<count;++i) {
        for(int k=0;k<2;++k)balance[k]+=mass[i]*(current[i].X[k]-old[i].X[k]);
        for(std::size_t j=regions[i].first;j<regions[i].second;++j) {
          const auto rates=nuclear.eval(thermal.T(j),thermal.rho(j),current[i]);
          for(int k=0;k<2;++k)balance[k]-=dt*weights[j]/thermal.M*rates.dXdt[k];
        }
      }
      if(std::max(std::abs(balance[0]),std::abs(balance[1]))>10*tolerance)
        throw std::runtime_error("burn_and_transport: integrated reaction/transport imbalance");
      auto result=previous.comp;
      for(std::size_t i=0;i<count;++i)for(std::size_t j=regions[i].first;j<regions[i].second;++j)result[j]=current[i];
      return result;
    }
  }
  throw std::runtime_error("burn_and_transport: implicit chemical solve did not converge");
}
} // namespace ember

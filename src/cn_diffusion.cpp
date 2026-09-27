#include "ember/cn_transport.hpp"
#include "cn_source.hpp"
#include "flux_chain.hpp"
#include "parallel_evaluate.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace ember {
namespace {
using V=CNSpeciesVector;using M=CNSpeciesMatrix;
using detail::independent_evaluations;
double norm(const V& v) {double n=0;for(double x:v)n=std::max(n,std::abs(x));return n;}
void check(const Composition& c,const Composition& reference) {
  if(c[Species::H2]!=0)
    throw std::domain_error("CN diffusion: initial D requires the physical seven-species solve");
  if(c.cn_mass_convention!=CNMassConvention::fixed_metal_proxy)
    throw std::domain_error("CN diffusion: explicit metal mass requires the six-species solve");
  if(!c.cn_molality)throw std::domain_error("CN diffusion: missing catalysts");
  (void)cn_physical_ledger(c,*c.cn_molality);
  for(std::size_t j=3;j<NSPEC;++j)if(c.X[j]!=reference.X[j])
    throw std::domain_error("CN diffusion: fixed uniform material lookup required");
}
}
CNSpeciesVector cn_transport_abundances(const Composition& c) {
  if(!c.cn_molality)throw std::invalid_argument("CN diffusion: missing catalysts");
  const auto& y=*c.cn_molality;return {c.X[0],c.X[1],12*y[0],13*y[1],14*y[2]};
}
Composition cn_transport_composition(const Composition& reference,const CNSpeciesVector& v) {
  auto c=reference;c.X[0]=v[0];c.X[1]=v[1];c.X[2]=1-c.Z()-v[0]-v[1];
  c.cn_molality=CNAbundances{v[2]/12,v[3]/13,v[4]/14};return c;
}
CNSpeciesFaceResponse trace_cn_flux(const SpeciesFaceResponse& hhe,
    const Composition& left,const Composition& right,CNMicroscopicApproximation mode,bool derivatives) {
  if(mode==CNMicroscopicApproximation::unselected)
    throw std::invalid_argument("CN diffusion: select a trace transport approximation explicitly");
  check(left,left);check(right,left);CNSpeciesFaceResponse out;
  for(std::size_t row=0;row<2;++row) {
    out.rate[row]=hhe.rate[row];
    if(derivatives)for(std::size_t col=0;col<2;++col) {
      out.dleft[row][col]=hhe.dleft[row][col];out.dright[row][col]=hhe.dright[row][col];
    }
  }
  if(mode==CNMicroscopicApproximation::zero_catalyst_drift)return out;
  if(mode!=CNMicroscopicApproximation::helium_velocity)
    throw std::invalid_argument("CN diffusion: unknown approximation");
  const auto a=cn_transport_abundances(left),b=cn_transport_abundances(right);
  const auto initial=initial_gs98_cn(left);
  const double original_cn=12*initial[0]+13*initial[1]+14*initial[2];
  // Physical He4 + CN = lookup He4 + original CN. Use that identity
  // directly; subtracting and adding the transported CN introduces spurious
  // roundoff dependence on otherwise independent catalyst coordinates.
  const double group_a=left.X[2]+original_cn,group_b=right.X[2]+original_cn;
  if(!(group_a>0 && group_b>0))throw std::domain_error("CN diffusion: empty helium/catalyst group");
  const double f=-(hhe.rate[0]+hhe.rate[1]),positive=std::max(f,0.),negative=std::min(f,0.);
  for(std::size_t k=2;k<5;++k) {
    const double qa=a[k]/group_a,qb=b[k]/group_b;
    // Upwind helium-group composition. At zero flux use the centered
    // generalized derivative of the continuous, piecewise-smooth flux.
    const double donor=f>0?qa:(f<0?qb:.5*(qa+qb));
    out.rate[k]=positive*qa+negative*qb;
    if(derivatives)for(std::size_t j=0;j<5;++j) {
      const double dfa=j<2?-(hhe.dleft[0][j]+hhe.dleft[1][j]):0;
      const double dfb=j<2?-(hhe.dright[0][j]+hhe.dright[1][j]):0;
      // group = 1-inert_metals-H-He3, independent of individual CN fractions.
      out.dleft[k][j]=dfa*donor+positive*(j<2?qa/group_a:(j==k?1/group_a:0));
      out.dright[k][j]=dfb*donor+negative*(j<2?qb/group_b:(j==k?1/group_b:0));
    }
  }
  return out;
}

CNTransportResult burn_cn_and_diffuse(const Model& thermal,const Model& previous,
    const PPCNNetwork& nuclear,const MixingRegions& regions,const CNSpeciesFlux& flux,double dt,
    const SpeciesTransportOptions& options,std::span<const double> common_mixing_rates) {
  if(thermal.m!=previous.m || thermal.comp.size()!=thermal.size() || previous.comp.size()!=thermal.size()
      || thermal.M!=previous.M || !flux || !std::isfinite(dt) || dt<=0
      || !std::isfinite(options.abundance_tolerance) || options.abundance_tolerance<=0
      || !std::isfinite(options.integrated_balance_tolerance) || options.integrated_balance_tolerance<0
      || !options.evaluation_threads || options.evaluation_threads>64
      || !options.max_iterations || !options.max_backtracks || regions.empty())
    throw std::invalid_argument("CN diffusion: invalid model, flux, step or options");
  const double balance_tolerance=options.integrated_balance_tolerance>0?
      options.integrated_balance_tolerance:options.abundance_tolerance;
  const double balance_scale=options.abundance_tolerance/balance_tolerance;
  if(!std::isfinite(balance_scale) || balance_scale<=0)
    throw std::invalid_argument("CN diffusion: unrepresentable tolerance ratio");
  if(!common_mixing_rates.empty() && common_mixing_rates.size()+1!=thermal.size())
    throw std::invalid_argument("CN diffusion: common mixing size differs");
  for(double g:common_mixing_rates)if(!std::isfinite(g) || g<0)
    throw std::invalid_argument("CN diffusion: invalid common mixing coefficient");
  const bool linear_mixing=std::any_of(common_mixing_rates.begin(),common_mixing_rates.end(),[](double g){return g>0;});
  const auto weights=nodal_mass_weights(thermal);const auto& reference=previous.comp.front();
  for(const auto& c:previous.comp)check(c,reference);
  std::vector<V> old;std::vector<double> mass;std::vector<std::size_t> faces;std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())throw std::invalid_argument("CN diffusion: invalid partition");
    next=end;V v{};double w=0;
    for(std::size_t i=begin;i<end;++i) {
      w+=weights[i];const auto x=cn_transport_abundances(previous.comp[i]);
      for(std::size_t k=0;k<5;++k)v[k]+=weights[i]*x[k];
    }
    for(double& x:v)x/=w;old.push_back(v);mass.push_back(w/thermal.M);
    if(end<thermal.size())faces.push_back(end-1);
  }
  if(next!=thermal.size())throw std::invalid_argument("CN diffusion: incomplete partition");
  const auto n=regions.size();auto current=old;
  if(options.initial_guess.empty()) {
    const auto predictor=burn_and_mix(thermal,previous,nuclear,regions,dt,options.abundance_tolerance*.1);
    for(std::size_t i=0;i<n;++i)current[i]=cn_transport_abundances(predictor[regions[i].first]);
  } else {
    if(options.initial_guess.size()!=thermal.size())throw std::invalid_argument("CN diffusion: guess size differs");
    for(const auto& c:options.initial_guess)check(c,reference);
    for(std::size_t i=0;i<n;++i) {
      current[i]={};
      for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell) {
        const auto v=cn_transport_abundances(options.initial_guess[cell]);
        for(std::size_t k=0;k<5;++k)current[i][k]+=weights[cell]/thermal.M/mass[i]*v[k];
      }
    }
  }
  if(options.seed_present_species) {
    V mean{};std::array<bool,5> absent{};
    for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<5;++k){mean[k]+=mass[i]*current[i][k];absent[k]|=current[i][k]==0;}
    for(auto& v:current)for(std::size_t k=0;k<5;++k)
      if(absent[k] && mean[k]>0)v[k]=.99*v[k]+.01*mean[k];
  }
  struct Evaluation {
    std::vector<M> E,A,B;std::vector<V> residual,storage;std::vector<CNSpeciesFaceResponse> faces;
    V balance{};double norm{};
  };
  const double factor=dt/thermal.M;
  auto evaluate=[&](const std::vector<V>& candidate,bool derivatives) {
    Evaluation out;out.E.resize(n);out.A.resize(n-1);out.B.resize(n-1);out.residual.resize(n);out.storage.resize(n);out.faces.resize(n-1);
    std::vector<Composition> composition(n);
    independent_evaluations(n,options.evaluation_threads,[&](std::size_t i) {
      auto c=cn_transport_composition(reference,candidate[i]);check(c,reference);composition[i]=c;
      for(std::size_t k=0;k<5;++k){out.E[i][k][k]=mass[i];out.residual[i][k]=mass[i]*(candidate[i][k]-old[i][k]);}
      for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell) {
        const auto source=detail::cn_full_source(thermal.T(cell),thermal.rho(cell),c,nuclear.pp(),nuclear.cn());
        for(std::size_t row=0;row<5;++row) {
          out.residual[i][row]-=factor*weights[cell]*source.source[row];
          if(derivatives)for(std::size_t col=0;col<5;++col)
            out.E[i][row][col]-=factor*weights[cell]*source.jacobian[row][col];
        }
      }
      out.storage[i]=out.residual[i];
    });
    for(std::size_t i=0;i<n;++i) {
      // Accumulate source/storage before adding fluxes: closed-end conservation
      // must not be hidden by cancellation of large opposing face terms.
      for(std::size_t k=0;k<5;++k)out.balance[k]+=out.residual[i][k];
    }
    independent_evaluations(n-1,options.evaluation_threads,[&](std::size_t i) {
      auto& f=out.faces[i];f=flux(faces[i],composition[i],composition[i+1],derivatives);
      detail::flux_finite(f.rate);
      if(derivatives){detail::flux_finite(f.dleft);detail::flux_finite(f.dright);}
    });
    for(std::size_t i=0;i+1<n;++i) {
      const auto& f=out.faces[i];
      const double g=linear_mixing?common_mixing_rates[faces[i]]:0.;
      for(std::size_t row=0;row<5;++row) {
        const double rate=f.rate[row]+g*(candidate[i][row]-candidate[i+1][row]);
        out.residual[i][row]+=factor*rate;out.residual[i+1][row]-=factor*rate;
        if(derivatives)for(std::size_t col=0;col<5;++col) {
          out.A[i][row][col]=factor*(f.dleft[row][col]+(row==col?g:0));
          out.B[i][row][col]=factor*(-f.dright[row][col]+(row==col?g:0));
        }
      }
    }
    for(std::size_t i=0;i<n;++i){detail::flux_finite(out.residual[i]);out.norm=std::max(out.norm,norm(out.residual[i])/mass[i]);}
    detail::flux_finite(out.balance);return out;
  };
  auto newton_correction=[&](const std::vector<V>& at,const Evaluation& value,const Evaluation& jacobian) {
    auto rhs=value.residual;std::vector<V> offset(n-1);
    if(!linear_mixing) {
      for(auto& row:rhs)for(double& x:row)x=-x;
      return detail::solve_flux_chain<5>(jacobian.E,jacobian.A,jacobian.B,rhs,offset);
    }
    // Solve for the absolute next state. The linear mixing flux has exactly
    // zero affine offset; forming it by subtracting its large face residuals
    // would lose finite mass/reaction terms before block elimination begins.
    for(std::size_t i=0;i<n;++i) {
      rhs[i]=detail::flux_product(jacobian.E[i],at[i]);
      for(std::size_t k=0;k<5;++k)rhs[i][k]-=value.storage[i][k];
    }
    for(std::size_t i=0;i+1<n;++i) {
      const auto left=detail::flux_product(jacobian.faces[i].dleft,at[i]);
      const auto right=detail::flux_product(jacobian.faces[i].dright,at[i+1]);
      for(std::size_t k=0;k<5;++k)offset[i][k]=factor*(value.faces[i].rate[k]-left[k]-right[k]);
    }
    auto destination=detail::solve_flux_chain<5>(jacobian.E,jacobian.A,jacobian.B,rhs,offset);
    for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<5;++k)destination[i][k]-=at[i][k];
    return destination;
  };
  CNTransportResult result;std::string last_rejection;
  for(std::size_t iteration=0;iteration<options.max_iterations;++iteration) {
    const auto state=evaluate(current,true);result.residual_history.push_back(state.norm);
    const auto correction=newton_correction(current,state,state);
    double change=0;for(const auto& row:correction)change=std::max(change,norm(row));
    result.residual=state.norm;result.abundance_correction=change;
    const double merit=std::max(change,balance_scale*norm(state.balance));
    if(merit<=options.abundance_tolerance) {
      result.iterations=iteration;result.integrated_balance=state.balance;result.composition=previous.comp;
      for(std::size_t i=0;i<n;++i)for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell)
        result.composition[cell]=cn_transport_composition(reference,current[i]);
      for(std::size_t i=0;i+1<n;++i) {
        auto rate=state.faces[i].rate;const double g=linear_mixing?common_mixing_rates[faces[i]]:0.;
        for(std::size_t k=0;k<5;++k)rate[k]+=g*(current[i][k]-current[i+1][k]);
        result.boundary_fluxes.push_back({faces[i],rate});
      }
      return result;
    }
    bool accepted=false;double best_merit=merit,accepted_damping=0;std::vector<V> best;
    auto search=[&](const std::vector<V>& direction) {
      double damping=1;std::vector<V> last;
      for(std::size_t trial=0;trial<options.max_backtracks;++trial,damping*=.5) {
        auto candidate=current;
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<5;++k)candidate[i][k]+=damping*direction[i][k];
        if(candidate==current)break;if(candidate==last)continue;last=candidate;
        try {
          const auto test=evaluate(candidate,false);
          const auto remaining=newton_correction(candidate,test,state);
          double trial_merit=balance_scale*norm(test.balance);for(const auto& row:remaining)trial_merit=std::max(trial_merit,norm(row));
          if(trial_merit<merit*(1-1e-4*damping) || trial_merit<=options.abundance_tolerance) {
            if(!accepted || trial_merit<best_merit){best=std::move(candidate);best_merit=trial_merit;accepted_damping=damping;accepted=true;}
            break;
          }
        } catch(const std::domain_error& e){last_rejection=e.what();}
      }
    };
    search(correction);
    if(!accepted || accepted_damping<.125) {
      auto limited=correction;bool changed=false;
      for(std::size_t i=0;i<n;++i) {
        double positive=0,negative=0;
        for(std::size_t k=0;k<5;++k) {
          if(limited[i][k]<-current[i][k]){limited[i][k]=-.9*current[i][k];changed=true;}
          if(limited[i][k]>0)positive+=limited[i][k];else negative+=limited[i][k];
        }
        const auto c=cn_transport_composition(reference,current[i]);const double helium=cn_physical_ledger(c,*c.cn_molality).helium4;
        if(positive+negative>helium && positive>0) {
          const double fraction=(.9*helium-negative)/positive;
          for(double& x:limited[i])if(x>0)x*=fraction;changed=true;
        }
      }
      if(changed)search(limited);
    }
    if(!accepted) {
      std::ostringstream error;error<<std::setprecision(4)<<"CN diffusion: line search failed at iteration "<<iteration
        <<", correction "<<change<<", balance "<<norm(state.balance)<<", last rejection "<<last_rejection;
      throw std::runtime_error(error.str());
    }
    current=std::move(best);
  }
  throw std::runtime_error("CN diffusion: iteration limit");
}
} // namespace ember

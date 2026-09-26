#include "ember/species_transport.hpp"
#include "flux_chain.hpp"
#include "ember/nuclear_cn.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <sstream>

namespace ember {
namespace {
using V=SpeciesVector;using M=SpeciesMatrix;
void finite(const V& a){for(double x:a)if(!std::isfinite(x))throw std::domain_error("species transport: non-finite vector");}
void finite(const M& a){for(const auto& row:a)finite(row);}
void composition_check(const Composition& c,const Composition& reference,bool allow_cn=false) {
  if(c.cn_molality && !allow_cn)
    throw std::invalid_argument("species transport: explicit CN transport not connected");
  if(allow_cn && c.cn_molality.has_value()!=reference.cn_molality.has_value())
    throw std::invalid_argument("species reconstruction: inconsistent CN activation");
  if(c.cn_molality)(void)cn_physical_ledger(c,*c.cn_molality);
  if(c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=reference.metal_inventory
      || !std::isfinite(c.sum()) || std::abs(c.sum()-1)>1e-10)
    throw std::domain_error("species transport: normalized common baryonic composition required");
  for(std::size_t j=0;j<NSPEC;++j) {
    if(!std::isfinite(c.X[j]) || c.X[j]<0)throw std::domain_error("species transport: negative or non-finite abundance");
    if(j>=3 && std::abs(c.X[j]-reference.X[j])>1e-12)
      throw std::domain_error("species transport: metal gradients require additional transported species");
  }
}
void reaction_check(const NuclearResponse& r,bool derivatives) {
  double sum=0,scale=0;
  for(std::size_t j=0;j<NSPEC;++j) {
    const double v=r.state.dXdt[j];
    if(!std::isfinite(v))throw std::domain_error("burn and diffuse: non-finite reaction source");
    if(j>=3 && v!=0)throw std::logic_error("burn and diffuse: reaction changes an unrepresented species");
    sum+=v;scale+=std::abs(v);
  }
  if(std::abs(sum)>1e-12*scale)throw std::logic_error("burn and diffuse: reaction does not conserve baryonic mass");
  if(derivatives)for(std::size_t k=0;k<3;++k) {
    sum=scale=0;
    for(std::size_t j=0;j<NSPEC;++j) {
      const double v=r.d_dXdt_dX[j][k];
      if(!std::isfinite(v))throw std::domain_error("burn and diffuse: non-finite reaction derivative");
      if(j>=3 && v!=0)throw std::logic_error("burn and diffuse: derivative changes an unrepresented species");
      sum+=v;scale+=std::abs(v);
    }
    if(std::abs(sum)>1e-12*scale)throw std::logic_error("burn and diffuse: reaction derivative does not conserve baryonic mass");
  }
}
}

std::vector<V> solve_species_flux_chain(std::span<const M> E,std::span<const M> A,
    std::span<const M> B,std::span<const V> b,std::span<const V> d) {
  return detail::solve_flux_chain<2>(E,A,B,b,d);
}

SpeciesFluxReconstruction reconstruct_species_fluxes(const Model& current,
    const Model& previous,const Nuclear& nuclear,const MixingRegions& regions,
    std::span<const SpeciesBoundaryFlux> boundaries,double dt) {
  if(current.m!=previous.m || current.size()!=previous.size()
      || current.comp.size()!=current.size() || previous.comp.size()!=previous.size()
      || current.M!=previous.M || !(current.M>0) || !std::isfinite(current.M)
      || !(dt>0) || !std::isfinite(dt) || regions.empty()
      || boundaries.size()!=regions.size()-1)
    throw std::invalid_argument("species reconstruction: invalid models, regions, boundaries or time");
  const auto weights=nodal_mass_weights(current);
  for(const auto& c:current.comp)composition_check(c,previous.comp.front(),true);
  for(const auto& c:previous.comp)composition_check(c,previous.comp.front(),true);
  std::size_t next=0;
  for(std::size_t r=0;r<regions.size();++r) {
    const auto [begin,end]=regions[r];
    if(begin!=next || end<=begin || end>current.size())
      throw std::invalid_argument("species reconstruction: invalid region partition");
    next=end;
    for(std::size_t i=begin+1;i<end;++i)
      if(current.comp[i]!=current.comp[begin])
        throw std::invalid_argument("species reconstruction: nonuniform new mixed region");
    if(r+1<regions.size()) {
      if(boundaries[r].face!=end-1)
        throw std::invalid_argument("species reconstruction: boundary index differs from partition");
      finite(boundaries[r].rate);
    }
  }
  if(next!=current.size())throw std::invalid_argument("species reconstruction: incomplete region partition");
  SpeciesFluxReconstruction result;
  result.face_rates.resize(current.size()-1);result.cell_balances.resize(current.size());
  result.region_balances.resize(regions.size());
  std::vector<std::array<long double,2>> change(current.size());
  std::array<long double,2> integrated{};
  const long double factor=static_cast<long double>(dt)/current.M;
  for(std::size_t r=0;r<regions.size();++r) {
    const auto [begin,end]=regions[r];
    const V left=r?boundaries[r-1].rate:V{};
    const V right=r+1<regions.size()?boundaries[r].rate:V{};
    std::array<long double,2> accumulated{};
    for(std::size_t i=begin;i<end;++i) {
      NuclearResponse reaction;reaction.state=nuclear.eval(current.T(i),current.rho(i),current.comp[i]);
      reaction_check(reaction,false);
      for(std::size_t k=0;k<2;++k) {
        change[i][k]=(static_cast<long double>(weights[i])/current.M)*
          ((static_cast<long double>(current.comp[i].X[k])-previous.comp[i].X[k])
           -static_cast<long double>(dt)*reaction.state.dXdt[k]);
        accumulated[k]+=change[i][k];integrated[k]+=change[i][k];
        if(i+1<end)result.face_rates[i][k]=static_cast<double>(left[k]-accumulated[k]/factor);
      }
      if(i+1<end)finite(result.face_rates[i]);
    }
    if(end<current.size())result.face_rates[end-1]=right;
    for(std::size_t k=0;k<2;++k)
      result.region_balances[r][k]=static_cast<double>(accumulated[k]+factor*(right[k]-left[k]));
    finite(result.region_balances[r]);
  }
  for(std::size_t i=0;i<current.size();++i) {
    const V left=i?result.face_rates[i-1]:V{};
    const V right=i+1<current.size()?result.face_rates[i]:V{};
    for(std::size_t k=0;k<2;++k)
      result.cell_balances[i][k]=static_cast<double>(change[i][k]+factor*
        (static_cast<long double>(right[k])-left[k]));
    finite(result.cell_balances[i]);
  }
  for(std::size_t k=0;k<2;++k)result.integrated_balance[k]=static_cast<double>(integrated[k]);
  finite(result.integrated_balance);return result;
}

SpeciesTransportResult burn_and_diffuse(const Model& thermal,const Model& previous,
    const Nuclear& nuclear,const MixingRegions& regions,const SpeciesFlux& flux,double dt,
    const SpeciesTransportOptions& options) {
  if(thermal.m!=previous.m || thermal.size()!=previous.size() || thermal.comp.size()!=thermal.size()
      || previous.comp.size()!=previous.size() || !flux || !(dt>0) || !std::isfinite(dt)
      || !(options.abundance_tolerance>0) || !std::isfinite(options.abundance_tolerance)
      || !options.max_iterations || !options.max_backtracks || regions.empty()
      || thermal.M!=previous.M || !(thermal.M>0) || !std::isfinite(thermal.M))
    throw std::invalid_argument("burn and diffuse: invalid models, callback, time or options");
  const auto weights=nodal_mass_weights(thermal);
  for(const auto& c:previous.comp)composition_check(c,previous.comp.front());
  std::vector<Composition> old;std::vector<double> mass;std::vector<std::size_t> faces;
  std::size_t next=0;
  for(auto [begin,end]:regions) {
    if(begin!=next || end<=begin || end>thermal.size())throw std::invalid_argument("burn and diffuse: invalid region partition");
    next=end;Composition c=previous.comp[begin];c.X.fill(0);double weight=0;
    for(std::size_t i=begin;i<end;++i) {
      weight+=weights[i];for(std::size_t j=0;j<NSPEC;++j)c.X[j]+=weights[i]*previous.comp[i].X[j];
    }
    for(auto& x:c.X)x/=weight;
    c.X[2]=1-c.Z()-c.X[0]-c.X[1];old.push_back(c);mass.push_back(weight/thermal.M);
    if(end<thermal.size())faces.push_back(end-1);
  }
  if(next!=thermal.size())throw std::invalid_argument("burn and diffuse: incomplete region partition");
  auto current=old;
  if(options.initial_guess.empty()) {
    const auto predictor=burn_and_mix(thermal,previous,nuclear,regions,dt,options.abundance_tolerance*.1);
    for(std::size_t i=0;i<regions.size();++i) {
      current[i]=predictor[regions[i].first];
      // These carriers are fixed in the represented problem. Region averaging
      // must not introduce one-ulp differences interpreted as metal gradients.
      for(std::size_t j=3;j<NSPEC;++j)current[i].X[j]=previous.comp.front().X[j];
      current[i].X[2]=1-current[i].Z()-current[i].X[0]-current[i].X[1];
    }
  } else {
    if(options.initial_guess.size()!=thermal.size())throw std::invalid_argument("burn and diffuse: initial guess size differs");
    for(const auto& c:options.initial_guess)composition_check(c,previous.comp.front());
    for(std::size_t r=0;r<regions.size();++r) {
      current[r].X.fill(0);double weight=0;
      for(std::size_t i=regions[r].first;i<regions[r].second;++i) {
        weight+=weights[i];
        for(std::size_t j=0;j<NSPEC;++j)current[r].X[j]+=weights[i]*options.initial_guess[i].X[j];
      }
      for(auto& x:current[r].X)x/=weight;
      // The conserved metal slots are fixed, rather than changed by rounding
      // while averaging a guess. Only the two represented species can evolve.
      for(std::size_t j=3;j<NSPEC;++j)current[r].X[j]=previous.comp.front().X[j];
      current[r].X[2]=1-current[r].Z()-current[r].X[0]-current[r].X[1];
    }
  }
  bool interior_adjusted=false;
  if(options.seed_present_species) {
    std::array<long double,2> sum{};long double total=0;
    std::array<bool,2> zero{};
    for(std::size_t i=0;i<current.size();++i) {
      total+=mass[i];
      for(std::size_t k=0;k<2;++k){sum[k]+=static_cast<long double>(mass[i])*current[i].X[k];zero[k]|=current[i].X[k]==0;}
    }
    interior_adjusted=(zero[0] && sum[0]>0)||(zero[1] && sum[1]>0);
    if(interior_adjusted)for(auto& c:current) {
      // A numerical starting guess, never a source of material or a minimum
      // accepted abundance. The two original conservation equations remain.
      for(std::size_t k=0;k<2;++k)c.X[k]=.99*c.X[k]+.01*static_cast<double>(sum[k]/total);
      c.X[2]=1-c.Z()-c.X[0]-c.X[1];
    }
  }
  struct Evaluation {std::vector<M> E,A,B;std::vector<V> residual;std::vector<SpeciesFaceResponse> faces;double norm{};};
  const auto n=regions.size();const double factor=dt/thermal.M;
  auto evaluate=[&](const std::vector<Composition>& candidate,bool derivatives) {
    Evaluation a;a.E.resize(n);a.A.resize(n-1);a.B.resize(n-1);a.residual.resize(n);a.faces.resize(n-1);
    for(std::size_t i=0;i<n;++i) {
      composition_check(candidate[i],previous.comp.front());
      a.E[i][0][0]=a.E[i][1][1]=mass[i];
      for(std::size_t j=0;j<2;++j)a.residual[i][j]=mass[i]*(candidate[i].X[j]-old[i].X[j]);
      for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell) {
        NuclearResponse response;
        if(derivatives)response=nuclear.composition_response(thermal.T(cell),thermal.rho(cell),candidate[i]);
        else response.state=nuclear.eval(thermal.T(cell),thermal.rho(cell),candidate[i]);
        reaction_check(response,derivatives);
        for(std::size_t j=0;j<2;++j) {
          const double source=factor*weights[cell]*response.state.dXdt[j];
          a.residual[i][j]-=source;
          if(derivatives)for(std::size_t k=0;k<2;++k) {
            const double jac=factor*weights[cell]*(response.d_dXdt_dX[j][k]-response.d_dXdt_dX[j][2]);
            a.E[i][j][k]-=jac;
          }
        }
      }
    }
    for(std::size_t i=0;i+1<n;++i) {
      a.faces[i]=flux(faces[i],candidate[i],candidate[i+1],derivatives);
      finite(a.faces[i].rate);
      if(derivatives){finite(a.faces[i].dleft);finite(a.faces[i].dright);}
      for(std::size_t j=0;j<2;++j) {
        const double f=factor*a.faces[i].rate[j];a.residual[i][j]+=f;a.residual[i+1][j]-=f;
        if(derivatives)for(std::size_t k=0;k<2;++k) {
          a.A[i][j][k]=factor*a.faces[i].dleft[j][k];a.B[i][j][k]=-factor*a.faces[i].dright[j][k];
        }
      }
    }
    for(std::size_t i=0;i<n;++i) {
      finite(a.residual[i]);for(double x:a.residual[i])a.norm=std::max(a.norm,std::abs(x)/mass[i]);
    }
    if(!std::isfinite(a.norm))throw std::domain_error("burn and diffuse: non-finite residual");
    return a;
  };
  SpeciesTransportResult result;result.interior_guess_adjusted=interior_adjusted;
  std::string last_domain_error;
  for(std::size_t iteration=0;iteration<options.max_iterations;++iteration) {
    auto state=evaluate(current,true);result.residual_history.push_back(state.norm);
    auto rhs=state.residual;
    for(auto& row:rhs)for(auto& x:row)x=-x;
    const auto correction=solve_species_flux_chain(state.E,state.A,state.B,rhs,std::vector<V>(n-1));
    double change=0;V balance{};
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<2;++j) {
      change=std::max(change,std::abs(correction[i][j]));balance[j]+=state.residual[i][j];
    }
    result.abundance_correction=change;result.residual=state.norm;
    if(change<=options.abundance_tolerance
        && std::max(std::abs(balance[0]),std::abs(balance[1]))<=options.abundance_tolerance) {
      result.iterations=iteration;result.residual=state.norm;result.abundance_correction=change;result.composition=previous.comp;
      for(std::size_t i=0;i<n;++i) {
        for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell)result.composition[cell]=current[i];
        for(std::size_t j=0;j<2;++j)result.integrated_balance[j]+=state.residual[i][j];
      }
      for(std::size_t i=0;i+1<n;++i)result.boundary_fluxes.push_back({faces[i],state.faces[i].rate});
      return result;
    }
    // Use the same quantities to judge a trial and convergence. A stiff
    // tiny-cell rate residual can increase while the abundance error falls.
    // The current Jacobian converts each trial residual to an estimated
    // remaining abundance correction; the next iteration recomputes it.
    const double merit=std::max({change,std::abs(balance[0]),std::abs(balance[1])});
    bool accepted=false;double best_merit=merit,accepted_damping=0;
    std::vector<Composition> best;
    auto search=[&](const std::vector<V>& direction) {
      double damping=1;
      std::vector<Composition> last_candidate;
      const auto same=[](const std::vector<Composition>& a,const std::vector<Composition>& b) {
        if(a.size()!=b.size())return false;
        for(std::size_t i=0;i<a.size();++i)if(a[i].X!=b[i].X)return false;
        return true;
      };
      for(std::size_t trial=0;trial<options.max_backtracks;++trial,damping*=.5) {
        auto candidate=current;
        for(std::size_t i=0;i<n;++i) {
          for(std::size_t j=0;j<2;++j)candidate[i].X[j]+=damping*direction[i][j];
          candidate[i].X[2]=1-candidate[i].Z()-candidate[i].X[0]-candidate[i].X[1];
        }
        // Smaller corrections cannot change any abundance after this point.
        // Scalar/derivative rounding differences must not accept a stationary
        // iterate as progress. Also avoid evaluating an identical trial twice.
        if(same(candidate,current))break;
        if(same(candidate,last_candidate))continue;
        last_candidate=candidate;
        try {
          const auto test=evaluate(candidate,false);
          // Reuse the flux evaluation. Only the small block solve is added.
          auto trial_rhs=test.residual;
          for(auto& row:trial_rhs)for(auto& x:row)x=-x;
          const auto remaining=solve_species_flux_chain(state.E,state.A,state.B,
              trial_rhs,std::vector<V>(n-1));
          double remaining_norm=0;
          for(const auto& row:remaining)for(double x:row)
            remaining_norm=std::max(remaining_norm,std::abs(x));
          V trial_balance{};
          for(const auto& row:test.residual)for(std::size_t j=0;j<2;++j)trial_balance[j]+=row[j];
          const double trial_merit=std::max({remaining_norm,
              std::abs(trial_balance[0]),std::abs(trial_balance[1])});
          if(trial_merit<merit*(1-1e-4*damping) || trial_merit<=options.abundance_tolerance) {
            if(!accepted || trial_merit<best_merit) {
              best=std::move(candidate);best_merit=trial_merit;accepted_damping=damping;accepted=true;
            }
            break;
          }
        } catch(const std::domain_error& e) {last_domain_error=e.what();}
      }
    };
    search(correction);
    if(!accepted || accepted_damping<.125) {
      // One nearly depleted component can otherwise stop every other Newton
      // component. Try a bounded direction with the same correction/balance
      // merit and the FULL new Newton correction for convergence.
      // This is a step limiter, not clipping of an accepted composition.
      auto limited=correction;bool changed=false;
      for(std::size_t i=0;i<n;++i) {
        double positive=0,negative=0;
        for(std::size_t j=0;j<2;++j) {
          if(limited[i][j]<-current[i].X[j]) {limited[i][j]=-.9*current[i].X[j];changed=true;}
          if(limited[i][j]>0)positive+=limited[i][j];else negative+=limited[i][j];
        }
        if(positive+negative>current[i].X[2]) {
          const double scale=(.9*current[i].X[2]-negative)/positive;
          for(auto& value:limited[i])if(value>0)value*=scale;
          changed=true;
        }
      }
      if(changed)search(limited);
    }
    if(!accepted) {
      std::ostringstream message;message<<std::setprecision(4)<<"burn and diffuse: line search failed at iteration "
        <<iteration<<", residual "<<state.norm<<", Newton abundance correction "<<change
        <<", integrated balances "<<balance[0]<<' '<<balance[1]<<", last domain rejection: "<<last_domain_error;
      throw std::runtime_error(message.str());
    }
    current=std::move(best);
  }
  std::ostringstream message;message<<std::setprecision(4)<<"burn and diffuse: iteration limit "<<options.max_iterations
    <<", final residual "<<result.residual<<", Newton abundance correction "<<result.abundance_correction<<", residual history";
  for(double value:result.residual_history)message<<' '<<value;
  message<<", last domain rejection: "<<last_domain_error;
  throw std::runtime_error(message.str());
}
} // namespace ember

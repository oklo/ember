#include "ember/metal_cn_transport.hpp"
#include "metal_cn_source.hpp"
#include "flux_chain.hpp"
#include "parallel_evaluate.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace ember {
namespace {
using V=MetalCNVector;using M=MetalCNMatrix;
using detail::independent_evaluations;
void check(const Composition& c,const Composition& reference) {
  if(!c.cn_molality || c.cn_mass_convention!=CNMassConvention::explicit_metal_mass
      || c.basis!=AbundanceBasis::baryon_mass || c.metal_inventory!=MetalInventory::gs98)
    throw std::domain_error("metal CN transport: explicit physical composition required");
  (void)cn_physical_ledger(c,*c.cn_molality);
  if(c.Z()>0 && reference.Z()>0)for(std::size_t k=3;k<METAL_END;++k)
    if(std::abs(c.X[k]/c.Z()-reference.X[k]/reference.Z())>1e-12)
      throw std::domain_error("metal CN transport: material patterns differ");
}
}
MetalCNVector metal_cn_abundances(const Composition& c) {
  check(c,c);const auto& y=*c.cn_molality;
  const double carbon12=12*y[0],carbon13=13*y[1],nitrogen14=14*y[2];
  return {c.X[0],c.X[1],carbon12,carbon13,nitrogen14,c.Z()-(carbon12+carbon13+nitrogen14),c[Species::H2]};
}
Composition metal_cn_composition(const Composition& reference,const MetalCNVector& v) {
  check(reference,reference);double total=0;
  for(double x:v){if(!std::isfinite(x) || x<0)throw std::domain_error("metal CN transport: invalid independent fraction");total+=x;}
  if(total>1)throw std::domain_error("metal CN transport: negative reference helium");
  const double Z=v[2]+v[3]+v[4]+v[5];auto c=reference;
  c.X[0]=v[0];c.X[1]=v[1];c.X[2]=1-total;c[Species::H2]=v[METAL_CN_D];
  const auto pattern=reference.Z()>0?reference:solar_scaled(0.,1.);
  for(std::size_t k=3;k<METAL_END;++k)c.X[k]=Z*pattern.X[k]/pattern.Z();
  c.cn_molality=CNAbundances{v[2]/12,v[3]/13,v[4]/14};check(c,reference);return c;
}
MetalCNFaceResponse common_metal_cn_flux(const MetalSpeciesFaceResponse& flux,
    const Composition& left,const Composition& right,bool derivatives) {
  if(left[Species::H2]!=0 || right[Species::H2]!=0)
    throw std::domain_error("metal CN flux: microscopic initial-D transport requires an isotope-aware provider");
  const auto a=metal_cn_abundances(left),b=metal_cn_abundances(right);check(right,left);
  MetalCNFaceResponse out;
  for(std::size_t row=0;row<2;++row) {
    out.rate[row]=flux.rate[row];
    if(derivatives)for(std::size_t col=0;col<6;++col) {
      const auto material=col<2?col:2;
      out.dleft[row][col]=flux.dleft[row][material];out.dright[row][col]=flux.dright[row][material];
    }
  }
  const double Za=left.Z(),Zb=right.Z(),f=flux.rate[2];
  if(Za<=0 || Zb<=0) {
    if(Za==0 && Zb==0 && f==0)return out;
    throw std::domain_error("metal CN flux: positive endpoint metal groups required");
  }
  const double positive=std::max(f,0.),negative=std::min(f,0.);
  for(std::size_t row=2;row<6;++row) {
    const double qa=a[row]/Za,qb=b[row]/Zb,donor=f>0?qa:(f<0?qb:.5*(qa+qb));
    out.rate[row]=positive*qa+negative*qb;
    if(derivatives)for(std::size_t col=0;col<6;++col) {
      const auto material=col<2?col:2;
      out.dleft[row][col]=flux.dleft[2][material]*donor
        +positive*((row==col?1.:0.)-(col>=2?qa:0.))/Za;
      out.dright[row][col]=flux.dright[2][material]*donor
        +negative*((row==col?1.:0.)-(col>=2?qb:0.))/Zb;
    }
  }
  return out;
}

MetalCNTransportResult burn_metal_cn_and_diffuse(const Model& thermal,const Model& previous,
    const PPCNNetwork& nuclear,const MixingRegions& regions,const MetalCNFlux& flux,double dt,
    const SpeciesTransportOptions& options,std::span<const double> common_mixing_rates) {
  if(thermal.m!=previous.m || thermal.comp.size()!=thermal.size() || previous.comp.size()!=thermal.size()
      || thermal.M!=previous.M || !flux || !std::isfinite(dt) || dt<=0
      || !std::isfinite(options.abundance_tolerance) || options.abundance_tolerance<=0
      || !options.evaluation_threads || options.evaluation_threads>64
      || !options.max_iterations || !options.max_backtracks || regions.empty())
    throw std::invalid_argument("CN diffusion: invalid model, flux, step or options");
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
    next=end;const auto anchor=metal_cn_abundances(previous.comp[begin]);V v{};double w=0;
    for(std::size_t i=begin;i<end;++i) {
      w+=weights[i];const auto x=metal_cn_abundances(previous.comp[i]);
      for(std::size_t k=0;k<METAL_CN_SIZE;++k)v[k]+=weights[i]*(x[k]-anchor[k]);
    }
    for(std::size_t k=0;k<METAL_CN_SIZE;++k)v[k]=anchor[k]+v[k]/w;old.push_back(v);mass.push_back(w/thermal.M);
    if(end<thermal.size())faces.push_back(end-1);
  }
  if(next!=thermal.size())throw std::invalid_argument("CN diffusion: incomplete partition");
  double maximum_D=0;
  for(const auto& v:old)maximum_D=std::max(maximum_D,v[METAL_CN_D]);
  const Nuclear& light=maximum_D>0?static_cast<const Nuclear&>(nuclear.light()):nuclear.pp();
  const double dunit=maximum_D>0 && maximum_D<1e-100?maximum_D:1.;
  const double dtol=dunit!=1.?1e-10:std::min(options.abundance_tolerance,
      std::max(1e-10*maximum_D,16*std::numeric_limits<double>::denorm_min()));
  auto physical=[&](V v){v[METAL_CN_D]*=dunit;return v;};
  auto composition=[&](const V& v){return metal_cn_composition(reference,physical(v));};
  auto measure=[&](const V& v) {
    double error=0;
    for(std::size_t k=0;k<METAL_CN_SIZE;++k)
      error=std::max(error,k==METAL_CN_D?std::abs(v[k])/dtol*options.abundance_tolerance:std::abs(v[k]));
    return error;
  };
  for(auto& v:old)v[METAL_CN_D]/=dunit;
  const auto n=regions.size();auto current=old;
  if(options.initial_guess.empty()) {
    // Old regional abundances are a conservative starting iterate. The
    // implicit solve updates burning and all seven independent masses together.
  } else {
    if(options.initial_guess.size()!=thermal.size())throw std::invalid_argument("CN diffusion: guess size differs");
    for(const auto& c:options.initial_guess)check(c,reference);
    for(std::size_t i=0;i<n;++i) {
      current[i]={};
      for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell) {
        const auto v=metal_cn_abundances(options.initial_guess[cell]);
        for(std::size_t k=0;k<METAL_CN_SIZE;++k)current[i][k]+=weights[cell]/thermal.M/mass[i]*v[k];
      }
    }
    for(auto& v:current)v[METAL_CN_D]/=dunit;
  }
  if(options.seed_present_species) {
    V mean{};std::array<bool,METAL_CN_SIZE> absent{};
    for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<METAL_CN_SIZE;++k){mean[k]+=mass[i]*current[i][k];absent[k]|=current[i][k]==0;}
    for(auto& v:current)for(std::size_t k=0;k<METAL_CN_SIZE;++k)
      if(absent[k] && mean[k]>0)v[k]=.99*v[k]+.01*mean[k];
  }
  struct Evaluation {
    std::vector<M> E,A,B;std::vector<V> residual,storage;std::vector<MetalCNFaceResponse> faces;
    V balance{};double norm{};
  };
  const double factor=dt/thermal.M;
  auto evaluate=[&](const std::vector<V>& candidate,bool derivatives) {
    Evaluation out;out.E.resize(n);out.A.resize(n-1);out.B.resize(n-1);out.residual.resize(n);out.storage.resize(n);out.faces.resize(n-1);
    std::vector<Composition> compositions(n);
    independent_evaluations(n,options.evaluation_threads,[&](std::size_t i) {
      auto c=composition(candidate[i]);check(c,reference);compositions[i]=c;
      for(std::size_t k=0;k<METAL_CN_SIZE;++k){out.E[i][k][k]=mass[i];out.residual[i][k]=mass[i]*(candidate[i][k]-old[i][k]);}
      for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell) {
        const auto source=detail::metal_cn_source(thermal.T(cell),thermal.rho(cell),c,light,nuclear.cn());
        for(std::size_t row=0;row<METAL_CN_SIZE;++row) {
          const double row_unit=row==METAL_CN_D?dunit:1.;
          // For trace D, lambda remains representable after lambda*XD has
          // underflowed. Its D-screening correction is below binary64
          // precision here; this is the same implicit destruction equation.
          const double rate=row==METAL_CN_D && dunit!=1.
              ?source.jacobian[row][row]*candidate[i][row]:source.source[row]/row_unit;
          out.residual[i][row]-=factor*weights[cell]*rate;
          if(derivatives)for(std::size_t col=0;col<METAL_CN_SIZE;++col) {
            const double column_unit=col==METAL_CN_D?dunit:1.;
            const double jacobian=row==METAL_CN_D && col==METAL_CN_D
                ?source.jacobian[row][col]:source.jacobian[row][col]*column_unit/row_unit;
            out.E[i][row][col]-=factor*weights[cell]*jacobian;
          }
        }
      }
      out.storage[i]=out.residual[i];
    });
    for(std::size_t i=0;i<n;++i) {
      // Accumulate source/storage before adding fluxes: closed-end conservation
      // must not be hidden by cancellation of large opposing face terms.
      for(std::size_t k=0;k<METAL_CN_SIZE;++k)out.balance[k]+=out.residual[i][k];
    }
    independent_evaluations(n-1,options.evaluation_threads,[&](std::size_t i) {
      auto& f=out.faces[i];f=flux(faces[i],compositions[i],compositions[i+1],derivatives);
      f.rate[METAL_CN_D]/=dunit;
      if(derivatives)for(std::size_t row=0;row<METAL_CN_SIZE;++row)for(std::size_t col=0;col<METAL_CN_SIZE;++col) {
        if(row==METAL_CN_D && col==METAL_CN_D)continue;
        const double column_unit=col==METAL_CN_D?dunit:1.,row_unit=row==METAL_CN_D?dunit:1.;
        f.dleft[row][col]=f.dleft[row][col]*column_unit/row_unit;
        f.dright[row][col]=f.dright[row][col]*column_unit/row_unit;
      }
      detail::flux_finite(f.rate);
      if(derivatives){detail::flux_finite(f.dleft);detail::flux_finite(f.dright);}
    });
    for(std::size_t i=0;i+1<n;++i) {
      const auto& f=out.faces[i];
      const double g=linear_mixing?common_mixing_rates[faces[i]]:0.;
      for(std::size_t row=0;row<METAL_CN_SIZE;++row) {
        const double rate=f.rate[row]+g*(candidate[i][row]-candidate[i+1][row]);
        out.residual[i][row]+=factor*rate;out.residual[i+1][row]-=factor*rate;
        if(derivatives)for(std::size_t col=0;col<METAL_CN_SIZE;++col) {
          out.A[i][row][col]=factor*(f.dleft[row][col]+(row==col?g:0));
          out.B[i][row][col]=factor*(-f.dright[row][col]+(row==col?g:0));
        }
      }
    }
    for(std::size_t i=0;i<n;++i){detail::flux_finite(out.residual[i]);out.norm=std::max(out.norm,measure(out.residual[i])/mass[i]);}
    detail::flux_finite(out.balance);return out;
  };
  auto newton_correction=[&](const std::vector<V>& at,const Evaluation& value,const Evaluation& jacobian,
      std::vector<V>* conserved_flux=nullptr) {
    auto rhs=value.residual;std::vector<V> offset(n-1);
    if(!linear_mixing) {
      for(auto& row:rhs)for(double& x:row)x=-x;
      return detail::solve_flux_chain<METAL_CN_SIZE>(jacobian.E,jacobian.A,jacobian.B,rhs,offset);
    }
    // Solve for the absolute next state. The linear mixing flux has exactly
    // zero affine offset; forming it by subtracting its large face residuals
    // would lose finite mass/reaction terms before block elimination begins.
    for(std::size_t i=0;i<n;++i) {
      rhs[i]=detail::flux_product(jacobian.E[i],at[i]);
      for(std::size_t k=0;k<METAL_CN_SIZE;++k)rhs[i][k]-=value.storage[i][k];
    }
    for(std::size_t i=0;i+1<n;++i) {
      const auto left=detail::flux_product(jacobian.faces[i].dleft,at[i]);
      const auto right=detail::flux_product(jacobian.faces[i].dright,at[i+1]);
      for(std::size_t k=0;k<METAL_CN_SIZE;++k)offset[i][k]=factor*(value.faces[i].rate[k]-left[k]-right[k]);
    }
    auto destination=detail::solve_flux_chain<METAL_CN_SIZE>(jacobian.E,jacobian.A,jacobian.B,rhs,offset,conserved_flux);
    for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<METAL_CN_SIZE;++k)destination[i][k]-=at[i][k];
    return destination;
  };
  MetalCNTransportResult result;std::string last_rejection;
  for(std::size_t iteration=0;iteration<options.max_iterations;++iteration) {
    const auto state=evaluate(current,true);result.residual_history.push_back(state.norm);
    std::vector<V> conserved_flux;
    const auto correction=newton_correction(current,state,state,linear_mixing?&conserved_flux:nullptr);
    double change=0;for(const auto& row:correction)change=std::max(change,measure(row));
    result.residual=state.norm;result.abundance_correction=change;
    const double merit=std::max(change,measure(state.balance));
    if(merit<=options.abundance_tolerance) {
      result.iterations=iteration;result.integrated_balance=physical(state.balance);result.composition=previous.comp;
      for(std::size_t i=0;i<n;++i)for(std::size_t cell=regions[i].first;cell<regions[i].second;++cell)
        result.composition[cell]=composition(current[i]);
      for(std::size_t i=0;i+1<n;++i) {
        auto rate=state.faces[i].rate;
        // Return the flux of the converged linearized conservation solve.
        // Its constituent gradients can be smaller than a stored abundance
        // ULP. The independent reconstruction still checks every cell at
        // the original tolerance; no boundary or source is repaired there.
        if(linear_mixing)for(std::size_t k=0;k<METAL_CN_SIZE;++k)rate[k]=conserved_flux[i][k]/factor;
        result.boundary_fluxes.push_back({faces[i],physical(rate)});
      }
      return result;
    }
    bool accepted=false;double best_merit=merit,accepted_damping=0;std::vector<V> best;
    auto search=[&](const std::vector<V>& direction) {
      double damping=1;std::vector<V> last;
      for(std::size_t trial=0;trial<options.max_backtracks;++trial,damping*=.5) {
        auto candidate=current;
        for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<METAL_CN_SIZE;++k)candidate[i][k]+=damping*direction[i][k];
        if(candidate==current)break;if(candidate==last)continue;last=candidate;
        try {
          const auto test=evaluate(candidate,false);
          const auto remaining=newton_correction(candidate,test,state);
          double trial_merit=measure(test.balance);for(const auto& row:remaining)trial_merit=std::max(trial_merit,measure(row));
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
        for(std::size_t k=0;k<METAL_CN_SIZE;++k) {
          if(limited[i][k]<-current[i][k]){limited[i][k]=-.9*current[i][k];changed=true;}
          const double step=limited[i][k]*(k==METAL_CN_D?dunit:1.);
          if(step>0)positive+=step;else negative+=step;
        }
        const auto c=composition(current[i]);const double helium=cn_physical_ledger(c,*c.cn_molality).helium4;
        if(positive+negative>helium && positive>0) {
          const double fraction=(.9*helium-negative)/positive;
          for(double& x:limited[i])if(x>0)x*=fraction;changed=true;
        }
      }
      if(changed)search(limited);
    }
    if(!accepted) {
      std::ostringstream error;error<<std::setprecision(4)<<"CN diffusion: line search failed at iteration "<<iteration
        <<", correction "<<change<<", balance "<<measure(state.balance)<<", last rejection "<<last_rejection;
      throw std::runtime_error(error.str());
    }
    current=std::move(best);
  }
  throw std::runtime_error("CN diffusion: iteration limit");
}
std::vector<Composition> burn_metal_cn_and_transport(const Model& thermal,const Model& previous,
    const PPCNNetwork& nuclear,const MixingRegions& regions,std::span<const double> diffusivity,
    double dt,double tolerance) {
  if(!diffusivity.empty() && diffusivity.size()+1!=thermal.size())
    throw std::invalid_argument("metal CN mixing: diffusivity dimensions differ");
  std::vector<double> rates(diffusivity.size());
  for(std::size_t i=0;i<rates.size();++i) {
    if(!std::isfinite(diffusivity[i]) || diffusivity[i]<0)
      throw std::invalid_argument("metal CN mixing: invalid diffusivity");
    const double r=.5*(thermal.r(i)+thermal.r(i+1)),rho=.5*(thermal.rho(i)+thermal.rho(i+1));
    const double area=4*std::acos(-1.)*r*r*rho;
    rates[i]=area*area*diffusivity[i]/(thermal.m[i+1]-thermal.m[i]);
  }
  const MetalCNFlux no_drift=[](std::size_t,const Composition&,const Composition&,bool) {
    return MetalCNFaceResponse{};
  };
  SpeciesTransportOptions options;options.abundance_tolerance=tolerance;
  return burn_metal_cn_and_diffuse(thermal,previous,nuclear,regions,no_drift,dt,options,rates).composition;
}

MetalFluxReconstruction reconstruct_metal_fluxes(const Model& current,
    const Model& previous,const PPCNNetwork& nuclear,const MixingRegions& regions,
    std::span<const MetalCNBoundaryFlux> boundaries,double dt) {
  if(current.m!=previous.m || current.size()!=previous.size()
      || current.comp.size()!=current.size() || previous.comp.size()!=previous.size()
      || current.M!=previous.M || !(current.M>0) || !std::isfinite(current.M)
      || !(dt>0) || !std::isfinite(dt) || regions.empty()
      || boundaries.size()!=regions.size()-1)
    throw std::invalid_argument("species reconstruction: invalid models, regions, boundaries or time");
  for(const auto* model:{&current,&previous})for(const auto& c:model->comp)
    if(c[Species::H2]!=0)
      throw std::invalid_argument("metal heat reconstruction: initial D needs an isotope-aware heat provider");
  using R=MetalSpeciesVector;
  auto collapsed=[](const MetalCNVector& x)->R{return {x[0],x[1],x[2]+x[3]+x[4]+x[5]};};
  auto abundance=[](const Composition& c,std::size_t k){return k<2?c.X[k]:c.Z();};
  const auto weights=nodal_mass_weights(current);
  for(const auto& c:current.comp)check(c,previous.comp.front());
  for(const auto& c:previous.comp)check(c,previous.comp.front());
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
      detail::flux_finite(boundaries[r].rate);
    }
  }
  if(next!=current.size())throw std::invalid_argument("species reconstruction: incomplete region partition");
  MetalFluxReconstruction result;
  result.face_rates.resize(current.size()-1);result.cell_balances.resize(current.size());
  result.region_balances.resize(regions.size());
  std::vector<std::array<long double,3>> change(current.size());
  std::array<long double,3> integrated{};
  const long double factor=static_cast<long double>(dt)/current.M;
  for(std::size_t r=0;r<regions.size();++r) {
    const auto [begin,end]=regions[r];
    const R left=r?collapsed(boundaries[r-1].rate):R{};
    const R right=r+1<regions.size()?collapsed(boundaries[r].rate):R{};
    std::array<long double,3> accumulated{};
    for(std::size_t i=begin;i<end;++i) {
      NuclearResponse reaction;reaction.state=nuclear.eval(current.T(i),current.rho(i),current.comp[i]);
      detail::flux_finite(reaction.state.dXdt);
      double metal_source=0;for(std::size_t j=3;j<METAL_END;++j)metal_source+=reaction.state.dXdt[j];
      for(std::size_t k=0;k<3;++k) {
        change[i][k]=(static_cast<long double>(weights[i])/current.M)*
          ((static_cast<long double>(abundance(current.comp[i],k))-abundance(previous.comp[i],k))
           -static_cast<long double>(dt)*(k<2?reaction.state.dXdt[k]:metal_source));
        accumulated[k]+=change[i][k];integrated[k]+=change[i][k];
        if(i+1<end)result.face_rates[i][k]=static_cast<double>(left[k]-accumulated[k]/factor);
      }
      if(i+1<end)detail::flux_finite(result.face_rates[i]);
    }
    if(end<current.size())result.face_rates[end-1]=right;
    for(std::size_t k=0;k<3;++k)
      result.region_balances[r][k]=static_cast<double>(accumulated[k]+factor*(right[k]-left[k]));
    detail::flux_finite(result.region_balances[r]);
  }
  for(std::size_t i=0;i<current.size();++i) {
    const R left=i?result.face_rates[i-1]:R{};
    const R right=i+1<current.size()?result.face_rates[i]:R{};
    for(std::size_t k=0;k<3;++k)
      result.cell_balances[i][k]=static_cast<double>(change[i][k]+factor*
        (static_cast<long double>(right[k])-left[k]));
    detail::flux_finite(result.cell_balances[i]);
  }
  for(std::size_t k=0;k<3;++k)result.integrated_balance[k]=static_cast<double>(integrated[k]);
  detail::flux_finite(result.integrated_balance);return result;
}

} // namespace ember

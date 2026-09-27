// Coupled 0.5-Msun regression with analytic EOS/opacity and a grey
// boundary. This tests the integration, not a physical atmosphere or track.
#include "ember/stellar_seed.hpp"
#include "ember/eos_composite.hpp"
#include "ember/metal_cn_transport.hpp"
#include "ember/deuterium_burning.hpp"
#include "ember/controller.hpp"
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace ember;
namespace {
class ControlOpacity final:public Opacity {
 public:
  OpacityState eval(double T,double rho,const Composition& c)const override {
    const double absorption=1e23*rho*std::pow(T,-3.5),base=absorption+.2*(1+c.X[0]);
    const double z=std::log(T/3e5)/.5,bump=20*std::exp(-z*z);
    return {base*(1+bump),-3.5*absorption/base-4*z*bump/(1+bump),absorption/base};
  }
  const char* name()const override{return "analytic opacity for coupled regression";}
};
void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
}
int main(int argc,char**) {
  const auto start=std::clock();
  try {
    CompositeEos eos;ControlOpacity opacity;GreyAtmosphere atmosphere(eos,opacity);
    PPCNNetwork nuclear;
    auto c=solar_scaled(.7,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.001;c.X[2]-=.001;
    c.cn_molality=initial_gs98_cn(c);c=explicit_cn_material(c);
    Physics physics{&eos,&opacity,&nuclear,1.9,ConvectiveCriterion::ledoux};
    physics.explicit_metal_mixing_only=true;
    auto seed_physics=physics;
    const auto grid=argc>1?LuminosityGrid::volume_faces:LuminosityGrid::mass_nodes;
    auto initial=example::stellar_seed(256,.5*constants::Msun,.6*constants::Rsun,c,nuclear,atmosphere,3.,0.,grid);
    RelaxationOptions ro;ro.max_iterations=200;
    const auto relaxed=relax(initial,seed_physics,atmosphere,ro);
    require(relaxed.converged,"seed: "+relaxed.message);initial=relaxed.model;
    const auto controller_initial=initial;
    // Declared numerical coupling test: inject trace D into the relaxed
    // thermal control. This is not a formation model or a production restart.
    for(auto& q:initial.comp){q.X[0]-=2e-8;q[Species::H2]=2e-8;}
    EvolutionOptions options;options.relaxation=ro;options.abundance_tolerance=1e-14;
    options.max_coupling_iterations=80;
    constexpr double dt=1000*31557600.;
    const auto full=evolve_step(initial,physics,atmosphere,dt,options);
    require(full.converged,"combined step: "+full.message);
    const auto half=evolve_step(initial,physics,atmosphere,dt/2,options);
    require(half.converged,"first half: "+half.message);
    const auto end=evolve_step(half.model,physics,atmosphere,dt/2,options);
    require(end.converged,"second half: "+end.message);
    double time_error=0,first_law=0,mass_error=0,initial_D=0,final_D=0;
    const auto weights=nodal_mass_weights(initial);
    for(std::size_t i=0;i<initial.size();++i) {
      const auto& a=full.model.y[i];const auto& q=end.model.y[i];
      time_error=std::max({time_error,std::abs(a.lnr-q.lnr),std::abs(a.lnrho-q.lnrho),std::abs(a.lnT-q.lnT)});
      require(full.model.comp[i].cn_molality.has_value(),"CN inventory lost");
      require(full.model.comp[i][Species::H2]<=initial.comp[i][Species::H2]*(1+1e-12),"D increased");
      initial_D+=weights[i]/initial.M*initial.comp[i][Species::H2];
      final_D+=weights[i]/initial.M*end.model.comp[i][Species::H2];
    }
    for(const auto* s:{&full,&half,&end}) {
      first_law=std::max(first_law,std::abs(s->luminosity_balance));
      mass_error=std::max(mass_error,std::abs(s->nuclear_mass_balance));
    }
    require(time_error<1e-4,"unresolved control interval");
    require(final_D<.99*initial_D,"total D did not burn");
    require(first_law<2e-7 && mass_error<2e-5,"coupled energy accounting");
    require(full.model.age==initial.age+dt && end.model.age==initial.age+dt,"clock mismatch");
    require(full.model.luminosity_grid==grid && end.model.luminosity_grid==grid,"luminosity grid changed");
    if(argc>2) {
      auto predicted=initial;
      for(std::size_t i=0;i<initial.size();++i)for(std::size_t k=0;k<NVAR;++k) {
        const auto v=static_cast<Var>(k);
        predicted.y[i][v]+=2*(half.model.y[i][v]-initial.y[i][v]);
      }
      auto predicted_options=options;predicted_options.initial_structure_guess=&predicted;
      const auto predicted_step=evolve_step(initial,physics,atmosphere,dt,predicted_options);
      require(predicted_step.converged,"predicted starting structure failed: "+predicted_step.message);
      double prediction_difference=0;
      for(std::size_t i=0;i<initial.size();++i) {
        for(const auto v:{Var::lnr,Var::lnrho,Var::lnT})
          prediction_difference=std::max(prediction_difference,std::abs(predicted_step.model.y[i][v]-full.model.y[i][v]));
        for(std::size_t k=0;k<initial.comp[i].X.size();++k)
          prediction_difference=std::max(prediction_difference,std::abs(predicted_step.model.comp[i].X[k]-full.model.comp[i].X[k]));
      }
      require(prediction_difference<1e-7,"starting prediction changed the converged solution");
      require(std::abs(predicted_step.luminosity_balance)<2e-7 && std::abs(predicted_step.nuclear_mass_balance)<2e-5,
              "prediction changed energy accounting");
      predicted.m[0]*=1.01;bool mismatch_rejected=false;
      try { (void)evolve_step(initial,physics,atmosphere,dt,predicted_options); }
      catch(const std::invalid_argument&) {mismatch_rejected=true;}
      require(mismatch_rejected,"starting prediction from another mesh was accepted");
      // Fault injection tests recovery policy independently of stellar accuracy.
      // Every accepted trial still passes the normal energy checks above.
      EvolutionControlOptions control;control.step=options;
      control.predict_structure=true;
      control.target_age=controller_initial.age+dt;control.maximum_dt=dt;
      control.audit_failure_is_fatal=false;control.maximum_consecutive_rejections=3;
      EvolutionControlHooks hooks;
      std::vector<double> physics_ages,option_ages;
      hooks.physics=[&](const Model& m)->const Physics& {physics_ages.push_back(m.age);return physics;};
      hooks.configure_step=[&](const Model& m,EvolutionOptions&) {option_ages.push_back(m.age);};
      hooks.species_difference=[](const Composition& a,const Composition& b) {
        double error=0;
        for(std::size_t k=0;k<a.X.size();++k)error=std::max(error,std::abs(a.X[k]-b.X[k]));
        return error;
      };
      hooks.assess=[](const Model&,std::span<const std::array<double,3>>) {};
      hooks.cpu_seconds=[] {return 0.;};
      std::size_t calls=0;bool fail_all=false;
      hooks.audit=[&](const Model&,const EvolutionStep& step,double) {
        const bool injected=fail_all || calls==0;++calls;
        EvolutionAudit audit; audit.pass=step.converged && !injected
          && std::abs(step.luminosity_balance)<2e-7 && std::abs(step.nuclear_mass_balance)<2e-5;
        return audit;
      };
      std::vector<EvolutionAttempt> attempts;
      hooks.attempted=[&](const EvolutionAttempt& attempt) {attempts.push_back(attempt);};
      EvolutionState state{controller_initial,dt,0,0};
      const auto recovered=evolve(state,atmosphere,control,hooks);
      require(recovered.requested_age_reached && state.rejected==1,
              "controller must recover from one failed audit: "+recovered.stop_reason
              +", rejected="+std::to_string(state.rejected)+", last error="
              +std::to_string(attempts.back().error_norm));
      require(attempts.size()>=2 && !attempts[0].accepted && !attempts[0].audit_pass,
              "failed full-step audit must reject the entire trial");
      require(attempts[1].start_age==controller_initial.age && attempts[1].dt==dt/2,
              "retry must start from the retained state at half the duration");
      require(physics_ages==option_ages && physics_ages.size()==2*attempts.size(),
              "each half interval must select physics and options from its own starting state");
      for(std::size_t i=0;i<attempts.size();++i)
        require(physics_ages[2*i]==attempts[i].start_age
            && physics_ages[2*i+1]==attempts[i].start_age+attempts[i].dt/2,
            "second half reused the first half's starting state");
      for(const auto& attempt:attempts)if(attempt.accepted)
        require(attempt.audit_pass && attempt.error_norm<=1,"failed trial was accepted");
      fail_all=true;attempts.clear();state={controller_initial,dt,0,0};
      const auto exhausted=evolve(state,atmosphere,control,hooks);
      require(!exhausted.requested_age_reached && state.accepted==0 && state.rejected==3,
              "persistent audit failure must stop at the rejection bound");
      require(state.model.age==controller_initial.age && state.model.comp==controller_initial.comp
          && state.metal_heat_rates.empty(),"failed trials changed the retained state");
      for(std::size_t i=0;i<controller_initial.size();++i)for(std::size_t k=0;k<NVAR;++k)
        require(state.model.y[i][static_cast<Var>(k)]==controller_initial.y[i][static_cast<Var>(k)],
                "failed trials changed retained structure");
      control.audit_failure_is_fatal=true;state={controller_initial,dt,0,0};
      const auto fatal=evolve(state,atmosphere,control,hooks);
      require(!fatal.requested_age_reached && state.rejected==1 && state.accepted==0,
              "explicit fatal-audit policy must still stop immediately");
      // An optional extrapolation must never turn a valid half-step result
      // into a failed evolution. Exercise actual controller fallbacks.
      fail_all=false;calls=1;control.audit_failure_is_fatal=false;
      state={controller_initial,dt,0,0};
      const auto ordinary=evolve(state,atmosphere,control,hooks);
      require(ordinary.requested_age_reached,"extrapolation control did not finish");
      const auto reference=state;
      control.richardson_extrapolation=true;
      bool missing_assessment=false;
      try {(void)evolve(state,atmosphere,control,hooks);}
      catch(const std::invalid_argument&) {missing_assessment=true;}
      require(missing_assessment,"extrapolation ran without a physical assessment");
      for(bool throws:{false,true}) {
        std::size_t assessed=0;
        hooks.assess_extrapolated=[&](const EvolutionState&,const EvolutionStep&,const EvolutionStep&,
            const EvolutionStep&,const Model&,const std::vector<std::array<double,3>>&,double)->std::string {
          ++assessed;if(throws)throw std::domain_error("injected unsupported candidate");
          return "injected decline";
        };
        calls=1;state={controller_initial,dt,0,0};
        const auto fallback=evolve(state,atmosphere,control,hooks);
        require(fallback.requested_age_reached && assessed>0 && fallback.richardson_accepted==0
          && fallback.richardson_declined>0,"extrapolation assessment did not fall back");
        require(state.model.comp==reference.model.comp && state.metal_heat_rates==reference.metal_heat_rates,
                "declining extrapolation changed abundances or heat history");
        for(std::size_t i=0;i<state.model.size();++i)for(std::size_t k=0;k<NVAR;++k)
          require(state.model.y[i][static_cast<Var>(k)]==reference.model.y[i][static_cast<Var>(k)],
                  "declining extrapolation changed the retained half-step result");
      }
      hooks.assess_extrapolated=[](const EvolutionState&,const EvolutionStep&,const EvolutionStep&,
          const EvolutionStep&,const Model&,const std::vector<std::array<double,3>>&,double){return std::string{};};
      calls=1;state={controller_initial,dt,0,0};
      const auto extrapolated=evolve(state,atmosphere,control,hooks);
      require(extrapolated.requested_age_reached && extrapolated.richardson_accepted>0,
              "controller did not accept an assessed admissible extrapolation");

    }
    std::cout<<std::setprecision(17)<<"{\"outcome\":\"passed\",\"scope\":\"analytic EOS/grey coupled 0.5 solar mass test with trace-D injection; not a physical track\""
      <<",\"face_luminosities\":"<<face_luminosities(initial)
      <<",\"full_half_structure_error\":"<<time_error<<",\"maximum_first_law_error\":"<<first_law
      <<",\"maximum_mass_defect_error\":"<<mass_error<<",\"mean_initial_D\":"<<initial_D
      <<",\"mean_final_D\":"<<final_D<<",\"cpu_seconds\":"<<double(std::clock()-start)/CLOCKS_PER_SEC<<"}\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

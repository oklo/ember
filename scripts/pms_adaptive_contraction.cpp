// Contracting PMS seed and first coupled interval with specified initial D/He3.
// No post-seed fuel injection. Requires a separately accepted composition-
// dependent atmosphere covering the actual model; there is no grey fallback.
// CLI: six material inputs, atmosphere, entropy-loss rate, years,
//      nine-species composition file, radius/Rsun, trial Teff, threads, work,
//      optional --seed-only, --resume CHECKPOINT, --stop-after N,
//      --max-log-g G, --max-teff K (planning stops inside verified tables).
// years is the target age from the original PMS seed, including on restart.
// Resume requires identical executable, material and atmosphere inputs.
// --refine-abundances-from CHECKPOINT --source-executable FILE explicitly
// reads a predecessor with tolerance 1e-14; all material/atmosphere inputs
// must still match. The new solve and checkpoints use 1e-15. Compatibility
// of the predecessor's equations requires an independent saved-step replay.
// Fully mixed pp + initial-D contraction. Partial mixing and metal transport
// remain explicitly outside the audit domain. --total-energy-error scales
// time discretization to max(nuclear power, photon power), retaining separate
// isotope accuracy and unchanged per-interval conservation checks.
// --upgrade-audit-from reads a pinned 1e-15 predecessor without altering its state.
#include "ember/stellar_seed.hpp"
#include "ember/eos_variable_metal.hpp"
#include "ember/boundary.hpp"
#include "ember/opacity_mixture.hpp"
#include "ember/atmosphere_grid.hpp"
#include "ember/conduction_table.hpp"
#include "ember/evolution.hpp"
#include "ember/deuterium.hpp"
#include "ember/eos_deuterium.hpp"
#include "ember/atmosphere_deuterium.hpp"
#include "ember/evolution_checkpoint.hpp"
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace ember;
#include <limits>

#include "pms_interval_audit.hpp"
#include "pms_atmosphere_extension.hpp"

struct ContractingSeed final : Nuclear {
  const Nuclear& nuclear;
  const double entropy_loss;
  ContractingSeed(const Nuclear& n, double rate) : nuclear(n), entropy_loss(rate) {}
  NuclearState eval(double T, double rho, const Composition& c) const override {
    auto s = nuclear.eval(T, rho, c);
    const double cooling = T * entropy_loss, total = s.eps + cooling;
    s.dlneps_dlnT = (s.eps * s.dlneps_dlnT + cooling) / total;
    s.dlneps_dlnRho *= s.eps / total;
    s.eps = total;
    return s;
  }
  const char* name() const override { return "initial-model entropy loss only"; }
};

struct ReportedAtmosphere final : Atmosphere {
  const Atmosphere& source;
  explicit ReportedAtmosphere(const Atmosphere& a):source(a){}
  AtmosphereState eval(double T,double g,const Composition& c) const override {
    try {return source.eval(T,g,c);}
    catch(const std::domain_error& e) {
      std::ostringstream message;message<<std::setprecision(17)<<e.what()
        << "; surface query Teff=" << T << " logg=" << std::log10(g)
        << " XH=" << c[Species::H1] << " XHe3=" << c[Species::He3]
        << " XD=" << c[Species::H2];
      throw std::domain_error(message.str());
    }
  }
  const char* name() const override {return source.name();}
};

int main(int argc, char** argv) {
  if(argc<15)return 2;
  bool seed_only=false;
  std::string resume,source_executable,source_atmosphere;
  bool refine_abundances=false, upgrade_audit=false, extend_atmosphere=false, total_energy_error=false;
  double maximum_step_years=5000, structure_tolerance=1e-4;
  constexpr double abundance_tolerance=1e-15;
  std::size_t stop_after=500;
  double maximum_log_g=3.235,maximum_teff=std::numeric_limits<double>::infinity();
  for(int i=15;i<argc;++i) {
    const std::string flag=argv[i];
    if(flag=="--seed-only")seed_only=true;
    else if((flag=="--resume" || flag=="--refine-abundances-from" || flag=="--upgrade-audit-from" || flag=="--extend-atmosphere-from") && i+1<argc && resume.empty()) {
      resume=argv[++i];refine_abundances=flag=="--refine-abundances-from";upgrade_audit=flag=="--upgrade-audit-from";
      extend_atmosphere=flag=="--extend-atmosphere-from";
    }else if(flag=="--source-executable" && i+1<argc && source_executable.empty())source_executable=argv[++i];
    else if(flag=="--source-atmosphere" && i+1<argc && source_atmosphere.empty())source_atmosphere=argv[++i];
    else if(flag=="--stop-after" && i+1<argc) {
      const std::string n=argv[++i];
      if(n.empty() || n.find_first_not_of("0123456789")!=std::string::npos)return 2;
      stop_after=std::stoul(n);if(stop_after==0 || stop_after>500)return 2;
    }else if(flag=="--total-energy-error")total_energy_error=true;
    else if((flag=="--max-log-g" || flag=="--max-teff" || flag=="--max-step-years" || flag=="--structure-error") && i+1<argc) {
      const std::string number=argv[++i];std::size_t used=0;double value;
      try {value=std::stod(number,&used);}catch(const std::exception&){return 2;}
      if(used!=number.size() || !std::isfinite(value) || value<=0)return 2;
      if(flag=="--max-log-g")maximum_log_g=value;
      else if(flag=="--max-teff")maximum_teff=value;
      else if(flag=="--max-step-years"){if(value<10 || value>1e7)return 2;maximum_step_years=value;}
      else {if(value>1e-4)return 2;structure_tolerance=value;}
    }else return 2;
  }
  if(seed_only && !resume.empty())return 2;
  if((refine_abundances || upgrade_audit || extend_atmosphere)!=!source_executable.empty())return 2;
  if(extend_atmosphere!=!source_atmosphere.empty())return 2;
  std::ostringstream report; report << std::setprecision(17);
  try {
    const double entropy_loss = std::stod(argv[8]), years = std::stod(argv[9]);
    if (!std::isfinite(entropy_loss) || !std::isfinite(years) ||
        entropy_loss <= 0 || years <= 0) throw std::invalid_argument("positive finite controls required");
    VariableMetalHelmholtzEos table_eos(argv[1], HelmholtzTableEos::Mixture::allow_documented_proxy);
    DeuteriumApproxEos eos(table_eos);
    MixtureOpacity molecular(argv[2]), warm(argv[3]), bridge(argv[4]), hot(argv[5]);
    BlendedOpacity mid(warm, bridge, 5.05, 5.10), upper(mid, hot, 5.6, 5.7);
    BlendedOpacity raw(molecular, upper, 4.4, 4.47);
    auto radiation = std::make_shared<ElementalOpacity>(raw);
    TabulatedConduction table(argv[6]);
    CombinedOpacity opacity(radiation, std::make_shared<HotConduction>(table));
    CompositionAtmosphereGrid table_atmosphere(eos, argv[7], CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
    TraceDeuteriumAtmosphere isotope_atmosphere(eos,table_atmosphere);
    ReportedAtmosphere atmosphere(isotope_atmosphere);
    PPDeuterium nuclear;
    ContractingSeed seed_source(nuclear, entropy_loss);
    Physics physics{&eos, &opacity, &nuclear, 1.9, ConvectiveCriterion::ledoux};
    auto seed_physics = physics; seed_physics.nuclear = &seed_source;
    Composition c{}; c.basis = AbundanceBasis::baryon_mass;
    c.metal_inventory = MetalInventory::gs98;
    std::ifstream composition_input(argv[10]);
    for (auto& x:c.X) {
      composition_input >> x;
      if(!composition_input || !std::isfinite(x) || x<0)throw std::invalid_argument("invalid initial isotope input");
    }
    std::string extra;
    if(composition_input>>extra || std::abs(c.sum()-1)>2e-12 || c[Species::H2]<=0)
      throw std::invalid_argument("initial isotope input is not normalized or lacks D");
    const double radius=std::stod(argv[11]),teff=std::stod(argv[12]);
    const auto threads=std::stoul(argv[13]);
    if(!std::isfinite(radius+teff) || radius<=0 || teff<=0 || threads<1 || threads>4)
      throw std::invalid_argument("invalid physical starting controls");
    const std::filesystem::path work(argv[14]);
    if(!std::filesystem::is_directory(work))throw std::invalid_argument("missing output directory");
    for(const auto* name:{"seed","full","half1","final"})
      if(std::filesystem::exists(work/(std::string(name)+".checkpoint")))throw std::invalid_argument("checkpoint output already exists");
    driver::Identities identities{{"executable",driver::file_identity(argv[0])}};
    for(int i=1;i<=7;++i)identities[std::filesystem::canonical(argv[i]).string()]=driver::file_identity(argv[i]);
    identities[std::filesystem::canonical(argv[10]).string()]=driver::file_identity(argv[10]);
    const driver::Selections selections{"PMS contraction with initial D","PPDeuterium","GS98 variable EOS", "composition atmosphere", "Ledoux instantaneous convection"};
    double checkpoint_step_years=1000;
    std::size_t checkpoint_rejected=0;
    const auto checkpoint=[&](const char* name,const Model& m,std::size_t accepted) {
      driver::Checkpoint state{m,checkpoint_step_years*31557600.,accepted,checkpoint_rejected};
      driver::write_checkpoint(work/(std::string(name)+".checkpoint"),state,selections,abundance_tolerance,identities);
    };
    Model current;
    std::size_t prior_accepted=0,prior_rejected=0;
    const auto start = std::clock();
    if(resume.empty()) {
    auto guess = example::stellar_seed(512, .1*constants::Msun, radius*constants::Rsun,
        c, seed_source, atmosphere, 1.5, teff);

    // Initial guess only: integrate the selected discrete stellar equations
    // inward through the envelope, then blend into the polytropic interior.
    // Global relaxation and the physical acceptance criteria remain unchanged.
    const Model polytrope = guess;
    const auto surface = atmosphere.eval(teff, constants::G*guess.M/
        std::pow(guess.r(guess.size()-1),2), c);
    guess.y.back().lnT=std::log(surface.T);
    guess.y.back().lnrho=std::log(surface.rho);
    const double initial_L = guess.y.back().L;
    std::ofstream envelope_log(work/"envelope_initializer.csv");
    envelope_log<<std::setprecision(17)<<"i,q,iterations,residual\n";
    auto norm4=[](const std::array<double,4>& a) {
      double n=0.;for(auto v:a)n=std::max(n,std::abs(v));return n;
    };
    for(std::size_t hi=guess.size()-1;hi>0 && guess.m[hi-1]/guess.M>.90;--hi) {
      const std::size_t lo=hi-1;const double dm=guess.m[hi]-guess.m[lo];
      const auto e=eos.eval(guess.T(hi),guess.rho(hi),c);
      const double Plo=e.P+constants::G*.5*(guess.m[hi]+guess.m[lo])*dm/
        (4*M_PI*std::pow(guess.r(hi),4));
      const double Tlo=guess.T(hi)*std::pow(Plo/e.P,e.grad_ad);
      guess.y[lo]={guess.y[hi].lnr-dm/(4*M_PI*std::pow(guess.r(hi),3)*guess.rho(hi)),
        std::log(eos.rho_from_PT(Tlo,Plo,c,guess.rho(hi))),std::log(Tlo),
        guess.y[hi].L-dm*seed_source.eval(guess.T(hi),guess.rho(hi),c).eps};
      bool converged=false;
      for(std::size_t it=0;it<40;++it) {
        const auto z=zone_residual(guess,lo,seed_physics,0.);
        std::array<double,4> f{};std::array<std::array<double,5>,4> mat{};
        for(std::size_t k=0;k<4;++k) {
          const double scale=k==2?dm/initial_L:dm;f[k]=z.f[k]*scale;
          for(std::size_t v=0;v<4;++v)
            mat[k][v]=z.dfdy_lo[k][v]*scale*(v==3?initial_L:1.);
          mat[k][4]=-f[k];
        }
        const double before=norm4(f);
        if(before<1e-10) {
          envelope_log<<lo<<','<<guess.m[lo]/guess.M<<','<<it<<','<<before<<'\n';
          converged=true;break;
        }
        for(std::size_t k=0;k<4;++k) {
          std::size_t pivot=k;
          for(std::size_t j=k+1;j<4;++j)if(std::abs(mat[j][k])>std::abs(mat[pivot][k]))pivot=j;
          if(std::abs(mat[pivot][k])<1e-15)throw std::runtime_error("singular envelope initializer");
          std::swap(mat[pivot],mat[k]);
          const double scale=mat[k][k];for(std::size_t j=k;j<5;++j)mat[k][j]/=scale;
          for(std::size_t i=0;i<4;++i)if(i!=k) {
            const double factor=mat[i][k];for(std::size_t j=k;j<5;++j)mat[i][j]-=factor*mat[k][j];
          }
        }
        double damping=1.;for(std::size_t v=0;v<4;++v)
          if(std::abs(mat[v][4])>.2)damping=std::min(damping,.2/std::abs(mat[v][4]));
        const auto base=guess.y[lo];bool accepted=false;
        for(std::size_t trial=0;trial<24;++trial) {
          for(std::size_t v=0;v<4;++v)guess.y[lo][static_cast<Var>(v)]=
            base[static_cast<Var>(v)]+damping*mat[v][4]*(v==3?initial_L:1.);
          try {
            auto candidate=zone_equations(guess,lo,seed_physics,0.);
            for(std::size_t k=0;k<4;++k)candidate[k]*=k==2?dm/initial_L:dm;
            if(norm4(candidate)<before){accepted=true;break;}
          }catch(const std::exception&){}
          damping*=.5;
        }
        if(!accepted)throw std::runtime_error("envelope initializer line search failed at "+std::to_string(lo));
      }
      if(!converged)throw std::runtime_error("envelope initializer iteration limit at "+std::to_string(lo));
    }
    for(std::size_t i=0;i<guess.size();++i) {
      double f=std::clamp((guess.m[i]/guess.M-.90)/.08,0.,1.);f=f*f*(3-2*f);
      for(std::size_t v=0;v<4;++v)guess.y[i][static_cast<Var>(v)]=
        f*guess.y[i][static_cast<Var>(v)]+(1-f)*polytrope.y[i][static_cast<Var>(v)];
    }

    const auto diagnose=[&](const char* name,const Model& m) {
      std::ofstream out(work/(std::string(name)+".csv"));out<<std::setprecision(17);
      out<<"q,r_cm,rho,T,P,S,grad_ad,cp,L,kappa,eps,mass_residual,pressure_residual,energy_residual,transport_residual\n";
      for(std::size_t i=0;i<m.size();++i) {
        const auto e=eos.eval(m.T(i),m.rho(i),m.comp[i]);
        const auto k=opacity.eval(m.T(i),m.rho(i),m.comp[i]);
        const auto n=seed_source.eval(m.T(i),m.rho(i),m.comp[i]);
        std::array<double,4> z{};
        if(i+1<m.size()) {
          z=zone_equations(m,i,seed_physics,0.);
          const double dm=m.m[i+1]-m.m[i];
          const double unit=std::max({std::abs(guess.y[i].L),std::abs(guess.y[i+1].L),
                guess.y.back().L*m.m[i+1]/m.M});
          for(std::size_t j=0;j<4;++j)z[j]*=(j==2?dm/unit:dm);
        }
        out<<m.m[i]/m.M<<','<<m.r(i)<<','<<m.rho(i)<<','<<m.T(i)<<','<<e.P<<','<<e.S<<','
          <<e.grad_ad<<','<<e.cp<<','<<m.y[i].L<<','<<k.kappa<<','<<n.eps;
        for(auto v:z)out<<','<<v;out<<'\n';
      }
      const auto a=surface_residual(m.y.back(),m.M,m.comp.back(),eos,atmosphere);
      const auto c=central_residual(m,seed_physics);
      std::ofstream bound(work/(std::string(name)+"_boundary.json"));bound<<std::setprecision(17)
        <<"{\"surface\":["<<a.f[0]<<','<<a.f[1]<<"],\"central\":["<<c.f[0]<<','<<c.f[1]
        <<"],\"Teff\":"<<a.Teff<<",\"log_g\":"<<std::log10(a.gravity)<<"}\n";
    };
    diagnose("initial",guess);
    RelaxationOptions ro; ro.zone_threads = threads;
    auto seed = relax(guess, seed_physics, atmosphere, ro);
    diagnose("last",seed.model);
    report << "{\"seed_converged\":" << seed.converged << ",\"seed_message\":" << std::quoted(seed.message);
    report << ",\"seed_only\":" << seed_only << ",\"seed_iterations\":" << seed.iterations
      << ",\"seed_residual\":" << seed.residual << ",\"seed_history\":[";
    bool first_iteration=true;
    for(const auto& h:seed.history) {
      if(!first_iteration)report<<',';first_iteration=false;
      report << "{\"residual\":" << h.residual << ",\"correction\":" << h.correction
        << ",\"damping\":" << h.damping << ",\"linear_error\":" << h.linear_error << '}';
    }
    report << ']';
    if (!seed.converged) {
      const auto& m=seed.model;
      const double r=m.r(m.size()-1);
      report << ",\"last_radius_Rsun\":" << r/constants::Rsun
        << ",\"last_Teff_K\":" << std::pow(m.y.back().L/(4*M_PI*constants::sigma_SB*r*r),.25)
        << ",\"last_central_T_K\":" << m.T(0) << "}\n";
      std::cout << report.str(); return 1;
    }
    // Relaxation holds the actual initial isotope abundances fixed. No
    // deuterium or helium is inserted after constructing the seed.
    for(const auto& comp:seed.model.comp)if(comp!=c)throw std::logic_error("seed changed initial fuel");
    checkpoint("seed",seed.model,0);
    if(seed_only) {
      const double r=seed.model.r(seed.model.size()-1),lum=seed.model.y.back().L;
      report << ",\"radius_Rsun\":" << r/constants::Rsun
        << ",\"luminosity_Lsun\":" << lum/constants::Lsun
        << ",\"Teff_K\":" << std::pow(lum/(4*M_PI*constants::sigma_SB*r*r),.25)
        << ",\"central_T_K\":" << seed.model.T(0)
        << ",\"entropy_loss_seed_only\":" << entropy_loss
        << ",\"accepted_evolutionary_intervals\":0,\"cpu_seconds\":"
        << double(std::clock()-start)/CLOCKS_PER_SEC << "}\n";
      std::cout << report.str();return 0;
    }
    current=seed.model;
    }else {
      auto original_identities=identities;
      if(refine_abundances || upgrade_audit || extend_atmosphere)original_identities["executable"]=driver::file_identity(source_executable);
      std::size_t appended_atmosphere_states=0;
      if(extend_atmosphere) {
        appended_atmosphere_states=pms::check_atmosphere_extension_files(source_atmosphere,argv[7]);
        original_identities.erase(std::filesystem::canonical(argv[7]).string());
        original_identities[std::filesystem::canonical(source_atmosphere).string()]=driver::file_identity(source_atmosphere);
      }
      const auto restored=driver::read_checkpoint(resume,512,.1*constants::Msun,c,
          selections,refine_abundances?1e-14:abundance_tolerance,original_identities);
      current=restored.model;
      if(extend_atmosphere) {
        CompositionAtmosphereGrid old_table(eos,source_atmosphere,CompositionAtmosphereGrid::Mixture::allow_documented_proxy);
        TraceDeuteriumAtmosphere old_atmosphere(eos,old_table);
        const double r=current.r(current.size()-1),g=constants::G*current.M/(r*r);
        const double t=std::pow(current.y.back().L/(4*M_PI*constants::sigma_SB*r*r),.25);
        const auto fields=[](const AtmosphereState& s){return std::array{s.T,s.P,s.Pgas,s.rho,s.tau,s.dlnT_dlnTeff,s.dlnT_dlng,s.dlnP_dlnTeff,s.dlnP_dlng};};
        if(fields(old_atmosphere.eval(t,g,current.comp.back()))!=fields(atmosphere.eval(t,g,current.comp.back())))
          throw std::runtime_error("extended atmosphere changes the saved model's values or derivatives");
      }
      prior_accepted=restored.accepted;prior_rejected=restored.rejected;
      checkpoint_step_years=restored.next_dt/31557600.;
      checkpoint_rejected=prior_rejected;
      if(current.age>=years*31557600. || checkpoint_step_years<10 || checkpoint_step_years>1e7)
        throw std::invalid_argument("restart age or next timestep is outside the requested initial contraction");
      checkpoint("seed",current,prior_accepted);
      report<<"{\"resumed\":true,\"parent_checkpoint\":"<<std::quoted(resume)
        <<",\"parent_age_years\":"<<current.age/31557600.
        <<",\"abundance_tolerance\":"<<abundance_tolerance
        <<",\"explicit_abundance_refinement\":"<<refine_abundances;
      report<<",\"explicit_audit_upgrade\":"<<upgrade_audit;
      report<<",\"explicit_atmosphere_extension\":"<<extend_atmosphere;
      if(extend_atmosphere)report<<",\"source_atmosphere\":"<<std::quoted(source_atmosphere)<<",\"appended_atmosphere_states\":"<<appended_atmosphere_states;
      if(refine_abundances || upgrade_audit || extend_atmosphere)report<<",\"source_executable\":"<<std::quoted(source_executable);
    }
    RelaxationOptions evolution_ro; evolution_ro.zone_threads=threads;
    EvolutionOptions options; options.relaxation=evolution_ro;options.abundance_tolerance=abundance_tolerance;

    double elapsed=current.age/31557600.,step_years=std::min(checkpoint_step_years,maximum_step_years);
    const double initial_age_years=elapsed;
    const double target_age=years*31557600.;
    // Full and half steps add seconds in different floating-point orders.
    // Preserve the model's actual clock, but do not attempt another solve
    // for an endpoint difference bounded by their representation precision.
    const double age_roundoff=8*(std::nextafter(target_age,
        std::numeric_limits<double>::infinity())-target_age);
    const auto target_reached=[&](){return current.age>=target_age
        || target_age-current.age<=age_roundoff;};
    unsigned accepted=0,rejected=0,attempt=0,failed=0;
    std::string stop="requested duration";
    std::ofstream history(work/"history.jsonl"),attempts(work/"attempts.jsonl");
    history<<std::setprecision(17);attempts<<std::setprecision(17);
    auto history_row=[&]() {
      const auto w=nodal_mass_weights(current);double D=0,S=0,lnuc=0;
      for(std::size_t i=0;i<current.size();++i) {
        D+=w[i]*current.comp[i][Species::H2];S+=w[i]*eos.eval(current.T(i),current.rho(i),current.comp[i]).S;
        lnuc+=w[i]*nuclear.eval(current.T(i),current.rho(i),current.comp[i]).eps;
      }
      const double r=current.r(current.size()-1),L=current.y.back().L;
      history<<"{\"years\":"<<elapsed<<",\"accepted\":"<<(prior_accepted+accepted)<<",\"radius_Rsun\":"<<r/constants::Rsun
        <<",\"Teff_K\":"<<std::pow(L/(4*M_PI*constants::sigma_SB*r*r),.25)<<",\"luminosity_Lsun\":"<<L/constants::Lsun
        <<",\"central_T_K\":"<<current.T(0)<<",\"log_g\":"<<std::log10(constants::G*current.M/(r*r))
        <<",\"D_mass_g\":"<<D<<",\"mean_entropy\":"<<S/current.M<<",\"nuclear_fraction\":"<<lnuc/L<<"}\n";history.flush();
    };
    history_row();
    while(!target_reached() && attempt<500 && accepted<stop_after) {
      if(double(std::clock()-start)/CLOCKS_PER_SEC>240){stop="CPU budget";break;}
      const double surface_radius=current.r(current.size()-1);
      if(std::log10(constants::G*current.M/(surface_radius*surface_radius))>maximum_log_g){stop="planned atmosphere gravity boundary";break;}
      const double surface_teff=std::pow(current.y.back().L/(4*M_PI*constants::sigma_SB*surface_radius*surface_radius),.25);
      if(surface_teff>maximum_teff){stop="planned atmosphere temperature boundary";break;}
      const double dy=std::min(step_years,years-elapsed),ds=dy*31557600.;
      const auto full=evolve_step(current,physics,atmosphere,ds,options);
      const auto h1=full.converged?evolve_step(current,physics,atmosphere,ds/2,options):full;
      const auto h2=h1.converged?evolve_step(h1.model,physics,atmosphere,ds/2,options):h1;
      const auto af=audit_pms_interval(current,full,ds,nuclear),a1=audit_pms_interval(current,h1,ds/2,nuclear),
        a2=audit_pms_interval(h1.model,h2,ds/2,nuclear);
      double structure=0,composition_error=0,energy_error=0,error=1e30;
      const bool converged=full.converged && h1.converged && h2.converged;
      if(converged) {
        for(std::size_t i=0;i<current.size();++i) {
          structure=std::max({structure,std::abs(full.model.y[i].lnr-h2.model.y[i].lnr),
            std::abs(full.model.y[i].lnrho-h2.model.y[i].lnrho),std::abs(full.model.y[i].lnT-h2.model.y[i].lnT)});
          for(std::size_t k=0;k<NSPEC;++k)composition_error=std::max(composition_error,std::abs(full.model.comp[i].X[k]-h2.model.comp[i].X[k]));
        }
        const double average_nuclear=.5*(h1.nuclear_luminosity+h2.nuclear_luminosity);
        const double energy_scale=total_energy_error?std::max(std::abs(average_nuclear),
            .5*(std::abs(h1.model.y.back().L)+std::abs(h2.model.y.back().L))):std::abs(average_nuclear);
        energy_error=std::abs(full.nuclear_luminosity-average_nuclear)/std::max(energy_scale,std::numeric_limits<double>::min());
        error=std::max({structure/structure_tolerance,composition_error/1e-8,energy_error/.005});
      }
      const bool audit=af.pass && a1.pass && a2.pass;
      const bool take=converged && audit && error<=1;
      attempts<<"{\"attempt\":"<<attempt<<",\"start_years\":"<<elapsed<<",\"dt_years\":"<<dy
        <<",\"converged\":"<<converged<<",\"audit_pass\":"<<audit<<",\"accepted\":"<<take
        <<",\"error_norm\":"<<error<<",\"structure_error\":"<<structure<<",\"abundance_error\":"<<composition_error
        <<",\"nuclear_energy_error\":"<<energy_error<<",\"full_message\":"<<std::quoted(full.message)
        <<",\"half1_message\":"<<std::quoted(h1.message)<<",\"half2_message\":"<<std::quoted(h2.message)
        <<",\"fuel_errors\":["<<af.fuel_error<<','<<a1.fuel_error<<','<<a2.fuel_error
        <<"],\"mass_errors\":["<<af.mass_error<<','<<a1.mass_error<<','<<a2.mass_error
        <<"],\"representation_allowances\":["<<af.roundoff<<','<<a1.roundoff<<','<<a2.roundoff
        <<"],\"first_law_errors\":["<<af.first_law<<','<<a1.first_law<<','<<a2.first_law<<"]}\n";attempts.flush();
      if(take) {
        const auto prefix="interval-"+std::to_string(accepted);std::filesystem::create_directory(work/prefix);
        checkpoint((prefix+"/seed").c_str(),current,prior_accepted+accepted);
        checkpoint((prefix+"/full").c_str(),full.model,prior_accepted+accepted+1);
        checkpoint((prefix+"/half1").c_str(),h1.model,prior_accepted+accepted+1);
        checkpoint((prefix+"/final").c_str(),h2.model,prior_accepted+accepted+1);
        current=h2.model;elapsed=current.age/31557600.;++accepted;failed=0;history_row();
        step_years=std::min(maximum_step_years,dy*std::clamp(.8/std::sqrt(std::max(error,1e-8)),.5,1.5));
      } else {
        ++rejected;
        // A failed physical audit needs diagnosis, not progressively smaller
        // steps that make weak-burning abundance subtraction less accurate.
        if(converged && !audit){stop="physical interval audit";break;}
        if(!converged && ++failed>=3){stop="repeated nonlinear failure";break;}
        step_years=dy*.5;if(step_years<10){stop="minimum interval";break;}
      }
      checkpoint_step_years=step_years;
      checkpoint_rejected=prior_rejected+rejected;
      ++attempt;
    }
    if(attempt>=500)stop="attempt budget";
    if(accepted>=stop_after && !target_reached())stop="accepted interval budget";
    checkpoint("final",current,prior_accepted+accepted);
    report<<",\"years\":"<<elapsed<<",\"initial_age_years\":"<<initial_age_years
      <<",\"advanced_years\":"<<elapsed-initial_age_years<<",\"requested_years\":"<<years<<",\"accepted_intervals\":"<<accepted
      <<",\"rejected_intervals\":"<<rejected<<",\"stop_reason\":"<<std::quoted(stop)
      <<",\"selected_for_physical_PMS_track\":false,\"cpu_seconds\":"<<double(std::clock()-start)/CLOCKS_PER_SEC<<"}\n";
    std::cout<<report.str();return (target_reached() || stop=="accepted interval budget")?0:1;
  } catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::cout << "{\"failed\":true,\"message\":" << std::quoted(e.what()) << "}\n"; return 1;
  }
}

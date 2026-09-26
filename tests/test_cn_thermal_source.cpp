#include "ember/cn_burning.hpp"
#include "ember/constants.hpp"
#include "ember/evolution_checkpoint.hpp"
#include <algorithm>
#include <cstdio>

using namespace ember;
int main() {
  int checks=0,failures=0;
  auto check=[&](bool ok,const char* s){++checks;if(!ok){++failures;std::printf("FAIL %s\n",s);}};
  auto rejects=[](auto f){try{f();}catch(const std::exception&){return true;}return false;};
  PPCNNetwork nuclear;
  double maximum_derivative_error=0;
  for(double X:{.2,.0006})for(double f:{0.,.5,1.}) {
    auto c=solar_scaled(X,.02);c.basis=AbundanceBasis::baryon_mass;c.metal_inventory=MetalInventory::gs98;
    c.X[1]=.001;c.X[2]-=.001;auto initial=initial_gs98_cn(c);
    c.cn_molality=CNAbundances{initial[0]*(1-f),0,initial[2]+f*initial[0]};
    const double T=1.4e7,rho=3000;
    const auto n=nuclear.composition_response(T,rho,c);
    const auto p=nuclear.pp().eval(T,rho,c),cn=nuclear.cn().response(T,rho,c,*c.cn_molality).physical.state;
    check(n.state.eps==p.eps+cn.eps && n.state.eps_neutrino==p.eps_neutrino+cn.eps_neutrino,"physical heat used once");
    double sum=0;for(double v:n.state.dXdt)sum+=v;
    check(std::abs(sum)<1e-12*std::abs(n.state.dXdt[0]),"lookup conserves baryons");
    for(std::size_t j=3;j<NSPEC;++j)check(n.state.dXdt[j]==0,"fixed table metal carriers");
    for(int variable=0;variable<2;++variable) {
      const double h=1e-5;
      const auto up=nuclear.eval(T*(variable==0?std::exp(h):1),rho*(variable==1?std::exp(h):1),c);
      const auto dn=nuclear.eval(T*(variable==0?std::exp(-h):1),rho*(variable==1?std::exp(-h):1),c);
      const double error=std::abs((std::log(up.eps)-std::log(dn.eps))/(2*h)
        -(variable==0?n.state.dlneps_dlnT:n.state.dlneps_dlnRho));
      maximum_derivative_error=std::max(maximum_derivative_error,error);
      check(error<2e-6,"thermal derivatives at fixed local catalysts");
    }
    Model m;m.M=10;m.m={1,4,10};
    for(int i=0;i<3;++i){m.y.push_back({std::log(i+1.),std::log(rho),std::log(T),1});m.comp.push_back(c);}
    const double dt=1e12;const auto before=m;
    m.comp=burn_and_mix(m,m,nuclear,{{0,3}},dt,1e-14);
    const auto w=nodal_mass_weights(m);double mass=0,heat=0;
    for(std::size_t i=0;i<m.size();++i) {
      check(m.comp[i].cn_molality.has_value(),"burn keeps CN inventory");
      const auto src=nuclear.eval(m.T(i),m.rho(i),m.comp[i]);heat+=w[i]*(src.eps+src.eps_neutrino);
      double delta=nuclear.rest_energy_correction(m.comp[i])-nuclear.rest_energy_correction(before.comp[i]);
      for(std::size_t j=0;j<NSPEC;++j)delta+=(m.comp[i].X[j]-before.comp[i].X[j])
        *(nuclides[j].A/mass_numbers[j]-1)*constants::c*constants::c;
      mass-=w[i]*delta/dt;
    }
    check(std::abs(mass/heat-1)<2e-6,"burn dispatch preserves physical rest energy");
    check(m.comp[0].cn_molality!=before.comp[0].cn_molality,"catalysts actually evolve");
    check(rejects([&]{burn_and_mix(m,m,nuclear.pp(),{{0,3}},dt);}),"incompatible pp-only evolution rejected");
    check(rejects([&]{burn_and_transport(m,m,nuclear.pp(),{{0,1},{1,3}},{1,1},dt);}),"CN transport with wrong network rejected");
    const auto path=std::filesystem::temp_directory_path()/"ember-cn-checkpoint-test.restart";
    if(std::filesystem::exists(path))throw std::runtime_error("test checkpoint already exists");
    driver::Checkpoint state{m,dt,1,0};driver::Selections selections{"cn","a","b","c","d"};
    driver::Identities identities{{"executable","test-only"}};
    driver::write_checkpoint(path,state,selections,1.,identities);
    const auto restored=driver::read_checkpoint(path,3,10,c,selections,1.,identities);
    check(restored.model.comp==m.comp && restored.model.age==m.age,"CN checkpoint exact round trip");
    auto inactive=c;inactive.cn_molality.reset();
    check(rejects([&]{driver::read_checkpoint(path,3,10,inactive,selections,1.,identities);}),"CN checkpoint requires active network selection");
    std::filesystem::remove(path);
    for(auto& cell:m.comp)cell.cn_molality.reset();
    state.model=m;driver::write_checkpoint(path,state,selections,1.,identities);
    check(driver::read_checkpoint(path,3,10,inactive,selections,1.,identities).model.comp==m.comp,"old checkpoint format retained");
    std::filesystem::remove(path);
  }
  std::printf("%d checks, %d failures; maximum thermal derivative error %.4g\n",checks,failures,maximum_derivative_error);
  return failures?1:0;
}

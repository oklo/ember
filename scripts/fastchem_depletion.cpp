// SPDX-License-Identifier: GPL-3.0-or-later
// Optional offline TLUSTY/SYNSPEC adapter; never linked into Ember itself.
// FastChem equilibrium removal with zero grain opacity (settled-grain limit).
#include "fastchem.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
constexpr const char* symbols[]={"H","He","Li","Be","B","C","N","O","F","Ne",
 "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn",
 "Fe","Co","Ni","Cu","Zn"};
struct Depletion {
  std::unique_ptr<fastchem::FastChem> chemistry;
  bool active{};
  std::vector<int> atomic_index;
  std::vector<std::vector<int>> gas_stoich,cond_stoich;
  std::vector<double> abundance;
  // One source process has one fixed bulk composition. MOLEQ already keeps
  // a separate cbase reservoir; retain our own copy defensively so a future
  // caller cannot accidentally recycle depleted gas as the bulk reservoir.
  std::vector<double> bulk_reference;
  unsigned h{},electron{};
  explicit Depletion(const double* bulk,int n) : bulk_reference(bulk,bulk+n) {
    std::ifstream input("ember-condensates.cfg");
    std::string mode,source,file;
    if(!std::getline(input,mode) || !std::getline(input,source) || !std::getline(input,file)
        || (mode!="gas" && mode!="equilibrium"))throw std::runtime_error("invalid depletion configuration");
    active=mode=="equilibrium";
    auto* old=std::cout.rdbuf(std::cerr.rdbuf());
    chemistry=std::make_unique<fastchem::FastChem>(file,source+"/input/logK/logK.dat",
        source+"/input/logK/logK_condensates.dat",0);
    std::cout.rdbuf(old);
    for(const auto* key:{"accuracyChem","accuracyElementConservation","accuracyCond"})
      if(!chemistry->setParameter(key,1e-8))throw std::runtime_error("unsupported chemistry tolerance");
    h=chemistry->getElementIndex("H");electron=chemistry->getElementIndex("e-");
    if(h>=chemistry->getElementNumber() || electron>=chemistry->getElementNumber())
      throw std::runtime_error("missing H or electron");
    abundance=chemistry->getElementAbundances();
    atomic_index.assign(abundance.size(),-1);
    for(unsigned k=0;k<abundance.size();++k) {
      const auto name=chemistry->getElementSymbol(k);
      for(int z=0;z<std::min(n,30);++z)if(name==symbols[z])atomic_index[k]=z;
      if(k==electron)continue;
      const int z=atomic_index[k];
      if(z<0 || std::abs(abundance[k]/abundance[h]/bulk[z]-1)>1e-10)
        throw std::runtime_error("depletion and opacity element numbers disagree");
    }
    for(int z=0;z<n;++z)if(bulk[z]>1e-90 &&
        std::find(atomic_index.begin(),atomic_index.end(),z)==atomic_index.end())
      throw std::runtime_error("element missing from depletion chemistry");
    for(unsigned j=0;j<chemistry->getGasSpeciesNumber();++j)gas_stoich.push_back(chemistry->getGasSpeciesStoichiometry(j));
    for(unsigned j=0;j<chemistry->getCondSpeciesNumber();++j)cond_stoich.push_back(chemistry->getCondSpeciesStoichiometry(j));
  }
  void evaluate(double T,double P,double* gas,int n) {
    if(static_cast<std::size_t>(n)!=bulk_reference.size())
      throw std::runtime_error("changed depletion element count");
    std::copy(bulk_reference.begin(),bulk_reference.end(),gas);
    if(!std::isfinite(T+P) || T<100 || P<1e-7 || P>1e9)
      throw std::domain_error("depletion outside T/P support");
    if(!active)return;
    // Explicit fully vaporized continuation above the chemistry domain.
    // The source preparation/audit must verify the 6000 K join over the
    // complete pressure interval used by the atmosphere/opacity calculation.
    if(T>6000)return;
    fastchem::FastChemInput input;fastchem::FastChemOutput out;
    input.temperature={T};input.pressure={P/1e6};input.equilibrium_condensation=true;
    if(chemistry->calcDensities(input,out)!=fastchem::FASTCHEM_SUCCESS || out.fastchem_flag[0]
       || std::any_of(out.element_conserved[0].begin(),out.element_conserved[0].end(),[](auto v){return !v;}))
      throw std::runtime_error("depletion source failed convergence/conservation");
    std::vector<double> atoms(abundance.size()),solid(abundance.size());double particles=0;
    for(unsigned j=0;j<gas_stoich.size();++j) {
      const double number=out.number_densities[0][j];
      if(!std::isfinite(number) || number<0)throw std::runtime_error("invalid gas density");
      particles+=number;
      for(unsigned k=0;k<atoms.size();++k)atoms[k]+=number*gas_stoich[j][k];
    }
    for(unsigned j=0;j<cond_stoich.size();++j) {
      const double number=out.number_densities_cond[0][j];
      if(!std::isfinite(number) || number<0)throw std::runtime_error("invalid condensate density");
      for(unsigned k=0;k<solid.size();++k)solid[k]+=number*cond_stoich[j][k];
    }
    if(!(atoms[h]>0) || std::abs(particles*1.380649e-16*T/P-1)>2e-7)
      throw std::runtime_error("depletion pressure closure failed");
    bool condensed=false;
    for(unsigned k=0;k<atoms.size();++k)if(k!=electron) {
      const double ratio=(atoms[k]+solid[k])/(atoms[h]+solid[h]);
      if(std::abs(ratio/(abundance[k]/abundance[h])-1)>2e-7)
        throw std::runtime_error("depletion element closure failed");
      condensed=condensed || solid[k]>0;
    }
    if(condensed)for(unsigned k=0;k<atoms.size();++k)if(k!=electron)
      gas[atomic_index[k]]=std::max(atoms[k]/atoms[h],1e-99);
  }
};
}
extern "C" void ember_condense_(const double* T,const double* P,const double* bulk,
                                 double* gas,const int* n,int* error) {
  try {
    if(*n!=92)throw std::runtime_error("unsupported source element count");
    static Depletion state(bulk,*n);
    state.evaluate(*T,*P,gas,*n);*error=0;
  } catch(const std::exception& e) {
    std::cerr<<"EMBER DEPLETION ERROR: "<<e.what()<<'\n';*error=1;
  }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Offline adapter for FastChem 4 (Kitzmann & Stock); not linked into Ember.
// Input: T [K], gas pressure [bar], deepest atmosphere layer first.
#include "fastchem.h"
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <stdexcept>

template<class T> void array(const std::vector<T>& v) {
  std::cout << '[';
  for (std::size_t i=0;i<v.size();++i) std::cout << (i?",":"") << v[i];
  std::cout << ']';
}
int main(int argc,char** argv) {
  try {
    if(argc!=4) throw std::invalid_argument("source-directory abundance-file gas|equilibrium|rainout");
    const std::filesystem::path source=argv[1]; const std::string mode=argv[3];
    if(mode!="gas" && mode!="equilibrium" && mode!="rainout") throw std::invalid_argument("unknown condensation mode");
    // The source writes initialization diagnostics to cout even at verbosity
    // zero. Preserve them on stderr, keeping the result stream valid JSON.
    auto* result_stream=std::cout.rdbuf(std::cerr.rdbuf());
    fastchem::FastChem chemistry(argv[2],(source/"input/logK/logK.dat").string(),
                               (source/"input/logK/logK_condensates.dat").string(),0);
    for(const auto* key:{"accuracyChem","accuracyElementConservation","accuracyCond"})
      if(!chemistry.setParameter(key,1e-8)) throw std::runtime_error("unsupported chemistry tolerance");
    fastchem::FastChemInput input; fastchem::FastChemOutput output;
    input.equilibrium_condensation=mode=="equilibrium";
    input.rainout_condensation=mode=="rainout";
    double t,p;
    while(std::cin>>t>>p) {
      if(!std::isfinite(t+p) || t<100 || t>6000 || p<1e-13 || p>1000)
        throw std::domain_error("outside documented tested FastChem T/P range");
      if(mode=="rainout" && !input.pressure.empty() && p>=input.pressure.back())
        throw std::domain_error("rainout requires decreasing pressure (bottom to top)");
      input.temperature.push_back(t); input.pressure.push_back(p);
    }
    if(!std::cin.eof() || input.temperature.empty()) throw std::invalid_argument("malformed/empty T/P input");
    const auto status=chemistry.calcDensities(input,output);
    std::cout.rdbuf(result_stream);
    std::cout<<std::setprecision(17)<<"{\"status\":"<<status<<",\"mode\":"<<std::quoted(mode)<<",\"elements\":[";
    for(unsigned j=0;j<chemistry.getElementNumber();++j)
      std::cout<<(j?",":"")<<std::quoted(chemistry.getElementSymbol(j));
    std::cout<<"],\"abundances\":"; array(chemistry.getElementAbundances());
    std::cout<<",\"gas_species\":[";
    for(unsigned j=0;j<chemistry.getGasSpeciesNumber();++j) {
      std::cout<<(j?",":"")<<"{\"symbol\":"<<std::quoted(chemistry.getGasSpeciesSymbol(j))<<",\"stoichiometry\":";
      array(chemistry.getGasSpeciesStoichiometry(j));std::cout<<'}';
    }
    std::cout<<"],\"condensates\":[";
    for(unsigned j=0;j<chemistry.getCondSpeciesNumber();++j) {
      std::cout<<(j?",":"")<<"{\"symbol\":"<<std::quoted(chemistry.getCondSpeciesSymbol(j))<<",\"stoichiometry\":";
      array(chemistry.getCondSpeciesStoichiometry(j));std::cout<<'}';
    }
    std::cout<<"],\"rows\":[";
    for(std::size_t i=0;i<input.temperature.size();++i) {
      std::cout<<(i?",":"")<<"{\"T_K\":"<<input.temperature[i]<<",\"P_bar\":"<<input.pressure[i]
               <<",\"flag\":"<<output.fastchem_flag[i]<<",\"gas\":";
      array(output.number_densities[i]); std::cout<<",\"condensed\":";array(output.number_densities_cond[i]);
      std::cout<<",\"element_conserved\":";array(output.element_conserved[i]);
      std::cout<<",\"element_condensation_degree\":";array(output.element_cond_degree[i]);
      std::cout<<",\"total_element_density\":"<<output.total_element_density[i]<<'}';
    }
    std::cout<<"]}\n";
    return status==fastchem::FASTCHEM_SUCCESS?0:2;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

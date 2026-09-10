// Replay a saved physical profile for a conduction-only timing comparison.
// This never constructs a stellar restart or advances an evolutionary model.
#include "ember/conduction_table.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <vector>

int main(int argc,char** argv) {
  try {
    if(argc!=4)throw std::invalid_argument("usage: conduction-profile TABLE PROFILE REPETITIONS");
    const int repetitions=std::stoi(argv[3]);
    if(repetitions<1 || repetitions>10000)throw std::invalid_argument("invalid repetitions");
    ember::TabulatedConduction source(argv[1]);
    struct State {double T,rho;ember::Composition c;};
    std::vector<State> states;
    std::ifstream in(argv[2]);
    double m,r,rho,T,L,X,Y3,Y4;
    while(in>>m>>r>>rho>>T>>L>>X>>Y3>>Y4) {
      auto c=ember::solar_scaled(X,.02);c.basis=ember::AbundanceBasis::baryon_mass;
      c.metal_inventory=ember::MetalInventory::gs98;c.X[1]=Y3;c.X[2]=Y4;
      // Below the production conduction join no source evaluation is needed.
      if(T>3e5)states.push_back({T,rho,c});
    }
    if(!in.eof() || states.empty())throw std::invalid_argument("invalid or empty profile");
    double checksum=0;
    for(int repeat=0;repeat<repetitions;++repeat)for(const auto& state:states) {
      const auto k=source.eval(state.T,state.rho,state.c);
      checksum+=k.kappa;
      if(repeat==0)std::printf("%.17g %.17g %.17g %.17g %.17g\n",
          k.kappa,k.dlnk_dlnT,k.dlnk_dlnRho,k.dlnk_dX,k.dlnk_dY3);
    }
    if(!std::isfinite(checksum))throw std::runtime_error("nonfinite accumulated opacity");
    std::printf("checksum %.17g\n",checksum);
  } catch(const std::exception& error) {
    std::fprintf(stderr,"%s\n",error.what());return 1;
  }
}

// Inspect the assembled opacity family through the production classes.
// Modes: high and atomic-blend use atomic X/Z, Y3=0, and atomic density.
// Stellar mode uses baryonic X/Z/Y3 and includes isotope mapping and AESOPUS.
// Input: X Z Y3 T[K] rho[g/cm3]. Unsupported states remain explicit JSON rows.
#include "ember/conduction_table.hpp"
#include "ember/opacity_mixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  using namespace ember;
  try {
    if(argc != 4) throw std::invalid_argument("usage: refractive-opacity MODE FAMILY CONDUCTION_TABLE");
    const std::string mode = argv[1];
    const std::filesystem::path family = argv[2];
    if(mode != "high" && mode != "atomic-blend" && mode != "stellar")
      throw std::invalid_argument("invalid opacity mode");
    std::unique_ptr<MixtureOpacity> low, high;
    std::unique_ptr<BlendedOpacity> blend;
    std::unique_ptr<StellarMixtureOpacity> stellar;
    const Opacity* opacity = nullptr;
    if(mode == "stellar") {
      stellar = std::make_unique<StellarMixtureOpacity>(family);
      opacity = stellar.get();
    } else {
      high = std::make_unique<MixtureOpacity>(family / "tops_gs98_mixture_high.dat");
      opacity = high.get();
      if(mode == "atomic-blend") {
        low = std::make_unique<MixtureOpacity>(family / "tops_gs98_mixture_low.dat");
        blend = std::make_unique<BlendedOpacity>(*low, *high, 5.6, 5.7);
        opacity = blend.get();
      }
    }
    TabulatedConduction conduction(argv[3]);
    HotConduction selected_conduction(conduction);
    double x, z, y3, T, rho;
    std::size_t count = 0;
    std::cout << std::setprecision(17);
    while(std::cin >> x >> z >> y3 >> T >> rho) {
      if(!std::isfinite(x+z+y3+T+rho) || x<0 || z<0 || y3<0 || x+z+y3>1 ||
         T<=0 || rho<=0 || (mode!="stellar" && y3!=0))
        throw std::invalid_argument("invalid material input");
      ++count;
      auto c = solar_scaled(x,z);
      c.X[1] = y3;
      c.X[2] -= y3;
      if(mode == "stellar") {
        c.basis = AbundanceBasis::baryon_mass;
        c.metal_inventory = MetalInventory::gs98;
      }
      std::cout << "{\"query\":[" << x << ',' << z << ',' << y3 << ',' << T << ',' << rho << ']';
      try {
        const auto value = opacity->eval(T,rho,c);
        const auto range = opacity->density_range(T,c);
        if(!range || !std::isfinite(value.kappa) || value.kappa<=0)
          throw std::runtime_error("invalid opacity response");
        auto finite = [](double v) {
          if(std::isfinite(v)) std::cout << v;
          else std::cout << "null";
        };
        std::cout << ",\"covered\":true,\"kappa\":" << value.kappa
                  << ",\"dlnk_dlnT\":" << value.dlnk_dlnT
                  << ",\"dlnk_dlnrho\":" << value.dlnk_dlnRho
                  << ",\"dlnk_dX\":";
        finite(value.dlnk_dX);
        std::cout << ",\"dlnk_dZ\":"; finite(value.dlnk_dZ);
        std::cout << ",\"dlnk_dY3\":"; finite(value.dlnk_dY3);
        std::cout << ",\"density_range\":[" << range->min << ',' << range->max << ']';
        auto cb = c;
        double mass_scale = 1.;
        if(mode != "stellar") {
          const double hs = nuclides[0].A, hes = nuclides[2].A/4, zs = gs98_atomic_mass_scale();
          mass_scale = 1/(x/hs+(1-x-z)/hes+z/zs);
          cb = solar_scaled(mass_scale*x/hs,mass_scale*z/zs);
          cb.basis = AbundanceBasis::baryon_mass;
          cb.metal_inventory = MetalInventory::gs98;
        }
        try {
          const double kc = selected_conduction.eval(T,rho/mass_scale,cb).kappa/mass_scale;
          const double combined = 1/(1/value.kappa+1/kc);
          if(!(kc>0) || !std::isfinite(combined) || combined<=0)
            throw std::runtime_error("invalid conduction response");
          std::cout << ",\"conduction_covered\":true,\"conductive_opacity\":";
          finite(kc);
          std::cout << ",\"combined_opacity\":" << combined;
        } catch(const std::domain_error& error) {
          std::cout << ",\"conduction_covered\":false,\"conduction_error\":" << std::quoted(error.what());
        }
      } catch(const std::domain_error& error) {
        std::cout << ",\"covered\":false,\"error\":" << std::quoted(error.what());
      }
      std::cout << "}\n";
    }
    if(!std::cin.eof() || count==0) throw std::invalid_argument("malformed or empty material input");
  } catch(const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

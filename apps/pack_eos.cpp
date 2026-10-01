#include "ember/eos_variable_metal.hpp"
#include <iostream>
#include <string>

int main(int argc,char** argv) {
  if(argc!=3 && !(argc==4 && std::string(argv[3])=="--zero-helium3")) {
    std::cerr<<"usage: ember-pack-eos FAMILY NEW_BINARY_FILE [--zero-helium3]\n";
    return 2;
  }
  try {
    ember::VariableMetalHelmholtzEos::pack_binary(argv[1],argv[2],argc==4);
    std::cout<<"Packed EOS family into "<<argv[2]<<'\n';return 0;
  }catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';return 1;
  }
}

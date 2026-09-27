#include "ember/eos_variable_metal.hpp"
#include <iostream>

int main(int argc,char** argv) {
  if(argc!=3) {
    std::cerr<<"usage: ember-pack-eos TEXT_FAMILY NEW_BINARY_FILE\n";
    return 2;
  }
  try {
    ember::VariableMetalHelmholtzEos::pack_binary(argv[1],argv[2]);
    std::cout<<"Packed EOS family into "<<argv[2]<<'\n';return 0;
  }catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';return 1;
  }
}

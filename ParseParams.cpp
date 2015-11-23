// Parse .yaml. Copied with permission from Jeffrey Moore

#include "ParseParams.h"

void parse_params(std::string param_file, system_params &params) {
    std::cout << "Reading parameters from " << param_file << std::endl;
    YAML::Node node = YAML::LoadFile(param_file);
        
    for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
      std::string param_name = it->first.as<std::string>();
      std::string param_value = it->second.as<std::string>();
     
      // check for all parameters
#include "ParseParamsBody.cpp"
    }
    std::cout << std::endl;
}

#include "HardRodMain.h"

int main(int argc, char *argv[]){
  
  // Read input file
  std::string Inpt = argv[1];

  //Parse parameters
  system_params params;
  params.init();
  std::string ParamsFile;
  
  //ParamsFile = "Params.yaml"; // Should be an input
  ParamsFile = Inpt + ".yaml"; // Should be an input
  parse_params(ParamsFile, params);
 
  if(params.IsoDiffFlag == 1){ HrIsoBody(params); } //Iso
  else{ HrAniBody(params); } // Aniso

  return 0;
}

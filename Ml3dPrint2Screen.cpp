// Ml3dPrint.cpp
// Print a 3D array like matlab does (:,:,0), (:,:,1) ...
// M(x,y,z)

#include "Ml3dPrint2Screen.h"

void Ml3dPrint(double*** M, int Nx, int Ny, int Nz){

for( int k = 0; k < Nz; ++k ) {
  for( int i = 0; i < Nx; ++i ) {
    for( int j = 0; j < Ny; ++j) {
      std::cout << M[i][j][k] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << std::endl;
}

}

void Ml3dPrint( double* M, int Nx, int Ny, int Nz ){
for( int k = 0; k < Nz; ++k ) {
  for( int i = 0; i < Nx; ++i ) {
    for( int j = 0; j < Ny; ++j) {
      std::cout << M[i +  j * Nx + k * Nx * Ny] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << std::endl;
}


}

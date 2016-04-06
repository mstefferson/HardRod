// CheckBrokeDen.cpp
//
// Function makes sures Density isn't negative, too large, infinite, NaN

#include "CheckBrokenDen.h"
int CheckBrokenDen( int Nx, int Ny, int Nm, Ad3 &rho, double rhoMax ){
 
  bool ShitIsFucked = 0;

  for( int i = 0; i < Nx; ++i ) {
    for( int j = 0; j < Ny; ++j ) {
      for( int k = 0; k < Nm; ++k ) {

        // Negative
        if( (rho[i][j][k] < 0) ) { 
          ShitIsFucked = 1;
          std::cout << " Density is negative " << std::endl;
            }
            
        // Nan
        if( (rho[i][j][k] != rho[i][j][k]) ) {
          ShitIsFucked = 1;
          std::cout << " Density is NaN " << std::endl;
        }

        // Infinit
        if( ( std::isinf(rho[i][j][k]) == 1 ) ){
          ShitIsFucked = 1;
          std::cout << " Density if Inf " << std::endl;
         }

        // Too large
            
        if( ( rho[i][j][k] > rhoMax ) ){
          ShitIsFucked = 1;
          std::cout << " Density too large " << std::endl;
        }

        if( ShitIsFucked == 1 ) { break;}
      }
      if( ShitIsFucked == 1 ) { break;}
    }
    if( ShitIsFucked == 1 ) { break;}
  }
  return ShitIsFucked;
}




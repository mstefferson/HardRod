#include <iostream>
//#include "Mayer.h"
#include <math.h>
#include "Array.h"
#include "fftw++.h"
#include "MayerFnc.h"



int main(){

  int Nx = 8;
  int Ny = 8;
  int Nm = 8;

  double Lrod, Lx,Ly;

  Lrod = sqrt(2.0) * 2.0;
  Lx   = 10.0;
  Ly   = 10.0;

  //int MayerFnc[Nx][Ny][Nm];
  fftwpp::fftw::maxthreads=get_max_threads(); //Multithreads for fft
  
  //Size of vectors to fft
  unsigned int Nkx,Nky,Nkm;
  size_t align=sizeof(Complex);

  //Build density matrix
  Nkx = Nx/2 + 1; Nky = Ny/2 +1; Nkm = Nm / 2 + 1;

  Array::array3<int> MayerFnc( Nx, Ny, Nm,align);
  Array::array3<Complex> MayerFncFT( Nx, Ny,Nkm,align);

  MayerFncHardRod(Nx, Ny,Nm, Lx,Ly, Lrod, MayerFnc);


 // Print it
 double dphi = 2 * M_PI / Nm;
  for( int k  = 0; k < Nm; ++k ){
    std::cout << " k = " << k << " phi =  " << k* dphi << std::endl;
    for( int i = 0; i < Nx; ++i ) {
      for( int j = 0; j < Ny; ++j ) {
    std::cout << MayerFnc[i][j][k] << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
  }
  return 0;
}



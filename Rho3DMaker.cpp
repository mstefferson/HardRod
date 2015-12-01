// Rho3DMaker.cpp

#include "Rho3DMaker.h"

void Rho3DMaker::BuilderEq( int Nx, int Ny, int Nm, double Lx, double Ly, double c, 
        Array::array3<double> &rho, Array::array1<double> &f ){

  double Norm = c / (Lx * Ly);

  for(int i  = 0; i < Nx; ++i){
    for(int j = 0; j < Ny; ++j){
      for(int k =0; k < Nm; ++k){

        rho[i][j][k] = Norm * f[k];

      }
    }
  }
}


void Rho3DMaker::BuilderIs( int Nx, int Ny, int Nm, double Lx, double Ly, double c, 
        Array::array3<double> &rho ){

  double Norm = c / (Lx * Ly);

  for(int i  = 0; i < Nx; ++i){
    for(int j = 0; j < Ny; ++j){
      for(int k =0; k < Nm; ++k){

        rho[i][j][k] = Norm / ( 2 * M_PI );

      }
    }
  }

}
/*
void Rho3DMaker::BuilderEqFT( int Nx, int Ny, int Nm, int Nkm, double Lx, double Ly, double c, 
        Array::array3<Complex> &rhoFT, Array::array1<Complex> &fFT ){

  double Norm = c / (Lx * Ly);

  for(int i  = 0; i < Nx; ++i){
    for(int j = 0; j < Ny; ++j){
      for(int k =0; k < Nkm; ++k){

        rhoFT[i][j][k] = Norm * fFT[k];

      }
    }
  }
}

void Rho3DMaker::BuilderEqFT( int Nx, int Ny, int Nm, int Nkm, double Lx, double Ly, double c, 
        Array::array3<Complex> &rhoFT ){

  double Norm = c / (Lx * Ly);

  for(int i  = 0; i < Nx; ++i){
    for(int j = 0; j < Ny; ++j){
      for(int k =0; k < Nkm; ++k){

        rhoFT[i][j][k] = Norm * fFT[k]

      }
    }
  }
  rhoFT[0][0][0] = Norm * Nx * Ny * Nm;
}
*/

void Rho3DMaker::Perturb(int Nx, int Ny, int Nm,  int pX, int pY, int pM, double PertrbAmp,
        Array::array3<Complex> &rhoFT){

  for( int i = 0; i <= pX; ++i ) {
    for( int j = 0; j <= pY; ++j ) {
      for( int k = 0; k <= pM; ++k ) {

        if( i + j + k !=0 ){
        rhoFT[i][j][k] += PertrbAmp * Nx * Ny * Nm;
        }

      }
    }
  }
}

void Rho3DMaker::print(int Nx, int Ny, int Nm, Array::array3<double> &rho){
 double dphi = 2 * M_PI / Nm;
  for( int k  = 0; k < Nm; ++k ){
    std::cout << " k = " << k << " phi =  " << k* dphi << std::endl;
    for( int i = 0; i < Nx; ++i ) {
      for( int j = 0; j < Ny; ++j ) {
    std::cout << rho[i][j][k] << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
  }


}

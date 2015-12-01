#include <iostream>
//#include "Mayer.h"
#include <math.h>
#include "Array.h"
#include "fftw++.h"



void MayerFncHardRod(int Nx,int Ny,int Nm,double Lx,double Ly,double Lrod,
    Array::array3<int> MayerFnc);

//void MayerFncHardRod(int Nx,int Ny,int Nm,double Lx,double Ly,double Lrod,
   // int ***MayerFnc);

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
   
  double xTemp, yTemp, phiTemp, xTempSq, yTempSq, ylL, yuL, yMaxSq, DistSq;
  double dx    = Lx/Nx;
  double dy    = Ly/Ny;
  double dphi  = 2*M_PI / Nm;
  double epsilon = 0.00001;
  double TooFar = Lrod * Lrod + epsilon;
  double LrodHSq = TooFar  / 4;

  for( int i = 0; i < Nx ; ++i ){   
    if( i <= Nx/2 ) { xTemp = i*dx;}
    else{ xTemp = ( - Nx + i) * dx; }
    xTempSq = xTemp * xTemp;
    
    for( int j = 0; j < Ny; ++j ){
      if( j <= Ny/2 ) { yTemp = j*dy;}
      else{ yTemp = ( - Ny + j ) * dy; }
      yTempSq = yTemp * yTemp;
      DistSq = yTempSq + xTempSq;
      
      for( int k = 0; k < Nm; ++k ) {
        //std::cout << " i = " << i << " j = " << j << " k = " << k << std::endl;
        //std::cout << " x = " << xTemp << " y = " << yTemp << " phi = " << k * dphi << std::endl;
        
        if( DistSq <= TooFar ){        
          phiTemp = k * dphi;
          
          // phi = 0,pi
          if( k == 0 || k == Nm / 2 ) {
             if( (xTempSq <= TooFar) && yTemp == 0 ) { MayerFnc[i][j][k] = -1; }
             else{ MayerFnc[i][j][k] = 0; }
          }
          
          // phi = pi/2, 3pi/2
          if( k == Nm / 4 || k == 3 * Nm / 4 ){
//            std::cout << "  phi = " << phiTemp << " y = " << yTemp << " x = " << xTemp << std::endl;
//            std::cout << " i = " << i << " j = " << j << " k = " << k << std::endl;
             if( (yTempSq <= LrodHSq) && xTempSq <= LrodHSq) { MayerFnc[i][j][k] = -1; }
             else{ MayerFnc[i][j][k] = 0; }
//             std::cout << " Fm = " << MayerFnc[i][j][k] << std::endl << std::endl;
          } 
          
          // 0 < phi < pi/2
          if( k < floor( Nm / 4) ){
            yuL = tan( phiTemp ) * ( xTemp + Lrod / 2) + epsilon;
            ylL = tan( phiTemp ) * ( xTemp - Lrod / 2) - epsilon;
            yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp );

            if( (yTemp >= ylL) && (yTemp <= yuL) 
                && (yTempSq <= yMaxSq) ){ MayerFnc[i][j][k] = - 1; }
            else{ MayerFnc[i][j][k] = 0; }
           }

          // pi/2  < phi < pi
          if( k > floor( Nm / 4) && k < floor(Nm/2) ){
            yuL = tan( phiTemp ) * ( xTemp - Lrod / 2) + epsilon;
            ylL = tan( phiTemp ) * ( xTemp + Lrod / 2) - epsilon;
            yMaxSq =  LrodHSq *  sin( phiTemp ) * sin( phiTemp );

            if( (yTemp >= ylL) && (yTemp <= yuL) 
                && (yTempSq <= yMaxSq) ){ MayerFnc[i][j][k] = - 1; }
            else{ MayerFnc[i][j][k] = 0; }
          }
          
          // pi < phi < 3pi/2
          if( k > floor(Nm/2) && k < floor(3*Nm / 4) ){
            yuL = tan( phiTemp ) * ( xTemp + Lrod / 2) + epsilon;
            ylL = tan( phiTemp ) * ( xTemp - Lrod / 2) - epsilon;
            yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp);
            
            if( (yTemp >= ylL) && (yTemp <= yuL) 
                && (yTempSq <= yMaxSq) ){ MayerFnc[i][j][k] = - 1; }
            else{ MayerFnc[i][j][k] = 0; }
           }

          // 3pi/2 < phi < 2pi
          if(  k > floor(3*Nm / 4)  ) {
            yuL = tan( phiTemp ) * ( xTemp - Lrod / 2) + epsilon;
            ylL = tan( phiTemp ) * ( xTemp + Lrod / 2) - epsilon;
            yMaxSq =  LrodHSq * sin( phiTemp ) * sin( phiTemp );

            if( (yTemp >= ylL) && (yTemp <= yuL) 
                && (yTempSq <= yMaxSq) ){ MayerFnc[i][j][k] = - 1; }
            else{ MayerFnc[i][j][k] = 0; }
          }    
        } //end if dist too far         
        
        else { MayerFnc[i][j][k] = 0; }
         

        if( k == 1 ){
            std::cout << " x = " << xTemp << " y = " << yTemp << " phi = " << phiTemp << std::endl;
            std::cout << " xSqr = " << xTempSq << " ySqr = " << yTempSq << " DistSq = " << DistSq << std::endl;
            std::cout << " i = " << i << " j = " << j << " k = " << k << std::endl;
            std::cout << " yuL = " << yuL << " ylL = " << ylL << " yMaxSq = " << yMaxSq << 
              " Too Far = " << TooFar << " LrodHSq = " << LrodHSq << std::endl;
            std::cout << " Fm = "  << MayerFnc[i][j][k] << std::endl << std::endl;
        }

        // std::cout << " Fm = " << MayerFnc[i][j][k] << std::endl << std::endl;
      } //end k loop
    } //end y loop
  } //end x loop

 // Print it 
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



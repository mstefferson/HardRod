// MayerFnc
#include "MayerFnc.h"

void MayerFncHardRod(int Nx, int Ny,int Nm, double Lx,double Ly, double Lrod,
    Ad3 &MayerFnc){
 
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
        
        if( DistSq <= TooFar ){        
          phiTemp = k * dphi;
          
          // phi = 0,pi
          if( k == 0 || k == Nm / 2 ) {
             if( (xTempSq <= TooFar) && yTemp == 0 ) { MayerFnc[i][j][k] = -1; }
             else{ MayerFnc[i][j][k] = 0; }
          }
          
          // phi = pi/2, 3pi/2
          if( k == Nm / 4 || k == 3 * Nm / 4 ){
             if( (yTempSq <= LrodHSq) && xTempSq <= LrodHSq) { MayerFnc[i][j][k] = -1; }
             else{ MayerFnc[i][j][k] = 0; }
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
         
      } //end k loop
    } //end y loop
  } //end x loop
} // end fnc

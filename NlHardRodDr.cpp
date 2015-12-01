#include "NlHardRodDr.h"

NlHardRodDr::NlHardRodDr(int Nx, int Ny, int Nm, double* kx, double* ky, double* km, double Lx, double Ly){

  Nx_ = Nx;
  Ny_ = Ny;
  Nm_ = Nm;
  Nkm_ = Nm / 2 + 1;
  
  i_ = 1i;

  IntFlag_ = 1;
  IsoDiffFlag_ = 0;

  kx_ = kx;
  ky_ = ky;
  km_ = km;

  ConvCoeff_ = 2 * M_PI * Lx * Ly / (Nx * Ny * Nm) ;

}


///////////////////////// NL interactions ////////////////////////////////
void NlHardRodDr::NLinit( Ac3 &NlFT ){

  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nkm_; ++kk ){
        
        NlFT[ii][jj][kk] = 0;

      }
    }
  }
}

void NlHardRodDr::NLIntCalcC( Ac3 &rhoFT, Ad3 &rho, Ac3 &FmFT, 
    Ac3 &MuExFT, Ad3 &diMuTemp, Ac3 &diMuTempFT,
    Ad3 &ji, Ac3 &jiFT, Ac3 &NlFT, fftwpp::rcfft3d &Forward3, fftwpp::crfft3d &Backward3){

  for( int ii = 0; ii < Nx_; ++ii){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nkm_; ++kk ){

        MuExFT[ii][jj][kk] = -ConvCoeff_ * FmFT[ii][jj][kk] * rhoFT[ii][jj][kk];

      }
    }
  }

  // Flux and NL contribution from x

  for( int hh = 0; hh < 3; ++hh) {
    
    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nkm_; ++kk ){

          switch( hh ) {
            case 0: diMuTempFT[ii][jj][kk] = i_ * kx_[ii] * MuExFT[ii][jj][kk];
            case 1: diMuTempFT[ii][jj][kk] = i_ * ky_[jj] * MuExFT[ii][jj][kk];
            case 2: diMuTempFT[ii][jj][kk] = i_ * km_[kk] * MuExFT[ii][jj][kk];
          }

        }
      }
    }

    Backward3.fftNormalized( diMuTempFT, diMuTemp );
 

    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nm_; ++kk ){

          ji[ii][jj][kk] = rho[ii][jj][kk] * diMuTemp[ii][jj][kk];

        }
      }
    }

    Forward3.fft( ji, jiFT );

    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nkm_; ++kk ){

          switch( hh ) {
            case 0: NlFT[ii][jj][kk] = i_  * kx_[ii] * jiFT[ii][jj][kk];
            case 1: NlFT[ii][jj][kk] = NlFT[ii][jj][kk] + i_ * ky_[jj] * jiFT[ii][jj][kk];
            case 2: NlFT[ii][jj][kk] = NlFT[ii][jj][kk] + i_ * km_[kk] * jiFT[ii][jj][kk];
          }

        }
      }
    }

  } // loop over dimensions

} // Function




void NlHardRodDr::NlDrCalcC( Ad3 &ji, Ac3 &jiFT, Ac3 &NlFT ){

  //Driving not written yet

  if( IntFlag_ == 0 ){
    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nkm_; ++kk ){
          
          NlFT[ii][jj][kk] = 0;

        }
      }
    }
  }

}






#include "NlHardRodDr.h"

NlHardRodDr::NlHardRodDr(int Nx, int Ny, int Nm, double* kx, double* ky, double* km, double Lx, double Ly, double* phi, double vD){

  Nx_ = Nx;
  Ny_ = Ny;
  Nm_ = Nm;
  Nkm_ = Nm / 2 + 1;
  
  ConvCoeff_ = 2 * M_PI * Lx * Ly / (Nx * Ny * Nm) ;

  i_ = 1i;

  IntFlag_ = 1;
  IsoDiffFlag_ = 0;
 
  ikx_ = new Complex[Nx];
  iky_ = new Complex[Ny];
  ikm_ = new Complex[Nm];
   
  for( int i = 0; i < Nx; ++i ){ ikx_[i] = i_ * kx[i]; }
  for( int i = 0; i < Ny; ++i ){ iky_[i] = i_ * ky[i]; }
  for( int i = 0; i < Nm; ++i ){ ikm_[i] = i_ * km[i]; }

  if( vD != 0 ){
    vDcos_ =  new double[Nm];
    vDsin_ =  new double[Nm];
    for( int i = 0; i < Nm; i++ ){
      vDcos_[i] = vD * cos( phi [i] );
      vDsin_[i] = vD * sin( phi [i] );
    }
  }
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
            case 0: 
              //diMuTempFT[ii][jj][kk] = i_ * kx_[ii] * MuExFT[ii][jj][kk];
              diMuTempFT[ii][jj][kk] = ikx_[ii] * MuExFT[ii][jj][kk];
              break;
            case 1: 
              //diMuTempFT[ii][jj][kk] = i_ * ky_[jj] * MuExFT[ii][jj][kk];
              diMuTempFT[ii][jj][kk] = iky_[jj] * MuExFT[ii][jj][kk];
              break;
            case 2: 
              //diMuTempFT[ii][jj][kk] = i_ * km_[kk] * MuExFT[ii][jj][kk];
              diMuTempFT[ii][jj][kk] = ikm_[kk] * MuExFT[ii][jj][kk];
              break;
          }
        }
      }
    }

    Backward3.fftNormalized( diMuTempFT, diMuTemp );

// Calculate flux. Really, -flux since sign will cancel
    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nm_; ++kk ){

          ji[ii][jj][kk] = rho[ii][jj][kk] *  diMuTemp[ii][jj][kk];

        }
      }
    }

    Forward3.fft( ji, jiFT );
// Calculate -Grad of flux. Really, grad since sign cancels.
    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nkm_; ++kk ){

          switch( hh ) {
            case 0: 
              NlFT[ii][jj][kk] = ikx_[ii] * jiFT[ii][jj][kk];
              //NlFT[ii][jj][kk] = i_  * kx_[ii] * jiFT[ii][jj][kk];
              break;
            case 1: 
              NlFT[ii][jj][kk] = NlFT[ii][jj][kk] + iky_[jj] * jiFT[ii][jj][kk];
              //NlFT[ii][jj][kk] = NlFT[ii][jj][kk] + i_ * ky_[jj] * jiFT[ii][jj][kk];
              break;
            case 2: 
              NlFT[ii][jj][kk] = NlFT[ii][jj][kk] + ikm_[kk] * jiFT[ii][jj][kk];
              //NlFT[ii][jj][kk] = NlFT[ii][jj][kk] + i_ * km_[kk] * jiFT[ii][jj][kk];
              break;
          }

        }
      }
    }

  } // loop over dimensions

} // Function


void NlHardRodDr::NlDrCalcC( Ad3 &rho, 
    Ad3 &ji, Ac3 &jiFT, Ac3 &NlFT, fftwpp::rcfft3d &Forward3){

  //Driving not written yet

    for( int ii = 0; ii < Nx_; ++ii ){
      for( int jj = 0; jj < Ny_; ++jj ){
        for( int kk = 0; kk < Nm_; ++kk ){
          
          ji[ii][jj][kk] = vDcos_[kk] * rho[ii][jj][kk];

        }
      }
    }

  Forward3.fft(ji,jiFT);

  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nkm_; ++kk ){
        
        NlFT[ii][jj][kk] += ikx_[ii] * jiFT[ii][jj][kk];

      }
    }
  }      

  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nm_; ++kk ){
        
        ji[ii][jj][kk] = vDsin_[kk] * rho[ii][jj][kk];

      }
    }
  }

  Forward3.fft(ji,jiFT);

  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nkm_; ++kk ){
        
        NlFT[ii][jj][kk] += iky_[ii] * jiFT[ii][jj][kk];

      }
    }
  }      
}







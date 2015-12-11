void NLIntCalcC( Ac3 &rhoFT, Ad3 &rho, Ac3 &FmFT, Ac3 &MuExFT, Ac3 &diMuTemp, Ac3 &diMuTempFT,
    Ad3 &ji, Ac3 &jiFT, Ac3 NlFT, fftwpp::rcfft3d &Forward3, fftwpp::crfft3d Backward3){

for( int ii = 0; ii < Nx_; ++ii ){
  for( int jj = 0; jj < Ny_; ++jj ){
    for( int kk = 0; kk < Nkm_; ++kk ){

      MuExFT[ii][jj][kk] = ConvCoeff_ * FmFT[ii][jj][kk] * rhoFT[ii][jj][kk]

    }
  }
}

Backward3.fftNormalized( MuExFT, MuEx );

// Flux and NL contribution from x

for( int hh = 0; hh < 3, ++h) {

  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nkm_; ++kk ){

        switch( hh ) {
          case 0: diMuTempFT[ii][jj][kk] = 1i * kx[ii] * MuExFT[ii][jj][kk];
          case 1: diMuTempFT[ii][jj][kk] = 1i * ky[jj] * MuExFT[ii][jj][kk];
          case 2: diMuTempFT[ii][jj][kk] = 1i * km[kk] * MuExFT[ii][jj][kk];
        }

      }
    }
  }

  Backward3.fftNormalized( diMuTempFT, diMuTemp );
  
  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nm_; ++k ){

        ji[ii][jj][kk] = rho[ii][jj][kk] * diMuTemp[ii][jj][kk];

      }
    }
  }

  Forward3.fft( ji, jiFT );

  for( int ii = 0; ii < Nx_; ++ii ){
    for( int jj = 0; jj < Ny_; ++jj ){
      for( int kk = 0; kk < Nkm_; ++kk ){

        switch( hh ) {
          case 0: NlFT[ii][jj][kk] = 1i  * kx_[ii] * jiFT[ii][jj][kk];
          case 1: NlFT[ii][jj][kk] += 1i * ky_[jj] * jiFT[ii][jj][kk];
          case 2: NlFT[ii][jj][kk] += 1i * km_[kk] * jiFT[ii][jj][kk];
        }

      }
    }
  }

} // loop over dimensions


} // Function




// NlHardRodDr.h

#ifndef _HARDROD_NLHARDRODDR_H_
#define _HARDROD_NLHARDRODDR_H_

#include <cstdlib>
#include <math.h>
#include <complex>
#include "fftw++.h"
#include "Array.h"
#include "typedef.h"
#include "stdio.h"

class NlHardRodDr{

  private:
    int Nx_;
    int Ny_;
    int Nm_;
    int Nkm_;
    
    bool IntFlag_;
    bool IsoDiffFlag_;

    double ConvCoeff_;

    double* vDcos_;
    double* vDsin_;

    Complex* ikx_;
    Complex* iky_;
    Complex* ikm_;

  public:
    // Constructor
    NlHardRodDr(int Nx, int Ny, int Nm, double* kx, double* ky, double* km, double Lx, double Ly, double* phi, double vD);
   
    // Functions
    void NLIntCalcC( Ac3 &rhoFT, Ad3 &rho, Ac3 &FmFT, Ac3 &MuExFT, Ad3 &diMuTemp, Ac3 &diMuTempFT,
    Ad3 &ji, Ac3 &jiFT, Ac3 &NlFT, fftwpp::rcfft3d &Forward3, fftwpp::crfft3d &Backward3);

    void NLinit( Ac3 &NlFT );

    void NlDrCalcC( Ad3 &rho, 
    Ad3 &ji, Ac3 &jiFT, Ac3 &NlFT, fftwpp::rcfft3d &Forward3);


    // NL calc
};

#endif

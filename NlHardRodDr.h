// NlHardRodDr.h

#include <cstdlib>
#include <math.h>
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
    Complex i_;

    bool IntFlag_;
    bool IsoDiffFlag_;

    double ConvCoeff_;
    double* kx_;
    double* ky_;
    double* km_;

  public:
    // Constructor
    NlHardRodDr(int Nx, int Ny, int Nm, double* kx, double* ky, double* km, double Lx, double Ly);
   
    // Functions
    void NLIntCalcC( Ac3 &rhoFT, Ad3 &rho, Ac3 &FmFT, Ac3 &MuExFT, Ad3 &diMuTemp, Ac3 &diMuTempFT,
    Ad3 &ji, Ac3 &jiFT, Ac3 &NlFT, fftwpp::rcfft3d &Forward3, fftwpp::crfft3d &Backward3);

    void NLinit( Ac3 &NlFT );

    void NlDrCalcC( Ad3 &ji, Ac3 &jiFT, Ac3 &NL );


    // NL calc
};

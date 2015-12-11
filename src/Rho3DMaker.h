// Build and perturb a 3D from a 1D distribution

#ifndef _HARDROD_RHO3DMAKER_H_
#define _HARDROD_RHO3DMAKER_H_


#include "Array.h"
#include "math.h"
#include "fftw++.h" 
#include "typedef.h"
#include <iostream>
#include <math.h>

class Rho3DMaker{
  public:
    
    void BuilderEq( int Nx, int Ny, int Nm, double Lx, double Ly, double c, 
        Array::array3<double> &rho, Array::array1<double> &f );
    
    void BuilderIs( int Nx, int Ny, int Nm, double Lx, double Ly, double c, 
        Array::array3<double> &rho );
    
/*    void BuilderEqFT( int Nx, int Ny, int Nm, int Nkm, double Lx, double Ly, double c, 
        Array::array3<double> &rhoFT, Array::array1<double> &fFT );
    
    void BuilderEqFT( int Nx, int Ny, int Nm, int Nkm, double Lx, double Ly, double c, 
        Array::array3<double> &rhoFT);
  */  
    void Perturb( int Nx, int Ny, int Nm, int pX, int pY, int pM,
        double PertrbAmp, Array::array3<Complex> &rhoFT );
    
    void PerturbReal(int Nx, int Ny, int Nm, int pX, int pY, int pM,
        double* kx , double* ky, double* km, double* x,double *y, double* phi, 
        double PertrbAmp, Ad3 &rho);
 
    void print(int Nx, int Ny, int Nm, Array::array3<double> &rho);

};

#endif

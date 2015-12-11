//EqDist.h
// All things equilbrium density
#ifndef _HARDROD_EQDIST_H_
#define _HARDROD_EQDIST_H_

#include <iostream>
#include <cstdlib>
#include <math.h>
#include "spGrid.h"
#include "Array.h"

class EqDist{

  private:
    int N_; // Number of gridpoints
    int Nc_; // Number of Coeff
    double bc_; // scaled concentration
    double* phi_;
    double* feq_; // equilbrium distribution
    double* fis_; // isotropic ditribution

  public:
    // Constructor
    EqDist();
    EqDist(int N, double bc);

    double trapz(double* f, double* x, int N);
    double trapz_periodic(double* f, double* x, int N);
    void fisInit( int N, double* fis );
    void fisArrayInit( int N, Array::array1<double> &fis );
    void feqInit( int N, double bc, double *f); 
    void KernCoeffCalcHardRod2D( double* d2nVec, int Nc );
    void BestCoeffExpLeg2D( double* Coeff, int Nc, double* phi, double* bc );
    double* getfisVec();

};

#endif

//EqDist.h
// All things equilbrium density
#ifndef _HARDROD_EQDIST_H_
#define _HARDROD_EQDIST_H_

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include "spGrid.h"
#include "Array.h"
#include "typedef.h"

class EqDist{

  private:
    int Nm_; // Number of gridpoints
    int Nc_; // Number of Coeff
    double bc_; // scaled concentration
    bool DEBUG_;

    double* Coeff_;
    double* d2nVec_;
    double* phi_;
    double* feq_; // equilbrium distribution
    double* fis_; // isotropic ditribution

  public:
    // Constructor
    EqDist();
    EqDist(int N, double bc);

    double trapz_periodicPhi(const double f[]  );
    void Normalize(double f[]);
    void bcFix();
    
    void fisInit( );
    void feqInit();
    
    void fisInit( Ad1 &fint );
    void feqInit( Ad1 &fint); 
    
    void KernCoeffCalcHardRod2D( );
    void BestCoeffExpLeg2D( );
    void DistBuilderExpCos2Dsing();
    double MaxMagElement(double v[], int size);

    double* getfisVec();
    void printFis();
    void printFeq();

};

#endif

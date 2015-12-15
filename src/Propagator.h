// Propagate.h
#ifndef _HARDROD_PROPAGATE_H_
#define _HARDROD_PROPAGATE_H_
#include <iostream>
#include <math.h>
#include "Array.h"
#include "fftw++.h"
#include "typedef.h"

class Propagator
{
  public:
     
    Propagator();
    Propagator(int N1,double* k1, double D, double dt, bool IsoFlag);
    Propagator(int N1, int N2, double* k1, double* k2,double D, double dt, bool IsoFlag);
    Propagator(int N1, int N2, int N3, double* k1, double* k2, double* k3, double D,
               double dt, bool IsoFlag);
    Propagator( int N1, int N2, int N3, double* k1, double* k2, double* k3,
                         double Dpos, double Dr, double dt, bool Iso);
    // Constructor
    
    void LopIsoDiagMaker(double*  k1t, double D);
    void LopIsoDiagMaker(double*  k1t, double* k2t, double* k3t, double D);
    void LopIsoDiagMaker(double* k1t, double* k2t, double* k3t, double Dpos, double Dr);

    //Build 
    
    void PropIsoMaker1(double dt);
    void PropIsoMaker2(double dt);
    void PropIsoMaker3(double dt);
  
    //Build isotropic diffusion Propagator
  
    double* getLop();
    double* getProp();
    
    void PropPrint();
    
    //void PropMaster( int ABflag, Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NL, Ac3 &NLprev);
    
    void PropAB0c( Ac3 &rhoFTnext, Ac3 &rhoFT); // Diffusion
    void PropAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NL ); // AB 1
    void PropAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NL, Ac3 &NLprev); // AB 2

  private:
    bool IsoFlag_;
    int N1_;
    int N2_;
    int N3_;
    double dt_;
    double* LopFTd_; // Diagonal of operator
    double* Ud_;    // Diagonal of propator
};

#endif //_DIFFUSEFT_LOP_H_

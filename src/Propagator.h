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
  
    // Prefactor initialization

    void PreFacInit();
    //Build isotropic diffusion Propagator
  
    double* getLop();
    double* getProp();
    
    void PropPrint();
    
    void PropMaster(Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NL);
    void PropMaster(Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NL, Ac3 &NLprev);
    
    void PropAB0c( Ac3 &rhoFTnext, Ac3 &rhoFT); // Diffusion
    void PropAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT ); // AB 1
    void PropAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev); // AB 2
    void PropHAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT); // Hybrid AB 1
    void PropHAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev); // Hybrid AB 2
    void PropBHAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT); // Better Hyrid AB 2
    void PropBHAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev); //Better Hyrid AB 2

  private:

    int IsoFlag_;
    int StepFlag_;

    int N1_;
    int N2_;
    int N3_;

    double dt_;
    double NlPf_;
    double NlPfPrev_;
    double NlPfExp_;
    double* LopFTd_; // Diagonal of operator
    double* Ud_;    // Diagonal of propator
};

#endif //_DIFFUSEFT_LOP_H_

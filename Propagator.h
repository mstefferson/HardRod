// Propagate.h
#ifndef _DIFFUSEFT_LOP_H
#define _DIFFUSEFT_LOP_H
#include <iostream>
#include <math.h>

class Propagator
{
  public:
     
    Propagator();
    Propagator(int N1,double* k1, double D, double dt, bool IsoFlag);
    Propagator(int N1, int N2, double* k1, double* k2,double D, double dt, bool IsoFlag);
    Propagator(int N1, int N2, int N3, double* k1, double* k2, double* k3, double D,
               double dt, bool IsoFlag);
    // Constructor
    
    void LopDiagMaker(double*  k1t, double D);
    void LopDiagMaker(double*  k1t, double* k2t, double* k3t, double D);
    //Build 
    
    void PropIsoMaker1(double dt);
    void PropIsoMaker2(double dt);
    void PropIsoMaker3(double dt);
    //Build isotropic diffusion Propagator
  
    double* getLop();
    double* getProp();
    void PropPrint();
  private:
    bool IsoFlag_;
    int N1_;
    int N2_;
    int N3_;
    double* LopFTd_; // Diagonal of operator
    double* Ud_;    // Diagonal of propator
};

#endif //_DIFFUSEFT_LOP_H_

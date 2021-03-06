#ifndef _HARDROD_PARAMETERS_H_
#define _HARDROD_PARAMETERS_H_

struct system_params{

  int seed;

  int trial,
      Nx,
      Ny,
      Nm,
      pXmodes,
      pYmodes,
      pMmodes,
      StepFlag, 
      IsoDiffFlag,
      IcFlag,
      RandPerbAmpFlag,
      wrtOpFlag,
      wrtRhoFlag;

       
  double dt,
         trec,
         tend,
         bc,
         vD,         
         Lx,
         Ly,
         Lrod,
         Dpar,
         Dperp,
         Dr,
         PertrbAmp;

  //Initialize

  inline void init(){
#include "ParametersDefault.h"
  }

};

#endif //_HARDROD_PARAMETERS_H_

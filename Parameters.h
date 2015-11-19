#ifndef _HARDROD_PARAMETERS_H_
#define _HARDROD_HARDROD_H_

struct system_params{

  int seed;

  int trial,
      Nx,
      Ny,
      Nm,
      IsoDiffFlag,
      EqIcFlag,
      RandPerbAmpFlag,
      pXmodes,
      pYmodes,
      pMmodes;
      

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
         PerbAmp;

  //Initialize

  inline void init(){
#include "ParametersDefault.h"
  }

};

#endif //_HARDROD_PARAMETERS_H_

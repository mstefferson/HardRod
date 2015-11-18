#ifndef _HARDROD_PARAMETERS_H_
#define _HARDROD_HARDROD_H_

struct system_params{

  long seed;

  int trial,
      Nx,
      Ny,
      Nm,
      IsoDiffFlag,
      EqIcFlag,

  double dt,
         trec,
         tot,
         bc,
         vD,         
         Lx,
         Ly,
         Lrod,
         Dpar,
         Dperp,
         Dr;

  //Initialize

  inline void init(){
#include "ParametersDefault.h"
  }

};

#endif //_HARDROD_PARAMETERS_H_

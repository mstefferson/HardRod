// MayerFnc.h

#ifndef _HARDROD_MAYER_H_
#define _HARDROD_MAYER_H_

#include <math.h>
#include "Array.h"
#include "fftw++.h"
#include "typedef.h"

void MayerFncHardRod(int Nx,int Ny,int Nm,double Lx,double Ly,double Lrod,
    Ad3 &MayerFnc);

#endif 

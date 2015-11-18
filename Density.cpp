// Density.cpp

#include "Density.h"


void Density::DenInt(int N, double f[], double x[], double k[] )
{


  for( int i = 0; i < N; ++i )
  {
    f[i] = 2;
    for( int j = 1; j <= 2; ++j )
    {
    f[i] = f[i]  + 0.1 * cos( k[j] * x[i] );
    }
  }

}

  
void Density::fisoStepper(double* f, double* U, int N)
{

  for(int i = 0; int i < N; ++i )
  {
    f[i] = U[i] * f[i];
  }

}


// Density.h
#include "math.h"
#include "Array.h"

class Density
{

  public:
    
    void DenInt( int N, double f[], double x[], double k[]);
    void fisoStepper(Complex f[], double U[], int N);
    void fStepperNL();
    void NlCalc();
    void recorder();
    
  };


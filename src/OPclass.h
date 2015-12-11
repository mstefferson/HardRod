// OPclass.h

#include <math.h>
#include "Array.h"
#include "typedef.h"

class OPclass{

  private:


    int Nx_;
    int Ny_;
    int Nm_;

    double intFac_;
    double nxTemp_;
    double nyTemp_;
    double QxxTemp_;
    double QxyTemp_;
    double QyyTemp_;

  
    double* sin_;
    double* cos_;
    double* QxxInt_;
    double* QyyInt_;
    double* cossin_;

    double** C_;
    double** NO_;
    double** PO_;


  public:
    
    OPclass( int Nx, int Ny, int Nm, double* phi );
    void OPmaker(Ad3& rho);
    void ConcCalc(Ad3& rho);
    void PolarOrdCalc(Ad3& rho);
    void NemOrdCalc(Ad3& rho);
    void printC();
    void printPO();
    void printNO();
    void printTrigs();
};

// All things equilbrium density
 #include "EqDist.h"

EqDist::EqDist(){
  
  int N_ = 0; // Number of gridpoints
  int Nc_ = 0; // Number of Coeff
  double bc_ = 0.0; // scaled concentration
  
  phi_ = NULL;
  fis_ = NULL;
  feq_ = NULL; // equilbrium distribution

}

EqDist::EqDist(int N, double bc){

  int N_ = N; // Number of gridpoints
  int Nc_ = 10; // Number of Coeff
  double bc_ = bc; // scaled concentration
  
  spGrid phi;
  phi_ = new double[N_];
  phi.xVecMaker(N_,2*M_PI,phi_);

  fis_ = new double[N_];
  
  double* feq_ = NULL; // equilbrium distribution
  fisInit(N_,fis_); // isotropic 
}


void EqDist::fisInit( int N, double* fis ){
 
  for( int i = 0; i < N; ++i ){
    fis[i] = 1 / (2 * M_PI);
  }

}

void EqDist::fisArrayInit( int N, Array::array1<double> &fis ){
 
  for( int i = 0; i < N; ++i ){
    fis[i] = 1 / (2 * M_PI);
  }

}

void EqDist::feqInit( int N, double bc, double *f ){

  std::cout << " Not written yet " << std::endl;

}

//KernCoeffCalcHardRod2D.m
// Michael Stefferson
// Calculate the coefficients of the Legendre expansion of the kernal
// |sin \gamma| for a 2d hard rod interactions.


void EqDist::KernCoeffCalcHardRod2D( double* d2nVec, int Nc ){

// Build a vector of the kernal's coefficients.

for( int n = 0; n < Nc; ++n){
    d2nVec[n] = 4 / ( M_PI * ( 4* ( n - 1) * (n - 1) - 1) );
}

}

void EqDist::BestCoeffExpLeg2D( double* Coeff, int Nc, double* phi, double* bc ){

  std::cout << " Not written yet " << std::endl;

}

// trapezoid integration
// Adds another point from final point to initial for periodic fnc
double EqDist::trapz_periodic(double* f, double* x, int N){

  double Intgral = 0.0;
  
  Intgral += f[0];

  for( int i = 1; i < N ; ++i ){

    Intgral += f[i];

  }
  
  Intgral = ( x[N-1] - x[0] ) / (N) * Intgral;

  return Intgral;

}

// trapezoid integration
double EqDist::trapz(double* f, double* x, int N){

  double Intgral = 0.0;
  
  Intgral += f[0];

  for( int i = 1; i < (N - 1); ++i ){

    Intgral += 2 * f[i];

  }

  Intgral += f[N-1];

  Intgral = ( x[N-1] - x[0] ) / (2*N) * Intgral;

  return Intgral;

}

  
double* EqDist::getfisVec(){return fis_;}











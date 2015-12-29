// All things equilbrium density
 #include "EqDist.h"

EqDist::EqDist(){
  
  Nm_ = 0; // Number of gridpoints
  Nc_ = 0; // Number of Coeff
  bc_ = 0.0; // scaled concentration
  
  phi_ = NULL;
  fis_ = NULL;
  feq_ = NULL; // equilbrium distribution
  Coeff_ = NULL;
  d2nVec_ = NULL;

}

EqDist::EqDist(int Nm, double bc){

  DEBUG_ = true;

  Nm_ = Nm; // Number of gridpoints
  Nc_ = 10; // Number of Coeff
  bc_ = bc; // scaled concentration
  
  std::cout << "Initializing vectors" << std::endl;
  spGrid phi;
  phi_ = new double[Nm_];
  phi.xVecMaker(Nm_, 2*M_PI,phi_);

  fis_ = new double[Nm_];  //Isotropic
  feq_ = new double[Nm_];  //Equilibrium
  Coeff_ = new double[Nc_]; // Vector of coeff an, exp[ sum an cos ( 2n phi ) ]
  d2nVec_ = new double[Nc_]; // Kernal

 
  feqInit(); // Equilbrium
  fisInit(); // isotropic 
}

void EqDist::fisInit( ){
 
  for( int i = 0; i < Nm_; ++i ){
    fis_[i] = 1 / (2 * M_PI);
  }

}

void EqDist::feqInit( ){

  // Fix bc if it's too close to 1.5
  bcFix();
 
 // Build Kernal
  KernCoeffCalcHardRod2D();
  
  // Build Coeff
  BestCoeffExpLeg2D();

  // Build Eq
  DistBuilderExpCos2Dsing();

}

void EqDist::fisInit( Ad1 &f ){
 
  for( int i = 0; i < Nm_; ++i ){
    f[i] = fis_[i];
  }

}

void EqDist::feqInit( Ad1 &f ){

  for( int i = 0; i < Nm_; ++i ){
    f[i] = feq_[i];
  }

}

//KernCoeffCalcHardRod2D.m
// Michael Stefferson
// Calculate the coefficients of the Legendre expansion of the kernal
// |sin \gamma| for a 2d hard rod interactions.


void EqDist::KernCoeffCalcHardRod2D( ){

  // Build a vector of the kernal's coefficients.

  for( int n = 0; n < Nc_; n++){
      d2nVec_[n] = 4 / ( M_PI * ( 4* ( n - 1) * (n - 1) - 1) );
  }

  if( DEBUG_){
    std::cout << "Kernal = " << std::endl;
    for( int n = 0; n < Nc_; n++ ){
      std::cout << d2nVec_[n] << std::endl;
    }
  }

}

//Algorithm doesn't converge for bc too close to 1.5
void EqDist::bcFix(){

  if( (bc_ > 1.499) && (bc_ < 1.501) ){
      bc_  = 1.499;
  }

  if( DEBUG_ ) { std::cout << "bc_= " << bc_ << std::endl; }

}

// CoeffCalcExpLeg2D.m
// Michael Stefferson
//
// Calculates the coefficients for a 2D distribution of the form
// f(\theta) = exp( \sum a_2n cos(2*n*\theta) ) / Z
// Inputs:
// Nc: Number of max coefficients
// x : spatial coordin. vector. x = cos(\theta)

void EqDist::BestCoeffExpLeg2D( ){
  
  
  // Make the first guess not zero. Iteration gets stuck
  double ExpCosSum[Nm_];
  double CosExpCos[Nm_];
  double CoeffTemp[Nc_];
  double Coeff_[Nc_];


  // for while loop. How much coefficient need to be changing by
  // before stopping
  double epsilon = 1e-12;
  
  //Integral Stuff
  double NormEc;
  double NormCeC;

  // Make initial guess for the coeff the previous coeff.
  CoeffTemp[0] = 1;
  Coeff_[0] = rand();
  for(int i = 1; i < Nc_ ; i++){

    CoeffTemp[i] = 1;
    Coeff_[i]    = 0;
    
  }

  // Iterate the coefficients
  
  int counter = 0;
  double MaxElement;
  double delta[Nm_];

 
  MaxElement = 10.0;

  while( MaxElement > epsilon ){       

    for( int i = 0; i < Nm_ ; i++ ) {
        CoeffTemp[i] = Coeff_[i];
      }
      
      for( int j = 0; j < Nm_; j++ ){
        ExpCosSum[j] = 1.0 ;
        // Sum over all the Legendre polynomials
        for( int n = 0; n < Nc_; n++ ){
          ExpCosSum[j] *= exp( CoeffTemp[n] * cos( 2.0 * (n + 1.0) * phi_[j] ) );
        }
      }

      for( int n = 0; n < Nc_; n++ ){
        
           for( int j = 0; j < Nm_; j++ ){
              CosExpCos[j] += cos( 2.0 * (n + 1) * phi_[j] ) * ExpCosSum[j];
           }

           NormCeC = trapz_periodicPhi(CosExpCos);
           NormEc  = trapz_periodicPhi(ExpCosSum);
           Coeff_[n] = M_PI * bc_ * d2nVec_[n] * NormCeC / NormEc;

      }
      
      counter += 1;

    for( int i = 0; i < Nm_; i++ ){
      delta[i] = Coeff_[i] - CoeffTemp[i];
    }
 
    MaxElement = MaxMagElement(delta,Nm_);    
    } // end while
}

void EqDist::DistBuilderExpCos2Dsing(){
  
  for( int n = 0; n < Nc_; n++ ){
    for( int j = 0; j < Nm_; j++ ){
      
      feq_[j]  *=  exp( Coeff_[n] * cos(2.0 * (n+1) * phi_[j] ) );    

    }
  }

  // Normalize it again to be safe.
  Normalize(feq_);

}

// trapezoid integration
// Adds another point from final point to initial for periodic fnc

void EqDist::Normalize( double f[] ){
  
  double Norm;
  Norm =  trapz_periodicPhi(f);

  for(int i = 0; i < Nm_; ++i ) {
    
    f[i] = f[i] / Norm;

  }
}


double EqDist::trapz_periodicPhi(const double f[]){

  double Intgral = 0.0;
  
  for( int i = 0; i < Nm_ ; ++i ){

    Intgral += 2 * f[i];

  }
  
  Intgral =  2 * M_PI / Nm_ * Intgral;

  return Intgral;

}
  
double EqDist::MaxMagElement(double v[], int size){
  
  double Max = 0.0;
  double Temp;
  
  for( int i = 0; i < size; i++ ){
    
    Temp = v[i] * v[i];

    if( Temp > Max * Max ){
      Max = Temp;
    }

  }

  return Max;

}
      
void EqDist::printFis(){

  std::cout << "fis= " << std::endl;

  for( int i = 0; i < Nm_; i++ ){
    std::cout << fis_[i] << std::endl;
  }

  std::cout << std::endl;

}

void EqDist::printFeq(){

  std::cout << "feq= " << std::endl;

  for( int i = 0; i < Nm_; i++ ){
    std::cout << feq_[i] << std::endl;
  }

  std::cout << std::endl;

}

double* EqDist::getfisVec(){return fis_;}











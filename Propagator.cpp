// Propagator.cpp

#include "Propagator.h"

//Constructor

Propagator::Propagator()
{
  IsoFlag_ = 0;
  N1_      = 0;
  N2_      = 0;
  N3_      = 0;
  LopFTd_   = NULL;
  Ud_       = NULL;
} 

// 1D
Propagator::Propagator( int N1, double* k1, double D, double dt, bool Iso)
{
  IsoFlag_ = Iso;
  N1_ = N1;
  N2_ = 0;
  N3_ = 0;
  if( IsoFlag_ == 1)
  {
  
    LopFTd_ = new double [N1_];
    Ud_     = new double [N1_];

    LopDiagMaker(k1, D);
    PropIsoMaker1(dt);
  
  }
  else //not written
  {
  
    LopFTd_ = NULL;
    Ud_   = NULL;
  }

}

// 3D
Propagator::Propagator( int N1, int N2, int N3, double* k1, double* k2, double* k3,
                         double D, double dt, bool Iso)
{
  IsoFlag_ = Iso;
  N1_ = N1;
  N2_ = N2;
  N3_ = N3;
  if( IsoFlag_ == 1)
  {
  
    LopFTd_ = new double [ N1_ * N2_ * N3_];

    Ud_ = new double [ N1_ * N2_ * N3_];

    LopDiagMaker(k1, k2, k3, D);
    PropIsoMaker3(dt);
  
  }
  else //not written
  {
  
    LopFTd_ = NULL;
    Ud_   = NULL;
  }

}


void Propagator::LopDiagMaker(double* kt, double D)
{
  for( int i = 0; i < N1_; ++i)
  {
   LopFTd_[i] = -D * kt[i] * kt[i];
  }
 
}


void Propagator::LopDiagMaker(double* k1t, double* k2t, double* k3t, double D)
{
  for( int i = 0; i < N1_; ++i){
    for( int j = 0; j < N2_; ++ j){
      for ( int k = 0; k < N3_; ++k ){
     LopFTd_[i + j*N1_ +  k*N1_*N2_] = -D * ( 
                                      k1t[i] * k1t[i] + k2t[j] * k2t[j] + k3t[k] * k3t[k] );
      }
    }
  }
}


 void Propagator::PropIsoMaker1(double dt)
{
 for( int i = 0; i < N1_; ++i ){
   Ud_[i] = exp( LopFTd_[i] * dt);
 }
}

void Propagator::PropIsoMaker3(double dt)
{
 for( int i = 0; i < N1_; ++i ){
   for( int j = 0; j < N2_; ++j ){
     for( int k = 0; k < N3_; ++k ){
        Ud_[i + j*N1_ +  k*N1_*N2_] = exp( LopFTd_[i + j*N1_ +  k*N1_*N2_] * dt);
     }
   }
 }
}

void  Propagator::PropPrint(){
  std::cout << "Prop as cube:";
 for( int k = 0; k < N3_; ++k ){
   std::cout << "km = " << k << std::endl;
   for( int i = 0; i < N1_; ++i ){
     for( int j = 0; j < N2_; ++j ){
       std::cout << Ud_[i + j*N1_ +  k*N1_*N2_] << "\t";
     }
     std::cout << std::endl;
   }
   std::cout << std::endl << std::endl;
 }// end for loop
}

double* Propagator::getLop(){ return LopFTd_; }
double* Propagator::getProp(){ return Ud_; }



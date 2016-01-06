// Propagator.cpp

#include "Propagator.h"

//Constructor

Propagator::Propagator()
{
  IsoFlag_ = 0;
  N1_      = 0;
  N2_      = 0;
  N3_      = 0;
  dt_      = 0;
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
  dt_ = dt;
  if( IsoFlag_ == 1)
  {
  
    LopFTd_ = new double [N1_];
    Ud_     = new double [N1_];

    LopIsoDiagMaker(k1, D);
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
                         double D, double dt, bool Iso){

  IsoFlag_ = Iso;
  N1_ = N1;
  N2_ = N2;
  N3_ = N3;
  dt_ = dt;
  if( IsoFlag_ == 1)
  {
  
    LopFTd_ = new double [ N1_ * N2_ * N3_];

    Ud_ = new double [ N1_ * N2_ * N3_];

    LopIsoDiagMaker(k1, k2, k3, D);
    PropIsoMaker3(dt);
  
  }
  else //not written
  {
  
    LopFTd_ = NULL;
    Ud_   = NULL;
  }

}

Propagator::Propagator( int N1, int N2, int N3, double* k1, double* k2, double* k3,
                         double Dpos, double Dr, double dt, bool Iso){

  IsoFlag_ = Iso;
  N1_ = N1;
  N2_ = N2;
  N3_ = N3;
  dt_ = dt;
  if( IsoFlag_ == 1)
  {
  
    LopFTd_ = new double [ N1_ * N2_ * N3_];

    Ud_ = new double [ N1_ * N2_ * N3_];

    LopIsoDiagMaker(k1, k2, k3, Dpos, Dr);
    PropIsoMaker3(dt);
    PreFacInit();
  
  }
  else //not written
  {
  
    LopFTd_ = NULL;
    Ud_   = NULL;
  }

}


void Propagator::LopIsoDiagMaker(double* kt, double D)
{
  for( int i = 0; i < N1_; ++i)
  {
   LopFTd_[i] = -D * kt[i] * kt[i];
  }
 
}


void Propagator::LopIsoDiagMaker(double* k1t, double* k2t, double* k3t, double D)
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

void Propagator::LopIsoDiagMaker(double* k1t, double* k2t, double* k3t, double Dpos, double Dr)
{
  for( int i = 0; i < N1_; ++i){
    for( int j = 0; j < N2_; ++ j){
      for ( int k = 0; k < N3_; ++k ){
     LopFTd_[i + j*N1_ +  k*N1_*N2_] = - ( Dpos * ( k1t[i] * k1t[i] + k2t[j] * k2t[j] ) 
                                      + Dr * k3t[k] * k3t[k] );
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

void Propagator::PreFacInit(){

  switch(StepFlag_) {
    case 0: // Diffusion
      NlPf_ = 0;
      NlPfPrev_ = 0;
      NlPfExp_ = 0;
      break;
    case 1: // AB1
      NlPf_ = dt_;
      NlPfPrev_ = 0;
      NlPfExp_ = 0;
      break;
    case 2: // AB2
      NlPf_ = 3.0 * dt_ / 2.0;
      NlPfPrev_ = dt_ / 2.0;
      NlPfExp_ = 0;
      break;
    case 3: // HAB 1
      NlPf_ = dt_;
      NlPfPrev_ = 0;
      NlPfExp_ = 0;
      break;
    case 4: // HAB 2
      NlPf_ = 3.0 * dt_ / 2.0;
      NlPfPrev_ = dt_ / 2.0;
      NlPfExp_ = 0;
      break;
    case 5:  // BHAB 1
      NlPf_ = dt_ / 2.0;
      NlPfPrev_ = 0;
      NlPfExp_ = dt_ / 2.0;
      break;
    case 6: // BHAB 2
      NlPf_ = dt_;
      NlPfPrev_ = dt_ / 2.0;
      NlPfExp_ = dt_ / 2.0;
      std::cout << "Error: need NlFT previous" << std::endl;
      break;
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


//////////////////////// Propagator Functions ///////////////////////


void Propagator::PropMaster( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT) {
 
if( IsoFlag_ == 1 ){
  
  switch(StepFlag_){
    case 0: 
      PropAB0c( rhoFTnext, rhoFT); // Diffusion
      break;
    case 1: 
      PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1
      break;
    case 2: 
      std::cout << "Error: need NlFT previous" << std::endl;
      break;
    case 3:
      PropHAB1c( rhoFTnext, rhoFT, NlFT);
      break;
    case 4:
      std::cout << "Error: need NlFT previous" << std::endl;
      break;
    case 5: 
      PropBHAB1c( rhoFTnext, rhoFT, NlFT);
      break;
    case 6:
      std::cout << "Error: need NlFT previous" << std::endl;
      break;
  } // switch
}
else{ std::cout << "not written yet" << std::endl; }
}

void Propagator::PropMaster( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev) {
 
if( IsoFlag_ == 1 ){
  switch(StepFlag_){

    case 0: // Diffusion
      PropAB0c( rhoFTnext, rhoFT); 
      break;
    case 1: // AB 1
      PropAB1c( rhoFTnext, rhoFT, NlFT ); 
      break;
    case 2: // AB 2
      PropAB2c( rhoFTnext, rhoFT, NlFT, NlFTprev); 
      break;
    case 3: // HAB1
      PropHAB1c( rhoFTnext, rhoFT, NlFT);
      break;
    case 4: // HAB2
      PropHAB2c( rhoFTnext, rhoFT, NlFT, NlFTprev);
      break;
    case 5: // BHAB1
      PropBHAB1c( rhoFTnext, rhoFT, NlFT);
      break;
    case 6: // BHAB2
      PropBHAB2c( rhoFTnext, rhoFT, NlFT, NlFTprev);
      break;
  } // switch
}
else{ std::cout << "not written yet" << std::endl; }
}

void Propagator::PropAB0c( Ac3 &rhoFTnext, Ac3 &rhoFT ) {

  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] = Ud_[i + j*N1_ +  k * N1_ * N2_] * rhoFT[i][j][k];

      }
    }
  }
} 

void Propagator::PropAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT ) {

  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] = Ud_[i + j*N1_ +  k * N1_ * N2_] 
          * ( rhoFT[i][j][k] ) + NlPf_ * NlFT[i][j][k] ;
        
      }
    }
  }
}

void Propagator::PropAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev ) {
  
  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] = Ud_[i + j*N1_ +  k * N1_ * N2_] * ( rhoFT[i][j][k] ) 
          + NlPf_ *  NlFT[i][j][k] - NlPfPrev_ * NlFTprev[i][j][k]    ;
              
      }
    }
  }
}

// Currently, propagator not included in prefactor
void Propagator::PropHAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT ) {

  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] = Ud_[i + j*N1_ +  k * N1_ * N2_] 
          * ( rhoFT[i][j][k]  + NlPf_ * NlFT[i][j][k] );
        
      }
    }
  }
}

// Currently, propagator not included in prefactor
void Propagator::PropHAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev ) {
  
  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] = Ud_[i + j*N1_ +  k * N1_ * N2_] * 
          ( rhoFT[i][j][k]  + NlPf_ *  NlFT[i][j][k] 
            - Ud_[i + j*N1_ +  k * N1_ * N2_] *  NlPfPrev_ * NlFTprev[i][j][k] )   ;
              
      }
    }
  }
}

// Currently, propagator not included in prefactor
void Propagator::PropBHAB1c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT ) {

  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] = Ud_[i + j*N1_ +  k * N1_ * N2_] 
          * ( rhoFT[i][j][k]  + NlPf_ * NlFT[i][j][k] ) + NlPf_ * NlFT[i][j][k];
        
      }
    }
  }
}

// Currently, propagator not included in prefactor
void Propagator::PropBHAB2c( Ac3 &rhoFTnext, Ac3 &rhoFT, Ac3 &NlFT, Ac3 &NlFTprev ) {
  
  for(int i = 0; i <  N1_ ; ++i){
    for( int j = 0; j <  N2_ ; ++j){
      for( int k = 0; k < N3_; ++k){
    
        rhoFTnext[i][j][k] =  Ud_[i + j*N1_ +  k * N1_ * N2_] * 
          ( rhoFT[i][j][k] +  NlPfExp_ *  NlFT[i][j][k] )
          + NlPf_ * NlFT[i][j][k] - NlPfPrev_ *  NlFTprev[i][j][k];
              
      }
    }
  }
}

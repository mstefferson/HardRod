// tGrid.cpp

#include "tGrid.h"

// Constructor
tGrid::tGrid()
{
  dt_ = 0;
  tend_ = 0;
  trec_ = 0;
  Nt_    = 0;  
  Nrec_ = 0;
  Ncount_ = 0;    
  tV_ = NULL;  
  trecV_ = NULL;
 }

tGrid::tGrid(double dt, double trec, double tend)
{
  dt_     = dt;
  
  if( trec < dt ){ trec = dt;}

  Ncount_ = floor(trec/dt);   
  trec_   = Ncount_ * dt_;
  
  if( tend < trec ){ tend = trec; }

  Nrec_   = floor(tend/trec_);   
  tend_   = Nrec_ * trec_;
  
  Nt_     = Nrec_ * Ncount_;  

  tV_     = new double [Nt_+1];  
  trecV_  = new double [Nrec_+1];
  
  tVecMaker(Nt_, dt_, tV_);
  tVecMaker(Nrec_, trec_, trecV_);
 }
    
// function that prints everything
void tGrid::print()
{
  std::cout << "dt = " << dt_ << std::endl;
  std::cout << "trec = " << trec_ << std::endl;
  std::cout << "tend = " << tend_ << std::endl;
  std::cout << "Nt  = " << Nt_ << std::endl;
  std::cout << "Nrec = " << Nrec_ << std::endl;
  std::cout << "Ncount = " << Ncount_ << std::endl;

  std::cout << "tvec = " << std::endl;
  for( int i = 0; i <= Nt_ ; ++i )
  { 
    std::cout << tV_[i] << std::endl; 
  }
  std::cout << std::endl;

  std::cout << "trec = " << std::endl;
  for( int i = 0; i <= Nrec_; ++i )
  { 
    std::cout << trecV_[i] << std::endl; 
  }
  std::cout << std::endl;


}

// Builds time vector
void  tGrid::tVecMaker( int N, double timespace, double* tvec)
{
 // tend = N * timespace
  for( int i = 0; i < N+1; i++)
  {
tvec[i] = i * timespace;
  }

}
 
double  tGrid::getDt(){ return dt_; }
double  tGrid::getTend(){ return tend_; }
double  tGrid::getTrec(){ return trec_; }
int  tGrid::getNt(){ return Nt_; }
int  tGrid::getNrec(){ return Nrec_; }

int  tGrid::getNcount(){ return Ncount_; }
double*  tGrid::getTvec(){ return tV_; }
double* tGrid::getTrecVec(){ return trecV_; }

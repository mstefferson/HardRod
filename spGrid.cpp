// spGrid.cpp

#include "spGrid.h"

// Constructor
spGrid::spGrid() 
{
  N_      = 0;
  L_      = 0;
  dx_     = 0;
  pos_ = NULL;
  k_ = NULL;
  kpos_ = NULL;
}  

spGrid::spGrid(int N, double L) 
{
  N_ = N;
  L_ = L;
  dx_ = dxMaker(N_,L_);
  pos_  = new double[N_];
  k_    = new double [N_];
  kpos_ = new double [N_/2 + 1];
  kFT_  = new double [N_];
  xVecMaker();
  kVecMaker();
  kposVecMaker();
  kFtVecMaker();
}
void spGrid::print()
{  
  std::cout << "N = " << N_  << std::endl;
   
  std::cout << "L = " << L_  << std::endl;
   
  std::cout << "dx = " << dx_ << std::endl;
  std::cout << std::endl;

  std::cout << "x vec = " << std::endl;
  for( int i = 0; i < N_; i++ )
  { 
    std::cout << pos_[i] << std::endl; 
  }
  std::cout << std::endl;

  std::cout << "k vec = " << std::endl;
  for( int i = 0; i < N_; i++ )
  { 
    std::cout << k_[i] << std::endl; 
  }
  std::cout << std::endl;

  std::cout << "kpos_ vec = " << std::endl;
  for( int i = 0; i < N_/2+1; i++ )
  { 
    std::cout << kpos_[i] << std::endl; 
  }
  std::cout << std::endl;
  
  std::cout << "kFT_ vec = " << std::endl;
  for( int i = 0; i < N_; i++ )
  { 
    std::cout << kFT_[i] << std::endl; 
  }
  std::cout << std::endl;
}

// Builds position vecotr
double spGrid::dxMaker() 
{
  return L_ / N_;
}

double spGrid::dxMaker(int N,double L) 
{
  return L / N;
}


// Builds poisition vector
void spGrid::xVecMaker()
{
  pos_[0] = -L_/2;

  for(int i = 1; i < N_; i++)
  {
    pos_[i] = pos_[i-1] + dx_;
  }
}
void spGrid::xVecMaker(int N, double L, double* pos )
{
  double dx = L/N;

  pos[0] = -L/2;

  for(int i = 1; i < N; i++)
  {
    pos[i] = pos[i-1] + dx;
  }
}

// Builds k vecotr
void spGrid::kVecMaker()
{
  double dk = 2 * M_PI  / L_;
  k_[0] = -M_PI * N_ / L_;
  
  for( int i = 1; i < N_; i ++ )
  {
     k_[i] = k_[i-1] + dk;
  }

}

void spGrid::kVecMaker(int N, double L, double *k)
{
  double dk = 2 * M_PI  / L;
  k[0] = -M_PI * N / L;
  
  for( int i = 1; i < N; ++i )
  {
     k[i] = k[i-1] + dk;
  }

}


// Builds k >= 0 vector
void spGrid::kposVecMaker()
{
  for(int i = 0; i < N_/2; ++i)
  {
    kpos_[i] = k_[N_/2 +  i];
   }
  kpos_[N_/2] = -k_[0];

}

// Builds k >= 0 vector
void spGrid::kposVecMaker(int N, double *k, double *kpos)
{
  for(int i = 0; i < N/2; ++i)
  {
    kpos[i] = k[N/2 +  i];
    std::cout << k[i] << std::endl;
   }
  kpos[N/2] = -k[0];

}

void spGrid::kFtVecMaker()
{
   double dk = 2 * M_PI  / L_;
   kFT_[0] = 0;
   kFT_[N_/2] = M_PI * N_ / L_;

  for( int i = 1; i < N_/2; ++i )
  {
     kFT_[i] = i*dk;
     kFT_[N_-i] =  -i*dk;
  }
}
  

void spGrid::kFtVecMaker(int N, double L, double *k)
{
   double dk = 2 * M_PI  / L;
   k[0] = 0;
   k[N/2] = M_PI * N / L;

  for( int i = 1; i < N/2; ++i )
  {
     k[i] = i*dk;
     k[N-i] =  -i*dk;
  }
}

     // return position vecto
 double* spGrid::getPos(){ return pos_; }
    
     // return k vec
 double* spGrid::getK(){ return k_; }
    
    // return positive k values
 double* spGrid::getKpos(){ return kpos_; }
 double* spGrid::getKft(){ return kFT_; }
    

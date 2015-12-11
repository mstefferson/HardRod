// OPclass.cpp

#include "OPclass.h"

// Constructor
OPclass::OPclass( int Nx, int Ny, int Nm, double* phi ){
  
  Nx_ = Nx;
  Ny_ = Ny;
  Nm_ = Nm;
  
  intFac_ = 2 * M_PI / Nm_;

  nxTemp_ = 0;
  nyTemp_ = 0;
  QxyTemp_ = 0;
  QyyTemp_ = 0;
  QxxTemp_ = 0;

  sin_ = new double[Nm_];
  cos_ = new double[Nm_];
  sinsin_ = new double[Nm_];
  coscos_ = new double[Nm_];
  cossin_ = new double[Nm_];

  C_ = new double*[Nx_];
  PO_ = new double*[Nx_];
  NO_ = new double*[Nx_];

  for( int i = 0; i < Nx_; ++i){
    C_[i] =  new double[Ny_];
    PO_[i] =  new double[Ny_];
    NO_[i] =  new double[Ny_];
  }


  for( int i = 0; i < Nm_; ++i ){
    sin_[i] = sin( phi[i] );
    cos_[i] = cos( phi[i] );
    sinsin_[i] = sin( phi[i] ) * sin( phi[i] );
    cossin_[i] = cos( phi[i] ) * sin( phi[i] );
    coscos_[i] = cos( phi[i] ) * cos( phi[i] );
  }

}


// Make OP

void OPclass::OPmaker(Ad3& rho){

  ConcCalc(rho);
  PolarOrdCalc(rho);
  NemOrdCalc(rho);

}

// Concentration
void OPclass::ConcCalc(Ad3& rho){

    for( int i = 0; i < Nx_; ++i ){
      for( int j = 0; j < Ny_; ++j){
        C_[i][j] = 0;
        for( int k = 0; k < Nm_; ++k){

          C_[i][j] += rho[i][j][k];

        }
        C_[i][j] *= intFac_;
      }
    }
}

// Polar Order
void OPclass::PolarOrdCalc(Ad3& rho){

    for( int i = 0; i < Nx_; ++i ){
      for( int j = 0; j < Ny_; ++j ){
        PO_[i][j] = 0;
        nxTemp_ = 0;
        nyTemp_ = 0;
        for( int k = 0; k < Nm_; ++k){

          nxTemp_ += rho[i][j][k] * cos_[k];
          nyTemp_ += rho[i][j][k] * sin_[k];

        }
        nxTemp_ *= intFac_;
        nyTemp_ *= intFac_;
        PO_[i][j] = sqrt( nxTemp_ * nxTemp_ + nyTemp_ * nyTemp_) / C_[i][j];
      }
    }
  }

// Nematic Order
void OPclass::NemOrdCalc(Ad3& rho){

    for( int i = 0; i < Nx_; ++i ){
      for( int j = 0; j < Ny_; ++j ){
        NO_[i][j] = 0;
        QxxTemp_ = 0;
        QxyTemp_ = 0;
        QyyTemp_ = 0;

        for( int k = 0; k < Nm_; ++k ){

          QxxTemp_ += rho[i][j][k] * ( coscos_[k] - 1/2 );
          QxyTemp_ += rho[i][j][k] * cossin_[k];
          QyyTemp_ += rho[i][j][k] * ( sinsin_[k] - 1/2 );

        }

        QxxTemp_ *= intFac_;
        QxyTemp_ *= intFac_;
        QyyTemp_ *= intFac_;
       
        // Nem order equal to max eigvector of Nem Mtrx * 2 
        NO_[i][j] = ( ( QxxTemp_ + QyyTemp_ )  + 
        sqrt( ( ( QxxTemp_ - QyyTemp_ ) * ( QxxTemp_ - QyyTemp_ ) )  + 
            4 * QxyTemp_ * QxyTemp_ ) ) / C_[i][j];
        if( i == Nx_ - 1 && j == Ny_ - 1 ){
          std::cout << "Qxx = " << QxxTemp_ << std::endl;
          std::cout << "Qxy = " << QxyTemp_ << std::endl;
          std::cout << "Qyy = " << QyyTemp_ << std::endl;
          std::cout << "Nij = " << NO_[i][j] << std::endl;
        }
      }
    }
  }

void OPclass::printC(){
  for( int i = 0; i < Nx_; i++ ){
    for( int j = 0; j < Ny_; j++ ) {
      std::cout << C_[i][j] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void OPclass::printPO(){
  for( int i = 0; i < Nx_; i++ ){
    for( int j = 0; j < Ny_; j++ ) {
      std::cout << PO_[i][j] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void OPclass::printNO(){
  for( int i = 0; i < Nx_; i++ ){
    for( int j = 0; j < Ny_; j++ ) {
      std::cout << NO_[i][j] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

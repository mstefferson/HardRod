// HRwriter.cpp

#include "HRwriter.h"

HRwriter::HRwriter(int Nx, int Ny, int trial, double** C, double** PO, double** NO){

  Nx_ = Nx;
  Ny_ = Ny;
  trial_ =  trial;

  cTl_ << "Conc" << trial << ".txt";
  pTl_ << "PO" << trial << ".txt";
  nTl_ << "NO" << trial << ".txt";
  rhoTl_ << "Rho" << trial << ".txt";
  ampTl_ << "Amp" << trial << ".txt";
  paramsTl_ << "Params" << trial << ".txt";
  //cTl_ = "Conc" + std::to_string(trial) + ".txt";

  C_ = C;
  PO_ = PO;
  NO_ = NO;

  openFiles();
}

void HRwriter::openFiles(){

  concFile_.open( ( cTl_.str() ).c_str() );
  poFile_.open( ( pTl_.str() ).c_str() );
  noFile_.open( ( nTl_.str() ).c_str() );
  rhoFile_.open( ( rhoTl_.str() ).c_str() );
  ampFile_.open( ( ampTl_.str() ).c_str() );
  paramsFile_.open( ( paramsTl_.str() ).c_str() );

}

void HRwriter::closeFiles(){

  concFile_.close();
  poFile_.close();
  noFile_.close();
  rhoFile_.close();
  ampFile_.close();
  paramsFile_.close();

}

void HRwriter::writeOP(){

  writeC();
  writeNO();
  writePO();

}

void HRwriter::writeC(){

  for( int i = 0; i < Nx_; i++){
    for( int j = 0; j < Ny_; j++){
      concFile_ << C_[i][j] << "\t";
    }
    concFile_ << std::endl;
  }
  concFile_ << std::endl;
}

void HRwriter::writePO(){

  for( int i = 0; i < Nx_; i++){
    for( int j = 0; j < Ny_; j++){
      poFile_ << PO_[i][j] << "\t";
    }
    poFile_ << std::endl;
  }
  poFile_ << std::endl;
}

void HRwriter::writeNO(){

  for( int i = 0; i < Nx_; i++){
    for( int j = 0; j < Ny_; j++){
      noFile_ << NO_[i][j] << "\t";
    }
    noFile_ << std::endl;
  }
  noFile_ << std::endl;
}

void HRwriter::writeAmp(Ac3 &rhoFT){

  ampFile_ << rhoFT[0][0][2].real() << "\t";
  ampFile_ << rhoFT[0][0][2].imag() << "\t";
/*
 *  ampFile_ << rhoFT[1][0][2].real << "\t";
 *  ampFile_ << rhoFT[1][0][2].imag << "\t";
 *  ampFile_ << rhoFT[0][1][2].real << "\t";
 *  ampFile_ << rhoFT[0][1][2].imag << "\t";
 *  ampFile_ << rhoFT[1][1][2].real << "\t";
 *  ampFile_ << rhoFT[1][1][2].imag << "\t";
 *
 */
  ampFile_ << std::endl;

}


void HRwriter::writeRho(Ad3& rho){

  rhoFile_ << rho << std::endl;

}

void HRwriter::writeParams(system_params params){

  paramsFile_ << params.trial << std::endl;
  paramsFile_ << params.Nx << std::endl;
  paramsFile_ << params.Ny << std::endl;
  paramsFile_ << params.Nm << std::endl;
  paramsFile_ << params.pXmodes << std::endl;
  paramsFile_ << params.pYmodes << std::endl;
  paramsFile_ << params.pMmodes << std::endl;
  paramsFile_ << params.dt << std::endl;
  paramsFile_ << params.trec << std::endl;
  paramsFile_ << params.tend << std::endl;
  paramsFile_ << params.bc << std::endl;
  paramsFile_ << params.vD << std::endl;
  paramsFile_ << params.Lx << std::endl;
  paramsFile_ << params.Ly << std::endl;
  paramsFile_ << params.Dpar << std::endl;
  paramsFile_ << params.Dperp << std::endl;
  paramsFile_ << params.Dr << std::endl;
  paramsFile_ << params.PertrbAmp << std::endl;

}

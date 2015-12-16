// HRwriter.cpp

#include "HRwriter.h"

HRwriter::HRwriter(int Nx, int Ny, int Nm, int trial, double** C, double** PO, double** NO, double* rhoFix, 
     Complex* rhoFtFix1, Complex* rhoFtFix2, Complex* rhoFtFix3, Complex* rhoFtFix4,
     Complex* rhoFtFix5, Complex* rhoFtFix6, Complex* rhoFtFix7, Complex* rhoFtFix8){

  Nx_ = Nx;
  Ny_ = Ny;
  Nm_ = Nm;
  trial_ =  trial;
  recNum_ = 0;

  cTl_ << "Conc" << trial << ".txt";
  pTl_ << "PO" << trial << ".txt";
  nTl_ << "NO" << trial << ".txt";
  rhoTl_ << "Rho" << trial << ".txt";
  ampTl_ << "Amp" << trial << ".txt";
  paramsTl_ << "Params" << trial << ".txt";
  distTl_ << "Dist" << trial << ".txt";
  //cTl_ = "Conc" + std::to_string(trial) + ".txt";

  C_ = C;
  PO_ = PO;
  NO_ = NO;
  rhoFixP_   = rhoFix;
  rhoFtFix1_ = rhoFtFix1;
  rhoFtFix2_ = rhoFtFix2;
  rhoFtFix3_ = rhoFtFix3;
  rhoFtFix4_ = rhoFtFix4;
  rhoFtFix5_ = rhoFtFix5;
  rhoFtFix6_ = rhoFtFix6;
  rhoFtFix7_ = rhoFtFix7;
  rhoFtFix8_ = rhoFtFix8;

  openFiles();
}

void HRwriter::openFiles(){

  concFile_.open( ( cTl_.str() ).c_str() );
  poFile_.open( ( pTl_.str() ).c_str() );
  noFile_.open( ( nTl_.str() ).c_str() );
  rhoFile_.open( ( rhoTl_.str() ).c_str() );
  ampFile_.open( ( ampTl_.str() ).c_str() );
  paramsFile_.open( ( paramsTl_.str() ).c_str() );
  distFile_.open( ( distTl_.str() ).c_str() );

}

void HRwriter::closeFiles(){

  paramsFile_ << recNum_;
  
  concFile_.close();
  poFile_.close();
  noFile_.close();
  rhoFile_.close();
  ampFile_.close();
  distFile_.close();
  paramsFile_.close();

}

void HRwriter::writeOP(){

  writeC();
  writeNO();
  writePO();
  writeDist();
  writeAmp();
  recNum_ ++;
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

void HRwriter::writeDist(){

  for( int i = 0; i < Nm_ ; i++ ){
    distFile_ << rhoFixP_[i] << std::endl;
  }
  distFile_ << std::endl;

}


void HRwriter::writeAmp( ){


  ampFile_ << (*rhoFtFix1_).real() << "\t";
  ampFile_ << (*rhoFtFix1_).imag() << "\t";

  ampFile_ << (*rhoFtFix2_).real() << "\t";
  ampFile_ << (*rhoFtFix2_).imag() << "\t";

  ampFile_ << (*rhoFtFix3_).real() << "\t";
  ampFile_ << (*rhoFtFix3_).imag() << "\t";

  ampFile_ << (*rhoFtFix4_).real() << "\t";
  ampFile_ << (*rhoFtFix4_).imag() << "\t";

  ampFile_ << (*rhoFtFix5_).real() << "\t";
  ampFile_ << (*rhoFtFix5_).imag() << "\t";

  ampFile_ << (*rhoFtFix6_).real() << "\t";
  ampFile_ << (*rhoFtFix6_).imag() << "\t";

  ampFile_ << (*rhoFtFix7_).real() << "\t";
  ampFile_ << (*rhoFtFix7_).imag() << "\t";

  ampFile_ << (*rhoFtFix8_).real() << "\t";
  ampFile_ << (*rhoFtFix8_).imag() << "\t";

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

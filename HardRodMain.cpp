// HardRodMain.cpp
// Michael Stefferson
//
// Description: Test program of diffusion in 1D using fft
// Compile:

#include "HardRodMain.h"


int main(){
    fftwpp::fftw::maxthreads=get_max_threads(); //Multithreads for fft
  
  //Size of vectors to fft
  unsigned int Nkx,Nky,Nkm;
  size_t align=sizeof(Complex);
 
  std::cout << "hello" << std::endl;
  system_params params;

  params.init();
  
  std::string ParamsFile;
  ParamsFile = "Params.yaml";
  parse_params(ParamsFile, params);
 
  std::cout << "Nx after yaml test " << params.Nx << std::endl;
  std::cout << " trial after yaml test " << params.trial << std::endl;

  //Build density matrix
  Nkx = params.Nx/2 + 1; Nky = params.Ny/2 +1; Nkm = params.Nm / 2 + 1;

  Array::array3<double> c(params.Nx,params.Ny,params.Nm,align);
  Array::array3<Complex> cFT(params.Nx,params.Ny,Nkm,align);
  Array::array3<Complex> cFTnext(params.Nx,params.Ny,Nkm,align);
 
  Array::array3<int> MayerFnc( params.Nx, params.Ny, params.Nm,align);
  Array::array3<Complex> MayerFncFT( params.Nx, params.Ny,Nkm,align);
   
  MayerFncHardRod(params.Nx,params.Ny,params.Nm,params.Lx,params.Ly,params.Lrod,MayerFnc);

 // Print it
 double dphi = 2 * M_PI / params.Nm;
  for( int k  = 0; k < params.Nm; ++k ){
    std::cout << " k = " << k << " phi =  " << k* dphi << std::endl;
    for( int i = 0; i < params.Nx; ++i ) {
      for( int j = 0; j < params.Ny; ++j ) {
    std::cout << MayerFnc[i][j][k] << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
  }

  fftwpp::rcfft3d Forward3(params.Nm,c,cFT);
  fftwpp::crfft3d Backward3(params.Nm,cFT,c);

  // Mayer Fnc

  // Time
  
  tGrid time(params.dt,params.trec,params.tend);
  
  // Build grid 

  spGrid Gridx(params.Nx,params.Lx);
  spGrid Gridy(params.Ny,params.Ly);
  spGrid Gridphi(params.Nm, 2 * M_PI);
  
  double*  kx =  Gridx.getKft();
  //double* kx = Gridx.getKpos();
  double*  x = Gridx.getPos();
  //std::cout << "x" << std::endl;
  //Gridx.print();
  double*  ky =  Gridy.getKft();
  //double* ky = Gridy.getKpos();
  double*  y = Gridy.getPos();
  //std::cout << "y" << std::endl;
  //Gridy.print(); 
  double*  km =  Gridphi.getKpos();
  double*  phi = Gridphi.getPos();
  //std::cout << "phi" << std::endl;
  //Gridphi.print();

  //Propagator DiffP(Nkx,Nky,Nkm, kx,ky,km, D, dt, Iso);
  Propagator DiffP(params.Nx,params.Ny,Nkm, kx,ky,km, params.Dr, params.dt, params.IsoDiffFlag);
  double* Lop = DiffP.getLop();
  double* Prop = DiffP.getProp();

  //Ml3dPrint( Lop, Nx, Ny, Nkm );
  //Ml3dPrint( Prop, Nx, Ny, Nkm );
 
for( int ik = 0; ik < params.Nx; ++ik ){
  for( int jk = 0; jk <  params.Ny; ++jk ){
    for( int kk = 0; kk < Nkm; ++kk ){
       
      if( ik <= params.pXmodes && jk <= params.pYmodes && kk <= params.pMmodes ){
        cFTnext[ik][jk][kk] = params.PertrbAmp * ( params.Nx * params.Ny * params.Nm);
      }
      else{
        cFTnext[ik][jk][kk] = 0;
      }
    }
  }
}
  cFTnext[0][0][0] = 2.0 * (params.Nx*params.Ny*params.Nm);

  cFT = cFTnext;
  Backward3.fftNormalized(cFTnext,c);

   
  // Write to file
  
  std::ofstream diffFile;
  diffFile.open ("DiffOut.txt");
 

//  diffFile << cFT << std::endl;
  //std::cout << cFT << std::endl; 

// First Step

  for(int i = 0; i <  params.Nx; ++i){
    for( int j = 0; j <  params.Ny; ++j){
      for( int k = 0; k < Nkm; ++k){
    cFTnext[i][j][k] = Prop[i + j*Nkx +  k*Nkx*Nky] * cFT[i][j][k];
      }
    }
  }
  

//std::cout << "t = 0 " << std::endl << cFT << std::endl;
std::cout << "t = 0" << std::endl;
diffFile << c << std::endl;
//std::cout << "t = 0 " << std::endl << c << std::endl;
// Main loop

for(int t = 1; t <= time.getNt(); ++t)
{
  cFT = cFTnext;
  Backward3.fftNormalized(cFTnext,c);   
  
  for(int i = 0; i <  params.Nx; ++i){
    for( int j = 0; j <  params.Ny; ++j){
      for( int k = 0; k < Nkm; ++k){
    cFTnext[i][j][k] = Prop[i + j*Nkx +  k*Nkx*Nky] * cFT[i][j][k];
      }
    }
  }
  
  
  if( t % time.getNcount() == 0 )
  {
     std::cout << "t = " << t << std::endl; 
      // diffFile << cFT << std::endl;
     diffFile << c << std::endl;
    //std::cout << cFT << std::endl;
    //std::cout << c << std::endl;
  } // recording
} // time loop 
  
 
diffFile.close();

  return 0;
  }

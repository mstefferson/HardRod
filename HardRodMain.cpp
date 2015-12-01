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
  size_t align=sizeof(Complex); //Complex defined in fftw.h std::complex<double>

  //Parse parameters
  system_params params;
  params.init();
  std::string ParamsFile;
  ParamsFile = "Params.yaml";
  parse_params(ParamsFile, params);
  double c = params.bc * M_PI / ( params.Lrod * params.Lrod );
 
  //Build density matrix
  Nkx = params.Nx/2 + 1; 

  Ad3 rho( params.Nx, params.Ny, params.Nm, align);
  Ac3 rhoFT( params.Nx, params.Ny, Nkm, align);
  Ac3 rhoFTnext( params.Nx, params.Ny, Nkm, align );
  
  Ad1 fint( params.Nm,align);
  Ac1 fintFT( Nkm,align);

  // Mayer Fnc
  Ad3 Fm( params.Nx, params.Ny, params.Nm, align );
  Ac3 FmFT( params.Nx, params.Ny, Nkm, align );

  // Nonlinearities
  Ac3 NlFT( params.Nx, params.Ny, Nkm, align );
  Ac3 NlFTprev( params.Nx, params.Ny, Nkm, align );
  Ac3 MuExFT( params.Nx, params.Ny, Nkm, align );
  Ad3 ji( params.Nx, params.Ny, params.Nm, align );
  Ac3 jiFT( params.Nx, params.Ny, Nkm, align);
  Ad3 diMuTemp( params.Nx, params.Ny, params.Nm, align);
  Ac3 diMuTempFT( params.Nx, params.Ny, Nkm, align);
   
  // fftw class
  fftwpp::rcfft3d Forward3( params.Nm, rho, rhoFT );
  fftwpp::crfft3d Backward3( params.Nm, rhoFT, rho );

  // Time
  tGrid time(params.dt,params.trec,params.tend);
  
  // Build grid 
  spGrid Gridx(params.Nx,params.Lx), Gridy(params.Ny,params.Ly), Gridphi(params.Nm, 2*M_PI);
  spGrid GridGen;
  
  double*  kx =  Gridx.getKft();
  //Gridx.print();
  double*  ky =  Gridy.getKft();
  //Gridy.print(); 
  double*  km =  Gridphi.getKpos();
  Gridphi.print();

  // Initial Densities

  // Distribution
  EqDist EqVar(params.Nm,params.bc);
  EqVar.fisArrayInit(params.Nm, fint);
  // Build the density matrix
  Rho3DMaker RhoInit;
  RhoInit.BuilderEq( params.Nx, params.Ny, params.Nm,params.Lx, params.Ly, 
      c, rho, fint) ;
  Forward3.fft(rho,rhoFT);
  RhoInit.Perturb(  params.Nx, params.Ny, params.Nm, 
      params.pXmodes, params.pYmodes, params.pMmodes, params.PertrbAmp, rhoFT );
  rhoFTnext = rhoFT;
  Backward3.fftNormalized(rhoFTnext,rho);
  
  // Propagator
  Propagator DiffP( params.Nx, params.Ny, Nkm, kx, ky, km, 
      params.Dpar, params.Dr, params.dt, params.IsoDiffFlag);
  double* Lop = DiffP.getLop();
  double* Prop = DiffP.getProp();

  // Mayer function
  MayerFncHardRod(params.Nx,params.Ny,params.Nm,params.Lx,params.Ly,params.Lrod,Fm);
  Forward3.fft(Fm,FmFT);

  // NL
  NlHardRodDr Nlclass( params.Nx, params.Ny, params.Nm, kx, ky, km, params.Lx, params.Ly);
  Nlclass.NLIntCalcC( rhoFT, rho, FmFT, MuExFT, diMuTemp, diMuTempFT, ji, jiFT, 
      NlFT, Forward3, Backward3);

   
  // Write to file
  std::ofstream diffFile;
  diffFile.open ("DiffOut.txt");

  // First Step
  if( params.IsoDiffFlag = 1 ){
    switch(params.ABFlag){

      case 0: DiffP.PropAB0c( rhoFTnext, rhoFT); // Diffusion
              break;
      case 1: DiffP.PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1
              break;
      case 2: DiffP.PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1
              break;

       }
  }
  
  std::cout << "t = 0" << std::endl;
  diffFile << rho << std::endl;
  
  // Main loop

  for(int t = 1; t <= time.getNt(); ++t)
  {
    //Update
    rhoFT = rhoFTnext;
    if( params.ABFlag == 2 ){ NlFTprev = NlFT; }
    
    //Calculate NL
    Nlclass.NLIntCalcC( rhoFT, rho, FmFT, MuExFT, diMuTemp, diMuTempFT, ji, jiFT, 
      NlFT, Forward3, Backward3);


  if( params.IsoDiffFlag = 1 ){
    switch(params.ABFlag){

      case 0: DiffP.PropAB0c( rhoFTnext, rhoFT); // Diffusion
              break;
      case 1: DiffP.PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1
              break;
      case 2: DiffP.PropAB2c( rhoFTnext, rhoFT, NlFT, NlFTprev); // AB 2
              break;
       }
  }
    
    if( t % time.getNcount() == 0 )
    {
       std::cout << "t = " << t << std::endl; 
       diffFile << rho << std::endl;
    } // recording
  } // time loop 

  diffFile.close();
  return 0;
}//main

// HardRodMain.cpp
// Michael Stefferson
//
// Description: Test program of diffusion in 1D using fft
// Compile:

#include "HardRodMain.h"

int main(){
  
  fftwpp::fftw::maxthreads=get_max_threads(); //Multithreads for fft
 
  //Size of vectors to fft
  size_t align=sizeof(Complex); //Complex defined in fftw.h std::complex<double>

  //Parse parameters
  system_params params;
  params.init();
  std::string ParamsFile;
  ParamsFile = "Params.yaml"; // Should be an input
  parse_params(ParamsFile, params);
  
  // Other variables
  double c = params.bc * M_PI / ( params.Lrod * params.Lrod );
  double rhoMax = 10 * c * params.Nx * params.Ny;
  unsigned int Nkm = params.Nm/2 + 1; 
  int ShitIsFucked = 0;
  int recCounter  = 0;
  //Build density matrix

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
  double*  x =  Gridx.getPos();
//  Gridx.print();
  double*  ky =  Gridy.getKft();
  double*  y =  Gridy.getPos();
  //Gridy.print();
  double*  km =  Gridphi.getKpos();
  double*  phi =  Gridphi.getPos();
  //Gridphi.print();


  // Initial Densities

  // Distribution
  EqDist EqVar(params.Nm,params.bc);
  EqVar.fisArrayInit(params.Nm, fint);
  // Build the density matrix
  Rho3DMaker RhoInit;
  RhoInit.BuilderEq( params.Nx, params.Ny, params.Nm,params.Lx, params.Ly, 
      c, rho, fint) ;
  
  Forward3.fft(rho,rhoFT);
  RhoInit.Perturb(params.Nx, params.Ny, params.Nm, 
      params.pXmodes, params.pYmodes, params.pMmodes, params.PertrbAmp, rhoFT );
  rhoFTnext = rhoFT;
  Backward3.fftNormalized( rhoFTnext,rho );
 
  //std::cout << "rhoFT" << std::endl << rhoFT;
  /*
   *RhoInit.PerturbReal(  params.Nx, params.Ny, params.Nm, 
   *    params.pXmodes, params.pYmodes, params.pMmodes, kx, ky, km, x, y, phi,
   *    params.PertrbAmp, rho );
   *Forward3.fft(rho,rhoFT);
   *rhoFTnext = rhoFT;
   */
  
 ShitIsFucked = CheckBrokenDen( params.Nx, params.Ny, params.Nm, rho, rhoMax );
  if( (ShitIsFucked == 1) ){ std::cout << "Broke from the start!" << std::endl;}
 /*
  * 
  *std::cout << "rho =" << std::endl << rho << std::endl;
  *std::cout << "rhoFT" << std::endl << rhoFT << std::endl;
  */

  // Propagator
  Propagator DiffP( params.Nx, params.Ny, Nkm, kx, ky, km, 
      params.Dpar, params.Dr, params.dt, params.IsoDiffFlag);
//  DiffP.PropPrint();

  // Mayer function
  MayerFncHardRod(params.Nx,params.Ny,params.Nm,params.Lx,params.Ly,params.Lrod,Fm);
  Forward3.fft(Fm,FmFT);

  // NL
  NlHardRodDr Nlclass( params.Nx, params.Ny, params.Nm, kx, ky, km, 
      params.Lx, params.Ly, phi, params.vD );
  Nlclass.NLIntCalcC( rhoFT, rho, FmFT, MuExFT, diMuTemp, diMuTempFT, ji, jiFT, 
      NlFT, Forward3, Backward3);
  Nlclass.NlDrCalcC(rho, ji, jiFT, NlFT, Forward3); 

  // Order parameters
  OPclass OPs( params.Nx, params.Ny, params.Nm, phi);
  OPs.OPmaker(rho);

  std::cout << "C = " << std::endl;
  OPs.printC();
  std::cout << "PO = " << std::endl;
  OPs.printPO();
  std::cout << "NO = " << std::endl;
  OPs.printNO();
    

//// Let program know where we are
  std::cout << "Initialized all variables " << std::endl;

  //std::cout << "Nl t =0" << std::endl << NlFT << std::endl;
  //std::cout << "Fm" << std::endl << Fm << std::endl;

  /*
   *std::cout << "NlFT" << std::endl;
   *for( int k = 0; k < Nkm; ++k ) {
   *  for( int i = 0; i < params.Nx; ++i ) {
   *    for( int j = 0; j < params.Ny; ++j) {
   *      std::cout << NlFT[i][j][k] << "\t";
   *    }
   *    std::cout << std::endl;
   *  }
   *  std::cout << std::endl << std::endl;
   *}
   */

  //std::cout << "FmFT" << std::endl << FmFT << std::endl;

  // Write to file
  std::ofstream diffFile;
  diffFile.open ("DiffOut.txt");
  std::cout << "t = 0" << std::endl;
  diffFile << rho << std::endl;
  recCounter = 1;
 

  // First Step
  if( (params.IsoDiffFlag == 1) ){
    switch(params.ABFlag){

      case 0: 
        DiffP.PropAB0c( rhoFTnext, rhoFT); // Diffusion
        std::cout << "Diffusion Progator" << std::endl;
        break;
      case 1: 
        DiffP.PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1
        std::cout << "Hybrid AB 1 " << std::endl;
        break;
      case 2: 
        DiffP.PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1        
        std::cout << "Hyrid AB 2 " << std::endl;
        break;

       }
  }
  
 //std::cout << "rhoFTnext" << std::endl << rhoFTnext; 
  // Main loop

  if( ShitIsFucked == 0 ) {
    std::cout << "starting time loop" << std::endl;
    for(int t = 1; t <= time.getNt(); ++t)
    {
      //Update
      rhoFT = rhoFTnext;
      Backward3.fftNormalized(rhoFTnext,rho);
      if( (params.ABFlag == 2) ){ NlFTprev = NlFT; }
      
      //Calculate NL
      Nlclass.NLIntCalcC( rhoFT, rho, FmFT, MuExFT, diMuTemp, diMuTempFT, ji, jiFT, 
        NlFT, Forward3, Backward3);


    if( (params.IsoDiffFlag == 1) ){
      switch(params.ABFlag){

        case 0: 
          DiffP.PropAB0c( rhoFTnext, rhoFT); // Diffusion
          break;
        case 1: 
          DiffP.PropAB1c( rhoFTnext, rhoFT, NlFT ); // AB 1
          break;
        case 2: 
          DiffP.PropAB2c( rhoFTnext, rhoFT, NlFT, NlFTprev); // AB 2
          break;
         }
    }
      
      if( (t % time.getNcount() == 0) )
      {
        recCounter++;
        std::cout << "t/t_end = " << (double) t/time.getNt()  << std::endl; 
        diffFile << rho << std::endl;
        ShitIsFucked = CheckBrokenDen( params.Nx, params.Ny, params.Nm, rho, rhoMax );
        if( (ShitIsFucked == 1) ){ std::cout << "Broke!" << std::endl; break; }
      } // recording
    } // time loop 
  } // if didn't start broken

  diffFile.close();
  std::cout << "Nrec = " << time.getNrec() + 1 << std::endl;
  std::cout << "Rec Counter = " << recCounter << std::endl;
  return 0;
}//main

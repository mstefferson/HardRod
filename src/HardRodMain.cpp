// HardRodMain.cpp
// Michael Stefferson
//
// Description: Test program of diffusion in 1D using fft
// Compile:

#include "HardRodMain.h"

int main(int argc, char *argv[]){
  
  fftwpp::fftw::maxthreads = get_max_threads(); //Multithreads for fft
 
  //Size of vectors to fft
  size_t align = sizeof(Complex); //Complex defined in fftw.h std::complex<double>

  // Read input file
  std::string Inpt = argv[1];

  //Parse parameters
  system_params params;
  params.init();
  std::string ParamsFile;
  
  //ParamsFile = "Params.yaml"; // Should be an input
  ParamsFile = Inpt + ".yaml"; // Should be an input
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
  
  Ad1 fint( params.Nm, align);
  Ac1 fintFT( Nkm, align);

  // Mayer 
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
  double*  ky =  Gridy.getKft();
  double*  km =  Gridphi.getKpos();
  double*  phi =  Gridphi.getPos();

  // Initial Densities
  std::cout << "In main, going to initialize densities" << std::endl;


  // Distribution
  EqDist EqVar(params.Nm,params.bc);
  //EqVar.printFis();
  //EqVar.printFeq();
  //EqVar.fisInit(fint);
  EqVar.feqInit(fint);

  // Build the density matrix
  Rho3DMaker RhoInit;
  RhoInit.BuilderEq( params.Nx, params.Ny, params.Nm,params.Lx, params.Ly, 
      c, rho, fint);
  
  Forward3.fft(rho,rhoFT);
  RhoInit.Perturb(params.Nx, params.Ny, params.Nm, 
      params.pXmodes, params.pYmodes, params.pMmodes, params.PertrbAmp, rhoFT );
  rhoFTnext = rhoFT;
  Backward3.fftNormalized( rhoFTnext,rho );
 
 ShitIsFucked = CheckBrokenDen( params.Nx, params.Ny, params.Nm, rho, rhoMax );
  if( (ShitIsFucked == 1) ){ std::cout << "Broke from the start!" << std::endl;}

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
  // Let program know where we are
  std::cout << "Initialized all variables " << std::endl;
  

  //Use Write Class
  
  HRwriter FileWrite( params.Nx, params.Ny, params.Nm, params.trial, OPs.getC(), OPs.getPO(), OPs.getNO(),
       &( rho[params.Nx/2 + 1][params.Ny/2 + 1][0] ), 
       &(rhoFT[0][0][1]), &(rhoFT[1][0][1]), &(rhoFT[0][1][1]), &(rhoFT[1][1][1]),
       &(rhoFT[0][0][2]), &(rhoFT[1][0][2]), &(rhoFT[0][1][2]), &(rhoFT[1][1][2]) );
  FileWrite.writeOP();
  FileWrite.writeParams(params);
  //FileWrite.writeRho(rho);
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
        
        OPs.OPmaker(rho);
    //    FileWrite.writeRho(rho);
        FileWrite.writeOP();

        ShitIsFucked = CheckBrokenDen( params.Nx, params.Ny, params.Nm, rho, rhoMax );
        if( (ShitIsFucked == 1) ){ std::cout << "Broke!" << std::endl; break; }
      } // recording
    } // time loop 
  } // if didn't start broken

  std::cout << "done running" << std::endl;
  FileWrite.closeFiles();

  std::cout << "Nrec = " << time.getNrec() + 1 << std::endl;
  std::cout << "Rec Counter = " << recCounter << std::endl;
  return 0;
}//main

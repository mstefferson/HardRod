// HardRodMain.cpp
// Michael Stefferson
//
// Description: Test program of diffusion in 1D using fft
// Compile:

#include "HardRodMain.h"


int main(){
  fftwpp::fftw::maxthreads=get_max_threads(); //Multithreads for fft
  
  //Size of vectors to fft
  unsigned int Nx,Ny,Nm;
  unsigned int Nkx,Nky,Nkm;
  size_t align=sizeof(Complex);

  system_params params;

  params.init();
 
  std::cout << params.Nx << std::endl;
  //Build density matrix
  Nx = 8; Ny = 8; Nm = 8;
  Nkx = Nx/2 + 1; Nky = Ny/2 +1; Nkm = Nm / 2 + 1;
  
  Array::array3<double> c(Nx,Ny,Nm,align);
  Array::array3<Complex> cFT(Nx,Ny,Nkm,align);
  Array::array3<Complex> cFTnext(Nx,Ny,Nkm,align);
 
  fftwpp::rcfft3d Forward3(Nm,c,cFT);
  fftwpp::crfft3d Backward3(Nm,cFT,c);

  // Time
  
  double dt   = 0.001; 
  double tend = 10.00;
  double trec = 1;
  tGrid time(dt,trec,tend);
  
  // Build grid 
  double Lx, Ly; 
  
  Lx = 10.0;
  Ly = Lx;

  spGrid Gridx(Nx,Lx);
  spGrid Gridy(Ny,Ly);
  spGrid Gridphi(Nm, 2 * M_PI);
  
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

  double D = 1.0;
  bool Iso = 1;
  
  
  //Propagator DiffP(Nkx,Nky,Nkm, kx,ky,km, D, dt, Iso);
  Propagator DiffP(Nx,Ny,Nkm, kx,ky,km, D, dt, Iso);
  double* Lop = DiffP.getLop();
  double* Prop = DiffP.getProp();

  //Ml3dPrint( Lop, Nx, Ny, Nkm );
  //Ml3dPrint( Prop, Nx, Ny, Nkm );
 
 // Initialize 
int Xmodes, Ymodes, Mmodes;
Xmodes = 1;
Ymodes = 2;
Mmodes = 1;

double PerturbAmp = 0.1;


for( int ik = 0; ik < Nx; ++ik ){
  for( int jk = 0; jk < Ny; ++jk ){
    for( int kk = 0; kk < Nkm; ++kk ){
       
      if( ik <= Xmodes && jk <= Ymodes && kk <= Mmodes ){
        cFTnext[ik][jk][kk] = PerturbAmp * (Nx * Ny * Nm);
      }
      else{
        cFTnext[ik][jk][kk] = 0;
      }
    }
  }
}
  cFTnext[0][0][0] = 2.0 * (Nx*Ny*Nm);

  cFT = cFTnext;
  Backward3.fftNormalized(cFTnext,c);

   
  // Write to file
  
  std::ofstream diffFile;
  diffFile.open ("DiffOut.txt");
 

//  diffFile << cFT << std::endl;
  //std::cout << cFT << std::endl; 

// First Step

  for(int i = 0; i < Nx; ++i){
    for( int j = 0; j < Ny; ++j){
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
  
  for(int i = 0; i < Nx; ++i){
    for( int j = 0; j < Ny; ++j){
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

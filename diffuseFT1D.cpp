// diffFT.cpp
// Michael Stefferson
//
// Description: Test program of diffusion in 1D using fft
// Compile:

#include "diffuseFT1D.h"


int main(){
  fftwpp::fftw::maxthreads=get_max_threads(); //Multithreads for fft
  
  //Size of vectors to fft
  unsigned int Nx  = 16;
  unsigned int Nkx = Nx/2 + 1;
  size_t align=sizeof(Complex);
  

  Array::array1<double> c(Nx,align);
  Array::array1<Complex> cFT(Nkx,align);
  Array::array1<Complex> cFTprev(Nkx,align);
  fftwpp::rcfft1d Forward(Nx,c,cFT);
  fftwpp::crfft1d Backward(Nx,cFT,c);
 
  
  double L = 5.5;
  spGrid Grid1d(Nx,L);
  double*  k =  Grid1d.getKpos();
  double*  x = Grid1d.getPos();
  double dt = 1.1; 
  double tend = 10.0;
  double trec = 2.4;
  tGrid time(dt,trec,tend);

  double D = 1.0;
  bool Iso = 1;
  Propagator DiffP(Nkx, k, D, dt, Iso);
  
  double* Lop = DiffP.getLop();
  double* Prop = DiffP.getProp();

// Write to file
  
  std::ofstream diffFile;
  diffFile.open ("DiffOut.txt");
  
  
  for( int i = 0; i < Nx ; ++i )
  {
    c[i] = 2;
    for( int j = 1; j <= 2; ++j )
    {
    c[i] = c[i]  + 0.1 * cos( k[j] * x[i] );
    }
    diffFile << c[i] << " " << std::endl;
  }
  diffFile << std::endl;
  Forward.fft(c,cFT);

  // Main loop

  for(int t = 1; t <= time.getNt(); ++t)
  {
    cFTprev = cFT;

    for(int i = 0; i < Nkx; ++i)
    {
      cFT[i] = Prop[i] * cFT[i];
    }
    if( t % time.getNcount() == 0 )
    {
      Backward.fftNormalized(cFT,c);
      
      for( int i = 0; i < Nx; ++i )
      {
        diffFile << c[i] << " " << std::endl;
      }
       diffFile<< std::endl ;
    } // recording
} // time loop 
  diffFile.close();

  return 0;
}

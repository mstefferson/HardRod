// HRwriter.h
#ifndef _HARDROD_HRWRITER_H_
#define _HARDROD_HRWRITER_H_

#include "Parameters.h"
#include "Array.h"
#include "fftw++.h"
#include "typedef.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

class HRwriter{

  private:
    int trial_;
    int Nx_;
    int Ny_;
    int Nm_;
    int recNum_;

    std::ostringstream cTl_;
    std::ostringstream pTl_;
    std::ostringstream nTl_;
    std::ostringstream rhoTl_;
    std::ostringstream distTl_;
    std::ostringstream ampTl_;
    std::ostringstream paramsTl_;

    std::ofstream concFile_;
    std::ofstream poFile_;
    std::ofstream noFile_;
    std::ofstream rhoFile_;
    std::ofstream distFile_;
    std::ofstream ampFile_;
    std::ofstream paramsFile_;

    double** C_;
    double** NO_;
    double** PO_;
    double* rhoFixP_;
    Complex* rhoFtFix1_;
    Complex* rhoFtFix2_;
    Complex* rhoFtFix3_;
    Complex* rhoFtFix4_;
    Complex* rhoFtFix5_;
    Complex* rhoFtFix6_;
    Complex* rhoFtFix7_;
    Complex* rhoFtFix8_;

  public:
    
    HRwriter(int N, int Ny, int Nm, int trial, double** C, double** PO, double** NO, double* rhoFix, 
        Complex* rhoFtfix1, Complex* rhoFtfix2, Complex* rhoFtfix3, Complex* rhoFtfix4, 
        Complex* rhoFtfix5, Complex* rhoFtfix6, Complex* rhoFtfix7, Complex* rhoFtfix8);
    void openFiles();
    void writeOP();
    void writeC();
    void writeNO();
    void writePO();
    void writeAmp();
    void writeDist();
    void writeRho(Ad3& rho);
    void writeParams(system_params params);
    void closeFiles();


};

#endif

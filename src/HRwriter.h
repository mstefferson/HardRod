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

    std::ostringstream cTl_;
    std::ostringstream pTl_;
    std::ostringstream nTl_;
    std::ostringstream rhoTl_;
    std::ostringstream ampTl_;
    std::ostringstream paramsTl_;

    std::ofstream concFile_;
    std::ofstream poFile_;
    std::ofstream noFile_;
    std::ofstream rhoFile_;
    std::ofstream ampFile_;
    std::ofstream paramsFile_;

    double** C_;
    double** NO_;
    double** PO_;

  public:
    
    HRwriter(int Nx, int Ny, int trial, double** C, double** PO, double** NO);
    void openFiles();
    void writeOP();
    void writeC();
    void writeNO();
    void writePO();
    void writeAmp(Ac3& rhoFT);
    void writeRho(Ad3& rho);
    void writeParams(system_params params);
    void closeFiles();


};

#endif

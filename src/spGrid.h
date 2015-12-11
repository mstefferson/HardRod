//spGrid.h
#ifndef _HARDROD_SPGRID_H_
#define _HARDROD_SPGRID_H_
#include <iostream>
#include <math.h>

class spGrid
{
  public: 
    // Constructor
    spGrid();
    spGrid(int N, double L );

    void print();
    // function that prints everything
    double dxMaker();   
    double dxMaker(int N, double L);   
     // Returns grid spacing

    void xVecMaker();
    void xVecMaker(int N, double L,double* pos);
    
     // Builds position vector

    void kVecMaker();
    void kVecMaker(int N, double L, double *k);
     // Builds k vector

    void kposVecMaker();
    void kposVecMaker(int N, double *k, double *kpos );
     // Builds k >= 0 vector

    void kFtVecMaker();
    void kFtVecMaker(int N, double L, double *k);
    // builds k vecotr matching FT

    double* getPos();
     // return position vector
    
    double* getK();
     // return k vec
    
    double* getKpos();
     // return positive k values
      
    double* getKft();
    // return k FT
  private:
    int N_;         // number of gridpoints
    double L_;      // box length
    double dx_;     // grid spacing
    double *pos_;  // grid vector
    double *k_;    // k vector [ k_(-N/2), k_(-N/2+1) ... k_(0) .. k_(N/2-1) ] 
    double *kFT_;  // k vector [ k_(0) k_(1) .. k_N/2  k_(-N/2 + 1) ... k(-1) ] 
    double *kpos_; // pos k vector [ k_(0) k_(1) .. k_(N/2 = -N/2) ] 

};

#endif // header gaurd 

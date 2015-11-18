// tGrid.h

#ifndef _HARDROD_TGRID_H_
#define _HARDROD_TGRID_H_
#include <iostream>
#include <math.h>

class tGrid
{
  public: 
    // Constructor
    tGrid();
    tGrid(double dt, double t_rec, double t_end);

    void print();
    // function that prints everything
    
    void tVecMaker( int N, double tend, double* tvec);
     // Builds time vector

    double getDt();
     // return dt
     
    double getTend();
     // return t_end
    
    double getTrec();
     // return t_rec

    int getNt();
     // return total number of time points
    
    int getNrec();
     // return number of recorded points
    
    int getNcount();
     // return number of recorded points
    
    double* getTvec();
     //  return tV_

    double* getTrecVec();
    // retun trecV_

  private:
    double dt_;       // time interval 
    double tend_;    // end time
    double trec_;    //  record time
    int Nt_;          // Number of time points
    int Nrec_;        // Number of recorded points
    int Ncount_;      // Number of time points before recording
    double* tV_;      // time vector
    double* trecV_;   // recorded time vec
    
};

#endif // header gaurd 

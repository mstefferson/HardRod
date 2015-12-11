double** trapzPeriodic(double &f, double Lx, int dim, int N1, int N2, int N3){


    if( dim == 1 ){

    double F[N2][N3];
    double intFac = Lx / N1;

    for( int j == 0; j < N2; j++ ){
      for( int k == 0; k < N3; k++){
        F[j][k] = 0;
        for( int i == 0; i < N1; j ++){

          F[j][k] += f[i][j][k];

        }
        F[j][k] *= intFac;
      }
    }
  }

  if( dim == 2 ){

    double F[N1][N3];
    double intFac = Lx /  N2;

    for( int i == 0; i < N1; i++ ){
      for( int k == 0; k < N3; k++){
        F[i][k] = 0;
        for( int j == 0; j < N2; j++){

          F[i][k] += f[i][j][k];

        }
        F[i][k] *= intFac;
      }
    }
  }

   if( dim == 3 ){

    double F[N1][N2];
    double intFac = Lx /  N3;

    for( int i == 0; i < N1; i++ ){
      for( int j == 0; j < N2; j++){
        F[i][j] = 0;
        for( int k == 0; k < N3; k++){

          F[i][j] += f[i][j][k];

        }
        F[i][k] *= intFac;
      }
    }
  }
  
   return F;
}

void trapzPeriodic(double &F, double &f, double Lx, int dim, int NiInt, int Ni1, int Ni2){


  double L = x[lenX] - x[0];
  double intFac = Lx /  NiInt;


  if( dim == 1 ){
    for( int i == 0; i < Ni1; i++ ){
      for( int j == 0; j < Ni2; j++){
        F[i][j] = 0;
        for( int k == 0; k < NiInt; k ++){
          F[i][j] += f[k][i][j];
        }
        F[i][j] *= intFac;
      }
    }
  }

  if( dim == 2 ){
    for( int i == 0; i < Ni1; i++ ){
      for( int j == 0; j < Ni2; j++){
        F[i][j] = 0;
        for( int k == 0; k < NiInt; k++){
          F[i][j] += f[i][k][j];
        }
        F[i][j] *= intFac;
      }
    }
  }

   if( dim == 3 ){
     for( int i == 0; i < Ni1; i++ ){
       for( int j == 0; j < Ni2; j++){
         F[i][j] = 0;
         for( int k == 0; k < NiInt; k ++){
           F[i][j] += f[i][j][k];
         }
        F[i][j] *= intFac;
       }
     }
   }
}
 // Integrate the product of two functions, one 1D, one 3D

 void trapzPeriodicPr2 (double &F, double &f, double &g, double Lx, int dim, int NiInt, int Ni1, int Ni2){


  double intFac = Lx / NiInt;



  if( dim == 1 ){
    for( int i == 0; i < Ni1; i++ ){
      for( int j == 0; j < Ni2; j++){
        F[i][j] = 0;
        for( int k == 0; k < NiInt; k ++){
          F[i][j] += f[k][i][j] * g[k];
        }
        F[i][j] *= intFac;
      }
    }
  }

  if( dim == 2 ){
    for( int i == 0; i < Ni1; i++ ){
      for( int j == 0; j < Ni2; j++){
        F[i][j] = 0;
        for( int k == 0; k < NiInt; k++){
          F[i][j] += f[i][k][j] * g[k];
        }
        F[i][j] *= intFac;
      }
    }
  }

   if( dim == 3 ){
     for( int i == 0; i < Ni1; i++ ){
       for( int j == 0; j < Ni2; j++){
         F[i][j] = 0;
         for( int k == 0; k < NiInt; k ++){
           F[i][j] += f[i][j][k] * g[k];
         }
        F[i][j] *= intFac;
       }
     }
   }
}

 
 // Integrate the product of two functions, two 1D, one 3D
 void trapzPeriodicPr3 (double &F, double &f, double &g, double &h, double Lx, int dim, int NiInt, int Ni1, int Ni2){


  double intFac = Lx /  NiInt;


  if( dim == 1 ){
    for( int i == 0; i < Ni1; i++ ){
      for( int j == 0; j < Ni2; j++){
        F[i][j] = 0;
        for( int k == 0; k < NiInt; k ++){
          F[i][j] += f[k][i][j] * g[k] * h[k];
        }
        F[i][j] *= intFac;
      }
    }
  }

  if( dim == 2 ){
    for( int i == 0; i < Ni1; i++ ){
      for( int j == 0; j < Ni2; j++){
        F[i][j] = 0;
        for( int k == 0; k < NiInt; k++){
          F[i][j] += f[i][k][j] * g[k] * h[k];
        }
        F[i][j] *= intFac;
      }
    }
  }

   if( dim == 3 ){
     for( int i == 0; i < Ni1; i++ ){
       for( int j == 0; j < Ni2; j++){
         F[i][j] = 0;
         for( int k == 0; k < NiInt; k ++){
           F[i][j] += f[i][j][k] * g[k] * h[k];
         }
        F[i][j] *= intFac;
       }
     }
   }
}




/* the ELASTICA problem:
 -------- elasticaC.c -> elastica.dll ------
 compile in R with: system("gcc -shared -o elasticaC.dll elasticaC.c")
 or with system("R CMD SHLIB elasticaC.c") */

#include <math.h>

//  The differential system:
  void fsub(int *n, double *x, double *z, double *f,
        double * RPAR, int * IPAR)  {
        
      f[0]=cos(z[2]);
      f[1]=sin(z[2]);
      f[2]=z[3]     ;
      f[3]=z[4]*cos(z[2]);
      f[4]=0;
  }

// The analytic Jacobian for the F-function:
  void dfsub(int * n, double *x, double *z, double * df,
      double *RPAR, int *IPAR)  {

      int j;
      for (j = 0; j< *n * *n; j++) df[j] = 0;

      df[*n *2]    = -sin(z[2]);
      df[*n *2 +1] = cos(z[2]);
      df[*n *3 +2] = 1.0;
      df[*n *2 +3] = -z[4]*sin(z[2]);
      df[*n *3 +3] = 1.0;
      df[*n *4 +3] = cos(z[2]);
  }

// The boundary conditions:
  void gsub(int *i, int *n, double *z, double *g,
      double *RPAR, int *IPAR)  {
      
/*   I == the boundary condition "number".
     The conditions at the left are enumerated first,
     then the ones at the right.
     The order of left or right conditions does not matter,
     but we must be consistent when defining the jacobian
     of the boundary conditions!

     We have specified 3 left bcs, and 5 bcs total.
     This means that:
     BC(1) = x(0) = 0
     BC(2) = y(0) = 0
     BC(3) = kappa(0) = 0
     BC(4) = y(0.5) = 0
     BC(5) = phi(0.5) = -pi/2  */

      if (*i==1) *g=z[0];
      else if (*i==2) *g=z[1];
      else if (*i==3) *g=z[3];
      else if (*i==4) *g=z[1];
      else if (*i==5) *g=z[2]+1.5707963267948966192313216916397514;
  }

// The analytic Jacobian for the G-function:
  void dgsub(int *i, int *n, double *z, double *dg,
      double *RPAR, int *IPAR)  {

      int j;
      for (j = 0; j< *n; j++) dg[j] = 0;

      if (*i == 1) dg[0] = 1.;
      else if (*i == 2) dg[1] = 1.;
      else if (*i == 3) dg[3] = 1.;
      else if (*i == 4) dg[1] = 1.;
      else if (*i == 5) dg[2] = 1.;
  }

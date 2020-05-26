#include <time.h>
#include <string.h>
#include "bvpSolve.h"
#include "externalptr.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   boundary value problem solvers.
   
   The C-wrappers that provide the interface between FORTRAN codes and R-code 
   are: C_bvp_deriv_func   : interface with R-code "derivfunc", passes derivatives  
        C_bvp_jac_func     : interface with R-code "jacfunc", passes jac_funcobian
        C_bvp_bound_func   : interface with R-code "bound_func", boundaries
        C_bvp_jacbound_func: interface with R-code "jacbound_func", jac_funcobian  of boundaries
        C_bvp_guess_func   : interface with R-code "guess_func", initial estimates
    November 2011: added coldae
		
		
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* definition of the calls to the fortran functions -

      Subroutine Colnew(Ncomp, M, Aleft, Aright, Zeta, Iset, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, 
     +     Fsub, Dfsub, Gsub, Dgsub, guess_func)
     
     
      Subroutine Appsln(Xx,Z,Fspace,Ispace)               */

void F77_NAME(colnew)(int*, int*, double *, double *, double *, int *, int *,
         double *, double *, int *, double *, int *, 
         void (*)(int *, double *, double *, double *, double *, int *),   /* fsub  */
		     void (*)(int *, double *, double *, double *, double *, int *),   /* dfsub */
			   void (*)(int *, int *, double *, double *, double *, int *),      /* gsub  */
		     void (*)(int *, int *, double *, double *, double *, int *),      /* dgsub */
         void (*)(double *, double *, double *, double *, int *),                           /* guess_func */
         double *, int *, int*) ;

void F77_NAME(colsys)(int*, int*, double *, double *, double *, int *, int *,
         double *, double *, int *, double *, int *, 
         void (*)(int *, double *, double *, double *, double *, int *),   /* fsub  */
		     void (*)(int *, double *, double *, double *, double *, int *),   /* dfsub */
			   void (*)(int *, int *, double *, double *, double *, int *),      /* gsub  */
		     void (*)(int *, int *, double *, double *, double *, int *),      /* dgsub */
         void (*)(double *, double *, double *, double *, int *),                           /* guess_func */
         double *, int *, int*) ;

void F77_NAME(coldae)(int*, int*, int *, double *, double *, double *, 
         int *, int *, double *, double *, int *, double *, int *, 
         void (*)(int *, double *, double *, double *, double *, double *, int *),   /* fsub  */
		     void (*)(int *, double *, double *, double *, double *, double *, int *),   /* dfsub */
			   void (*)(int *, int *, double *, double *, double *, int *),              /* gsub  */
		     void (*)(int *, int *, double *, double *, double *, int *),                /* dgsub */
         void (*)(double *, double *, double *, double *, double *, int *),              /* guess_func */
         double *, int *, int*) ;

void F77_NAME(appsln)(double *, double *, double *, int *);
void F77_NAME(sysappsln)(double *, double *, double *, int *);   /*also defined in init..?*/  
void F77_NAME(appsln_dae)(double *, double *, double *, double *, int *);

typedef void C_deriv_func_DAE_type    (int *, double *, double *, double *, double *, double *, int *);
typedef void C_jac_func_DAE_type      (int *, double *, double *, double *, double *, double *, int *);
typedef void C_guess_func_type2       (double *, double *, double *, double *, int *);
typedef void C_guess_func_DAE_type    (double *, double *, double *, double *, double *, int *);
C_deriv_func_type  *derfun_DAE = NULL;
C_jac_func_type *jacfundae = NULL;
C_deriv_func_DAE_type  *jderfundae = NULL;
int nalg;

/* KARLINE -> FRANCESCA: NUMERICAL FUNCTIONS, IF NOT GIVEN */

static void C_num_jac_func (int *n,  double *x, double *y, double *pd,
                            double * rpar, int * ipar)
{
  int i, j;
  double perturb;

  for (i = 0; i < mstar; i++) ycopy[i]   = y[i];

  jderfun(n, x, y, dy, rpar, ipar);
  for (i = 0; i < n_eq; i++) dycopy[i]  = dy[i];
  
  for (j = 0; j < mstar; j++) {
     if (y[j] > 1.)
       perturb = y[j]*1e-8;
     else
       perturb = 1e-8 ;  
     ycopy[j] = y[j] + perturb;
     
     jderfun(n, x, ycopy, dycopy, rpar, ipar);
     
     ycopy[j] = y[j];
     
     for (i = 0; i < n_eq; i++) 
       pd[j* n_eq + i] = (dycopy[i] - dy[i])/perturb;

  }
}

static void C_num_jac_func_DAE (int *n,  double *x, double *y, double *y2, double *pd,
                            double * rpar, int * ipar)
{
  int i, j;
  double perturb;

  for (i = 0; i < mstar - nalg; i++) ycopy[i]   = y[i];
  for (i = 0; i < nalg; i++) ycopy2[i]          = y2[i];

  jderfundae(n, x, y, y2, dy, rpar, ipar);
  for (i = 0; i < n_eq; i++) dycopy[i]  = dy[i];
  
  for (j = 0; j < mstar - nalg; j++) {
     
     if (y[j] > 1.)
       perturb = y[j]*1e-8;
     else
       perturb = 1e-8 ;  
     ycopy[j] = y[j] + perturb;
     
     jderfundae(n, x, ycopy, y2, dycopy, rpar, ipar);
     
     ycopy[j] = y[j];
     
     for (i = 0; i < n_eq; i++) 
       pd[j* n_eq + i] = (dycopy[i] - dy[i])/perturb;

  }
  
  for (j = 0; j < nalg; j++) {
     
     if (y2[j] > 1.)
       perturb = y2[j]*1e-8;
     else
       perturb = 1e-8 ;  
     ycopy2[j] = y2[j] + perturb;
     
     jderfundae(n, x, y, ycopy2, dycopy, rpar, ipar);
     
     ycopy2[j] = y2[j];
     
     for (i = 0; i < n_eq; i++) 
       pd[(j + mstar - nalg)* n_eq + i] = (dycopy[i] - dy[i])/perturb;

  }
}

static void C_num_bound_func (int *ii, int *n, double *y, double *gout,
                              double * rpar, int * ipar)
{
   int i, ib;
   i =  ii[0] - 1;        /*-1 to go from R to C indexing*/

   ib = iibb[i] - 1;

   gout[0] = y[ib] - bb[i];       
}

static void C_num_jacbound_func (int *ii, int *n, double *y, double *dg,
                                 double * rpar, int * ipar)
{
  int i;
  double perturb;

  for (i = 0; i < mstar-nalg; i++) ycopy[i]  = y[i];

  for (i = 0; i < mstar-nalg; i++) {
    jbndfun(ii, n, y, g, rpar, ipar);

    if (y[i] > 1.)
       perturb = y[i]*1e-8;
    else
       perturb = 1e-8;  

    ycopy[i] = y[i] + perturb;
    jbndfun(ii, n, ycopy, gcopy, rpar, ipar);
    ycopy[i] = y[i];
    dg[i] = (gcopy[0] - g[0])/perturb;
  }
}


/* -----------------------------------------------------------------------------
                        - when model in compiled code
----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
  wrapper above the derivate function that first estimates the
values of the forcing functions */

static void dll_bvp_deriv_func_forc (int *neq, double *x, double *y,
                         double *ydot, double *rpar, int *ipar)
{
  updatedeforc(x);
  derfun(neq, x, y, ydot, rpar, ipar);
}

static void dll_bvp_deriv_func_DAE_forc (int *neq, double *x, double *y, double *y2,
                         double *ydot, double *rpar, int *ipar)
{
  int i;
  updatedeforc(x);
  for (i = 0; i < mstar - nalg; i++) ycopy[i] = y[i];
  for (i = 0; i < nalg; i++) ycopy[mstar - nalg + i] = y2[i];
  
  derfun_DAE(neq, x, ycopy, ydot, rpar, ipar);
}                       

static void dll_bvp_jac_func_DAE_forc (int *neq, double *x, double *y, double *y2,
                         double *ydot, double *rpar, int *ipar)
{
  int i;
  updatedeforc(x);
  for (i = 0; i < mstar - nalg; i++) ycopy[i] = y[i];
  for (i = 0; i < nalg; i++) ycopy[mstar - nalg + i] = y2[i];
  
  jacfundae(neq, x, ycopy, ydot, rpar, ipar);
} 
                      
static void wrap_bvp_deriv_func_DAE (int *neq, double *x, double *y, double *y2,
                         double *ydot, double *rpar, int *ipar)
{
  int i;
  for (i = 0; i < mstar - nalg; i++) ycopy[i] = y[i];
  for (i = 0; i < nalg; i++) ycopy[mstar - nalg + i] = y2[i];
  
  derfun_DAE(neq, x, ycopy, ydot, rpar, ipar);   /* user DLL */
}                       


static void wrap_bvp_jac_func_DAE (int *neq, double *x, double *y, double *y2,
                         double *pd, double *rpar, int *ipar)
{
  int i;
  for (i = 0; i < mstar - nalg; i++) ycopy[i] = y[i];
  for (i = 0; i < nalg; i++) ycopy[mstar - nalg + i] = y2[i];
  
  jacfundae(neq, x, ycopy, pd, rpar, ipar); /* user-defined*/
}                       

/* -----------------------------------------------------------------------------
   interface between fortran function calls and R functions
   Note: passing of parameter values and "..." is done in R-function bvpcol
----------------------------------------------------------------------------- */

/* initialisation function                                                    */
static void C_bvp_guess_func (double *x, double *y,  double *ydot,
                              double *rpar, int *ipar)
{
  int i;
  double p;
  SEXP R_fcall, X, ans, R_fcall2, ans2;

  PROTECT(X = ScalarReal(*x));                      
  PROTECT(R_fcall = lang2(R_bvp_guess_func, X));    
  PROTECT(ans = eval(R_fcall, R_envir));            

  p = fmax(1e-7, *x*1e-7 );
  REAL(X)[0]   = *x+p;
  PROTECT(R_fcall2 = lang2(R_bvp_guess_func, X));   
  PROTECT(ans2 = eval(R_fcall2, R_envir));          

  /* both have the same dimensions... */
  for (i = 0; i < n_eq; i++) y[i] = REAL(ans)[i];
  for (i = 0; i < n_eq; i++) ydot[i] = (REAL(ans2)[i]-y[i])/p;

  UNPROTECT(5);

}

/* derivative function                                                        */
static void C_bvp_deriv_func (int * n, double *x, double *y,
                              double *ydot, double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, X, ans;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                    
  PROTECT(R_fcall = lang3(R_bvp_deriv_func,X,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));          

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  UNPROTECT(3);
}

/* jacobian                                                                   */
static void C_bvp_jac_func (int *n, double *x, double *y, double *pd,
                            double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, X, ans;
  
  for (i = 0; i < mstar; i++) REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                  
  PROTECT(R_fcall = lang3(R_bvp_jac_func,X,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));        

  for (i = 0; i < n_eq * mstar; i++)  pd[i] = REAL(ans)[i];
  UNPROTECT(3);
}


/* derivative function                                                        */
static void C_bvp_deriv_func_DAE (int * n, double *x, double *y, double *y2,
                              double *ydot, double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, X, ans;

  for (i = 0; i < mstar-nalg ; i++)  REAL(Y)[i]    = y[i];
  for (i = 0; i < nalg; i++) REAL(Y)[i+mstar-nalg] = y2[i];

  PROTECT(X = ScalarReal(*x));                    
  PROTECT(R_fcall = lang3(R_bvp_deriv_func,X,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));          

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  UNPROTECT(3);
}

/* jacobian                                                                   */
static void C_bvp_jac_func_DAE (int *n, double *x, double *y, double *y2, double *pd,
                            double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, X, ans;

  for (i = 0; i < mstar-nalg; i++) REAL(Y)[i]      = y[i];
  for (i = 0; i < nalg; i++) REAL(Y)[i+mstar-nalg] = y2[i];

  PROTECT(X = ScalarReal(*x));                    
  PROTECT(R_fcall = lang3(R_bvp_jac_func,X,Y));   
  PROTECT(ans = eval(R_fcall, R_envir));          
  
  for (i = 0; i < n_eq * mstar; i++)  pd[i] = REAL(ans)[i];

  UNPROTECT(3);
}

/* initialisation function                                                    */
static void C_bvp_guess_func_DAE (double *x, double *y, double *y2, double *ydot,
                              double *rpar, int *ipar)
{
  int i;
  double p;
  SEXP R_fcall, X, ans, R_fcall2, ans2;

  PROTECT(X = ScalarReal(*x));                      
  PROTECT(R_fcall = lang2(R_bvp_guess_func, X));    
  PROTECT(ans = eval(R_fcall, R_envir));            

  p = fmax(1e-7, *x*1e-7 );
  REAL(X)[0]   = *x+p;
  PROTECT(R_fcall2 = lang2(R_bvp_guess_func, X));   
  PROTECT(ans2 = eval(R_fcall2, R_envir));          

  /* both have the same dimensions... */
  for (i = 0; i <  mstar-nalg; i++) y[i] = REAL(ans)[i];
  for (i = 0; i < nalg; i++)       y2[i] = REAL(ans)[i+mstar-nalg] ;

  for (i = 0; i <  mstar-nalg; i++) ydot[i]           = (REAL(ans2)[i]-y[i])/p;
  for (i = 0; i < nalg; i++)       ydot[i+mstar-nalg] = (REAL(ans2)[i+mstar-nalg]-y2[i])/p;

  UNPROTECT(5);

}

/*  boundary condition                                                        */
static void C_bvp_bound_func (int *ii, int * n, double *y, double *gout,
                              double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, J, ans;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                  
  PROTECT(R_fcall = lang3(R_bvp_bound_func,J,Y));   
  PROTECT(ans = eval(R_fcall, R_envir));            
  /* only one element returned... */
  gout[0] = REAL(ans)[0];

  UNPROTECT(3);
}

/* jacobian of boundary condition                                             */
static void C_bvp_jacbound_func (int *ii, int *n, double *y, double *dg,
                                 double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, J, ans;
  for (i = 0; i < mstar; i++) REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                   
  PROTECT(R_fcall = lang3(R_bvp_jacbound_func,J,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));             

  for (i = 0; i < mstar ; i++)  dg[i] = REAL(ans)[i];
  UNPROTECT(3);
}

/*  boundary condition                                                        */
static void C_bvp_bound_func_DAE (int *ii, int * n, double *y, double *gout,
                              double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, J, ans;

  for (i = 0; i < mstar -nalg; i++)  REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                  
  PROTECT(R_fcall = lang3(R_bvp_bound_func,J,Y));   
  PROTECT(ans = eval(R_fcall, R_envir));            
  /* only one element returned... */
  gout[0] = REAL(ans)[0];

  UNPROTECT(3);
}

/* jacobian of boundary condition                                             */
static void C_bvp_jacbound_func_DAE (int *ii, int *n, double *y, double *dg,
                                 double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, J, ans;

  for (i = 0; i < mstar-nalg; i++) REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                   
  PROTECT(R_fcall = lang3(R_bvp_jacbound_func,J,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));             

  for (i = 0; i < mstar-nalg ; i++)  dg[i] = REAL(ans)[i];
  UNPROTECT(3);
}


/* -----------------------------------------------------------------------------
                  MAIN C-FUNCTION, CALLED FROM R-code
----------------------------------------------------------------------------- */

SEXP call_colnew(SEXP Ncomp, SEXP Xout, SEXP Aleft, SEXP Aright,
		SEXP Zeta, SEXP Mstar, SEXP M, SEXP Iset, SEXP Rwork, SEXP Iwork,
    SEXP Tol, SEXP Fixpnt, SEXP Rpar, SEXP Ipar,
		SEXP derivfunc, SEXP jacfunc, SEXP boundfunc,
    SEXP jacboundfunc, SEXP guessfunc, SEXP Initfunc, SEXP Parms, SEXP flist,
    SEXP Type, SEXP Absent, SEXP RRwork, SEXP rho)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, ICOUNT, RWORK;

  int  j, ii, ncomp, k, nx, ntol, nfixpnt, isForcing, type;
  double aleft, aright, *zeta, *fspace, *tol, *fixpnt, *z, *yz, *rpar;
  double xout;
  int *m, *ispace, *iset, *icount, *ltol, *ipar, iflag, isDll, FullOut;
  int *absent;
  double *rwork;

  C_deriv_func_type        *deriv_func = NULL;
  C_jac_func_type          *jac_func = NULL;
  C_deriv_func_DAE_type    *deriv_func_DAE = NULL;
  C_jac_func_DAE_type      *jac_func_DAE = NULL;
  C_bound_func_type        *bound_func = NULL;
  C_jacbound_func_type     *jacbound_func = NULL;
  C_guess_func_type2       *guess_func = NULL;
  C_guess_func_DAE_type    *guess_func_DAE = NULL;
  
/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  int nprot = 0;

  aleft  =REAL(Aleft)[0];
  aright =REAL(Aright)[0];

  ncomp = INTEGER(Ncomp)[0];     /* number of equations */
  type  = INTEGER(Type)[0];      /* 0=colnew, 1 = colsys, 2 = coldae */
  
  ii = LENGTH(Absent);
  absent = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) absent[j] = INTEGER(Absent)[j];

  ii = LENGTH(RRwork);
  rwork = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rwork[j] = REAL(RRwork)[j];


  n_eq  = INTEGER(Ncomp)[0];     /* number of equations -global variable */
  mstar = INTEGER(Mstar)[0];     /* number of variables */
  nalg = 0; /* Francesca Mazzia: needed by C_num_jacbound_func */
  if (type == 2) {
    nalg = INTEGER(Mstar)[1]; 
    ncomp = ncomp  - nalg;       /* number of differential equations */
  }
  /* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  m  = (int *) R_alloc(n_eq, sizeof(int));  /* order of diff eqns */
  for (j = 0; j < ncomp; j++) m[j] = INTEGER(M)[j];


  ii = LENGTH(Zeta);
  zeta   =(double *) R_alloc(ii, sizeof(double));
  for (j = 0; j < ii;j++) zeta[j] = REAL(Zeta)[j];

  ii = LENGTH(Iset) -1;    /* length of Iset, integer settings */
  iset  = (int *)    R_alloc(ii, sizeof(int));
  for (j = 0; j < ii; j++) iset[j] = INTEGER(Iset)[j];

  icount  = (int *)    R_alloc(6, sizeof(int));
  for (j = 0; j < 6; j++) icount[j] = 0;
  
  FullOut = INTEGER(Iset)[ii];
  
  ntol = LENGTH(Tol);
  tol   =(double *) R_alloc(ntol, sizeof(double));
  for (j = 0; j < ntol; j++) tol[j] = REAL(Tol)[j];
    
  ltol   =(int *) R_alloc(ntol, sizeof(int));
  for (j = 0; j < ntol; j++) ltol[j] = j+1;

  nfixpnt =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(nfixpnt, sizeof(double));
  for (j = 0; j < nfixpnt; j++) fixpnt[j] = REAL(Fixpnt)[j];

/*      error("mstar, ncomp,ltol %i %i %i %f",ncomp,mstar,ltol,tol[0]);*/

  ii = iset[5];
  ispace = (int *) R_alloc(ii, sizeof(int));
  for (j = 0; j < ii; j++) ispace[j] = 0;

  ii = iset[4];
  fspace = (double *) R_alloc(ii, sizeof(double));

  if (iset[8] > 1) {   /* continuation request: values in Rwork, Iwork*/
    ii = LENGTH(Rwork);
    for (j=0; j < ii; j++) fspace[j] = REAL(Rwork)[j];
    ii = LENGTH(Iwork);
    for (j=0; j < ii; j++) ispace[j] = INTEGER(Iwork)[j];
  }
  ii = LENGTH(Ipar);
  ipar = (int *) R_alloc(ii, sizeof(int));
  for (j=0; j<ii; j++) ipar[j] = INTEGER(Ipar)[j];

  ii = LENGTH(Rpar);
  rpar = (double *) R_alloc(ii, sizeof(double));
  for (j=0; j<ii; j++) rpar[j] = REAL(Rpar)[j];

/* initialise global R-variables... */

  if (isDll == 0) {
    PROTECT(Y = allocVector(REALSXP,mstar));    nprot++;
  }
  /* Initialization of Parameters and Forcings (DLL functions)  */
  isForcing = initForcings(flist);
  //initParms(Initfunc, Parms);
  if (Initfunc != NA_STRING) {
    if (inherits(Initfunc, "NativeSymbol"))  {
      init_func *initializer;
      
      PROTECT(bvp_gparms = Parms);     nprot++;
      initializer = (init_func *) R_ExternalPtrAddrFn_(Initfunc);
      initializer(Initbvpparms);
    }
  }
  
  R_envir = rho;
  ycopy  = (double *) R_alloc(mstar, sizeof(double)); 

  /* pointers to functions passed to FORTRAN */
  if (isDll) {   /* DLL addresses passed to fortran */
     if (type == 2) {
      deriv_func_DAE = (C_deriv_func_DAE_type *) wrap_bvp_deriv_func_DAE;
      derfun_DAE     = (C_deriv_func_type *)     R_ExternalPtrAddrFn_(derivfunc);
     } else { 
      deriv_func     = (C_deriv_func_type *)     R_ExternalPtrAddrFn_(derivfunc);
     }
     if (absent[0] == 0) { 
       if (type == 2) {
         jac_func_DAE = (C_jac_func_DAE_type *) wrap_bvp_jac_func_DAE;
         jacfundae = (C_jac_func_type *)       R_ExternalPtrAddrFn_(jacfunc);
       } else {
         jac_func      = (C_jac_func_type *)   R_ExternalPtrAddrFn_(jacfunc);
        }
     }
     
     if (absent[1] == 0)
        bound_func    = (C_bound_func_type *)  R_ExternalPtrAddrFn_(boundfunc);

     if (absent[2] == 0)   /* not given*/
        jacbound_func = (C_jacbound_func_type *) R_ExternalPtrAddrFn_(jacboundfunc);

	  /* here overruling deriv_func if forcing */
     if (isForcing) {
        if (type ==2) {
          deriv_func_DAE = dll_bvp_deriv_func_DAE_forc;
        if (absent[0] == 0) 
          jac_func_DAE   = (C_jac_func_DAE_type *) dll_bvp_jac_func_DAE_forc;
        } else {
          derfun =     (C_deriv_func_type *) R_ExternalPtrAddrFn_(derivfunc);
          deriv_func = (C_deriv_func_type *) dll_bvp_deriv_func_forc;
        }
      }

  } else {      /* interface functions between fortran and R */
     
    if (type ==2) 
      deriv_func_DAE = C_bvp_deriv_func_DAE;
    else
      deriv_func = C_bvp_deriv_func;
     
    R_bvp_deriv_func = derivfunc;

    if (absent[0] == 0) {
     
      if (type ==2) 
       jac_func_DAE = C_bvp_jac_func_DAE;
      else
       jac_func = C_bvp_jac_func;

      R_bvp_jac_func = jacfunc;
    }
    
    if (absent[1] == 0) {

      if (type ==2) 
       bound_func = C_bvp_bound_func_DAE;
      else
       bound_func = C_bvp_bound_func;
     
      R_bvp_bound_func = boundfunc;
    }

    if (absent[2] == 0) {
      if (type ==2) 
       jacbound_func = C_bvp_jacbound_func_DAE;
      else
       jacbound_func = C_bvp_jacbound_func;
      R_bvp_jacbound_func = jacboundfunc;
    }
  }

/* if numerical approximates should be used */    
    if (absent[0] == 1) {
        dy     = (double *) R_alloc(n_eq, sizeof(double));
        dycopy = (double *) R_alloc(n_eq, sizeof(double));
      if (type ==2) {
       jac_func_DAE = C_num_jac_func_DAE;
       jderfundae  = deriv_func_DAE;
       ycopy2  = (double *) R_alloc(mstar, sizeof(double)); 

      } else {
        jac_func = (C_jac_func_type *)      C_num_jac_func;
        jderfun  = deriv_func;
      }
    }
    
    if (absent[1] == 1) {
      bound_func = (C_bound_func_type *) C_num_bound_func;
      iibb = (int *) R_alloc(mstar-nalg, sizeof(int));
      for (j = 0; j < mstar-nalg; j++)
        iibb[j] = absent[3 + j];
      bb = (double *) R_alloc(mstar-nalg, sizeof(double));
      for (j = 0; j < mstar-nalg; j++)
        bb[j] = rwork[j];
    }
    
    if (absent[2] == 1) {
      jacbound_func = (C_jacbound_func_type *) C_num_jacbound_func;
      jbndfun = bound_func;
      g = (double *) R_alloc(1, sizeof(double));
      gcopy = (double *) R_alloc(1, sizeof(double));
      if (absent[0] != 1)
        ycopy  = (double *) R_alloc(mstar, sizeof(double)); 
    }
    

     if (type ==2) 
      guess_func_DAE = (C_guess_func_DAE_type *) C_bvp_guess_func_DAE;
     else
      guess_func = (C_guess_func_type2 *) C_bvp_guess_func;
	  R_bvp_guess_func = guessfunc;
      
/* Call the fortran function -
      Subroutine colnew(Ncomp, M, Aleft, Aright, Zeta, Iset, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, 
     +     Fsub, Dfsub, Gsub, Dgsub, guess_func)             */
   if (type == 0) 
	  F77_CALL(colnew) (&ncomp, m, &aleft, &aright, zeta, iset, ltol,
        tol, fixpnt, ispace, fspace, &iflag, 
        deriv_func, jac_func, bound_func, jacbound_func, 
        guess_func, rpar, ipar, icount);
   else if (type ==1)
	  F77_CALL(colsys) (&ncomp, m, &aleft, &aright, zeta, iset, ltol,
        tol, fixpnt, ispace, fspace, &iflag, 
        deriv_func, jac_func, bound_func, jacbound_func, 
        guess_func, rpar, ipar, icount);
   else
	  F77_CALL(coldae) (&ncomp, &nalg, m, &aleft, &aright, zeta, iset, ltol,
        tol, fixpnt, ispace, fspace, &iflag, 
        deriv_func_DAE, jac_func_DAE, bound_func, jacbound_func, 
        guess_func_DAE, rpar, ipar, icount);
   
/*             Call Appsln(Xx,Z,Fspace,Ispace)
C....   Iflag - The Mode Of Return From colnew/colsys.
C....         =  1  For Normal Return
C....         =  0  If The Collocation Matrix Is Singular For The Final
C....               Continuation Problem.
C....         = -1  If The Expected No. Of Subintervals Exceeds Storage
C....               Specifications.
C....         = -2  If The Nonlinear Iteration Has Not Converged For The
C....               Final Continuation Problem.
C....         = -3  If There Is An Input Data Error.

          */

  if (iflag == 0)
	{
    UNPROTECT(nprot);
    error("The collocation matrix is singular for the final continuation problem\n");
	}
  else if (iflag == -1)
	{
    UNPROTECT(nprot);
    error("The Expected No. Of Subintervals Exceeds Storage Specifications.\n");
	}
  else if (iflag == -2)
	{
    UNPROTECT(nprot);
    error("The Nonlinear Iteration Has Not Converged For The Final Continuation Problem.\n");
	}
  else  if (iflag == -3)
	{
    UNPROTECT(nprot);
    error("Illegal input to bvpcol\n");
	}
  else
	{
    nx = LENGTH(Xout);
    z  =(double *) R_alloc(mstar - nalg, sizeof(double));

    PROTECT(yout = allocMatrix(REALSXP,mstar+1,nx)); nprot++;
   if (type == 0) 
	  for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(appsln)(&xout,z,fspace,ispace);
                 for (j=0;j<mstar;j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
      }  /* end main x loop */
    else if (type == 1)
	   for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(sysappsln)(&xout,z,fspace,ispace);
                 for (j=0;j<mstar;j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
      }  /* end main x loop */
    else {
	   yz = (double *) R_alloc(nalg, sizeof(double));
	   for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(appsln_dae)(&xout, z, yz, fspace, ispace);
                 for (j=0;j < mstar-nalg; j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
                 for (j=0;j < nalg; j++) REAL(yout)[k*(mstar+1) + mstar -nalg +  j+1] = yz[ j];
      }  /* end main x loop */
    }
  if (type <= 1)
     ii = ncomp+7;
  else
	 ii = ncomp+8;
  PROTECT(ICOUNT = allocVector(INTSXP, 6));    nprot++;
  PROTECT(ISTATE = allocVector(INTSXP, ii+6)); nprot++;
  INTEGER(ISTATE)[0] = iflag;
  for (k = 0; k < 6; k++)  INTEGER(ICOUNT)[k] = icount[k];
  for (k = 0; k < 5; k++)  INTEGER(ISTATE)[k+1] = icount[k];
  for (k = 0; k < ii; k++) INTEGER(ISTATE)[k+6] = ispace[k];
  if (FullOut) 
	if (type == 2)
       ii = ispace[7];
	else
	   ii = ispace[6];
  else 
    ii = 1;

  PROTECT(RWORK = allocVector(REALSXP, ii));   nprot++;
  for (k = 0; k<ii; k++) REAL(RWORK)[k] = fspace[k];
  setAttrib(yout, install("icount"), ICOUNT);
  setAttrib(yout, install("istate"), ISTATE);
  setAttrib(yout, install("rstate"), RWORK); 
  }
/*                       ####   termination   ####                            */    
  UNPROTECT(nprot);
  return(yout);
}

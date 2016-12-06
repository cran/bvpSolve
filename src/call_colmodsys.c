#include <time.h>
#include <string.h>
#include "bvpSolve.h"
#include "externalptr.h"

/* -----------------------------------------------------------------------------
  definition of the calls to the fortran functions -
      Subroutine Colmod(Ncomp, M, Aleft, Aright, Zeta, Ipar, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, Eps, Epsmin,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess)
                                                                              */
void F77_NAME(colmod)(int*, int*, double *, double *, double *, int *, int *,
 double *, double *, int *, double *, int *, double *, double *,
 void (*)(double *, double *, double *, double *, double *, int *),        /* fsub  */
 void (*)(double *, double *, double *, int *, double *, double *, int *), /* dfsub */
 void (*)(int *, double *, double *, double *, double *, int *),           /* gsub  */
 void (*)(int *, double *, double *, double *, double *, int *),           /* dgsub */
 void (*)(double *, double *, double *, double *, double *, int *),        /* guess */
         double *, int *, int*);      

/*    Subroutine Appsln(Xx,Z,Fspace,Ispace)                                   */
void F77_NAME(mappsln)(double *, double *, double *, int *);


static void C_num_epsjac_func (double *x, double *y, double *pd,
                            int *n,  double * eps, double * rpar, int * ipar)
{
  int i, j;
  double perturb;

  for (i = 0; i < mstar; i++) ycopy[i]   = y[i];

  jepsderfun(x, y, dy, eps, rpar, ipar);
  for (i = 0; i < n_eq; i++) dycopy[i]  = dy[i];
  
  for (j = 0; j < mstar; j++) {
     if (y[j] > 1.)
       perturb = y[j]*1e-8;
     else
       perturb = 1e-8 ;  
     ycopy[j] = y[j] + perturb;
     
     jepsderfun(x, ycopy, dycopy, eps, rpar, ipar);
     
     ycopy[j] = y[j];
     
     for (i = 0; i < n_eq; i++) 
       pd[j* n_eq + i] = (dycopy[i] - dy[i])/perturb;

  }
}

static void C_num_epsbound_func (int *ii, double *y, double *gout,
                              double * eps, double * rpar, int * ipar)
{
   int i, ib;
   
   i =  ii[0] - 1;        /*-1 to go from R to C indexing*/

   ib = iibb[i] - 1;

   gout[0] = y[ib] - bb[i];       
}

static void C_num_epsjacbound_func (int *ii, double *y, double *dg,
                                 double * eps, double * rpar, int * ipar)
{
  int i;
  double perturb;

  for (i = 0; i < mstar; i++) ycopy[i]  = y[i];

  for (i = 0; i < mstar; i++) {
    jepsbndfun(ii, y, g, eps, rpar, ipar);

    if (y[i] > 1.)
       perturb = y[i]*1e-8;
    else
       perturb = 1e-8;  

    ycopy[i] = y[i] + perturb;
    jepsbndfun(ii, ycopy, gcopy, eps, rpar, ipar);
    ycopy[i] = y[i];
    dg[i] = (gcopy[0] - g[0])/perturb;
  }
}


/* -----------------------------------------------------------------------------
                        - when model in compiled code
----------------------------------------------------------------------------- */

/* wrapper above the derivate function that first estimates the
values of the forcing functions and puts new value of eps in par
The interface is slightly different in colmod - overruled in C-code           */

static void dll_colmod_deriv_func_forc (double *x, double *y,
                         double *ydot, double *eps, double *rpar, int *ipar)
{
  updatedeforc(x);
  epsval[0] = eps[0];
  rpar[ipar[0]] = eps[0];
  derfun(&n_eq, x, y, ydot, rpar, ipar);
}

static void dll_colmod_deriv_func (double *x, double *y,
                         double *ydot, double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];       /* value of parameter */
  rpar[ipar[0]] = eps[0];   /* value of parameter */
  derfun(&n_eq, x, y, ydot, rpar, ipar);
}

static void dll_colmod_jac_func (double *x, double *y,
                         double *pd, int * n, double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  rpar[ipar[0]] = eps[0];   /* value of parameter */
  jacfun(n, x, y, pd, rpar, ipar);
}

static void dll_colmod_bound_func (int *ii, double *y, double *gout,
                        double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  rpar[ipar[0]] = eps[0];   /* value of parameter */
  boundfun(ii, &n_eq, y, gout, rpar, ipar);
}

static void dll_colmod_jacbound_func (int *ii, double *y, double *dg,
                        double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  rpar[ipar[0]] = eps[0];     /* value of parameter */
  jacboundfun(ii, &n_eq, y, dg, rpar, ipar);
}

/* -----------------------------------------------------------------------------
   interface between fortran function calls and R functions
   Note: passing of parameter values and "..." is done in R-function bvpcol
----------------------------------------------------------------------------- */

/* initialisation function                                                    */
static void C_colmod_guess (double *x, double *y,
             double *ydot, double *eps, double *rpar, int *ipar)
{
  int i;
  double p;
  SEXP R_fcall, X, ans, R_fcall2, ans2;

  PROTECT(X = ScalarReal(*x));                     incr_N_Protect();
  PROTECT(R_fcall = lang2(R_cont_guess_func, X));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));           incr_N_Protect();

  p = fmax(1e-7, *x*1e-7 );
  REAL(X)[0]   = *x+p;
  PROTECT(R_fcall2 = lang2(R_cont_guess_func, X)); incr_N_Protect();
  PROTECT(ans2 = eval(R_fcall2, R_envir));         incr_N_Protect();
  
  /* both have the same dimensions... */
  for (i = 0; i < n_eq; i++) y[i] = REAL(ans)[i];
  for (i = 0; i < n_eq; i++) ydot[i] = (REAL(ans2)[i]-y[i])/p;

  my_unprotect(5);
}

/* derivative function                                                        */
static void C_colmod_derivs (double *x, double *y,
                        double *ydot, double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, X, ans;

                                REAL(EPS)[0] = *eps;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                         incr_N_Protect();
  PROTECT(R_fcall = lang4(R_cont_deriv_func,X,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));               incr_N_Protect();

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(3);
}


/* jacobian                                                                   */
static void C_colmod_jac (double *x, double *y, double *pd, int *neq,
                        double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, X, ans;

                              REAL(EPS)[0] = *eps;
  for (i = 0; i < mstar; i++) REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                       incr_N_Protect();
  PROTECT(R_fcall = lang4(R_cont_jac_func,X,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));             incr_N_Protect();

  for (i = 0; i < *neq * mstar; i++)  pd[i] = REAL(ans)[i];
  my_unprotect(3);
}

/*  boundary condition                                                        */
static void C_colmod_bound (int *ii, double *y, double *gout, double *eps,
                          double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, J, ans;
  
                             REAL(EPS)[0]  = *eps;

  for (i = 0; i < mstar ; i++)  REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                     incr_N_Protect();
  PROTECT(R_fcall = lang4(R_cont_bound_func,J,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));               incr_N_Protect();
  /* only on e element returned... */
  gout[0] = REAL(ans)[0];
  my_unprotect(3);
}

/* jacobian of boundary condition                                             */
static void C_colmod_jacbound (int *ii, double *y, double *dg, double *eps,
                             double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, J, ans;
  
                             REAL(EPS)[0]  = *eps;

  for (i = 0; i < mstar; i++) REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                        incr_N_Protect();
  PROTECT(R_fcall = lang4(R_cont_jacbound_func,J,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));                  incr_N_Protect();

  for (i = 0; i < mstar ; i++)  dg[i] = REAL(ans)[i];
  my_unprotect(3);
}

/* -----------------------------------------------------------------------------
  give name to data types
----------------------------------------------------------------------------- */


/* -----------------------------------------------------------------------------
                  MAIN C-FUNCTION, CALLED FROM R-code
----------------------------------------------------------------------------- */

SEXP call_colmodsys(SEXP Ncomp, SEXP Mstar, SEXP M, SEXP Xout, SEXP Aleft,
    SEXP Aright, SEXP Zeta, SEXP Iset, SEXP Ltol, SEXP Tol, SEXP Fixpnt,
		SEXP Rpar, SEXP Ipar, SEXP Epsini, SEXP Eps, SEXP derivfunc,
		SEXP jacfunc, SEXP boundfunc, SEXP jacboundfunc, SEXP guessfunc,
    SEXP Initfunc, SEXP Parms, SEXP flist,
    SEXP Rwork, SEXP Iwork, SEXP Absent, SEXP RRwork, SEXP rho)
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, ICOUNT, RWORK, EPSS;

  int  j, k, ii, nx, ncomp, isForcing, isDll;
  double *aleft, *aright, *zeta, *fspace, *tol, *fixpnt, *z;
  double epsini, epsmin, xout, *rpar;
  int *m, *ispace, *ipar, *iset, *ltol, iflag, *icount, FullOut;
  int *absent;
  double *rwork;

/* pointers to functions passed to FORTRAN                                    */
  C_deriv_func2_type    *deriv_func = NULL;
  C_jac_func2_type      *jac_func = NULL;
  C_jacbound_func2_type *jacbound_func = NULL;
  C_bound_func2_type    *bound_func = NULL;
  C_guess_func2_type    *guess_func = NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();

  ncomp = INTEGER(Ncomp)[0];   /* number of equations -global variable */
  n_eq  = INTEGER(Ncomp)[0];   /* number of equations -global variable */
  mstar = INTEGER(Mstar)[0];   /* number of variables */

  m  = (int *) R_alloc(n_eq, sizeof(int));  /* order of diff eqns */
  for (j = 0; j < n_eq; j++) m[j] = INTEGER(M)[j];

  ii = LENGTH(Aleft);
  aleft  =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) aleft[j] = REAL(Aleft)[j];

  ii = LENGTH(Aright);
  aright =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) aright[j] = REAL(Aright)[j];

  ii = LENGTH(Zeta);
  zeta   =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) zeta[j] = REAL(Zeta)[j];

  ii = LENGTH(Iset) - 1;    /* length of iset */
  iset  = (int *)    R_alloc(ii, sizeof(int));
    for (j = 0; j < ii;j++) iset[j] = INTEGER(Iset)[j];
  FullOut = INTEGER(Iset)[ii];

  ii = LENGTH(Tol);
  tol   =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) tol[j] = REAL(Tol)[j];
    
  ltol = (int *)    R_alloc(ii, sizeof(int));
    for (j = 0; j < ii;j++) ltol[j] = INTEGER(Ltol)[j];

/*      error("mstar, ncomp,ltol %i %i %i %f",ncomp,mstar,ltol,tol[0]);*/

  ii =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) fixpnt[j] = REAL(Fixpnt)[j];

  ii = iset[5];
  ispace = (int *) R_alloc(ii, sizeof(int));

  ii = iset[4];
  fspace = (double *) R_alloc(ii, sizeof(double));

  if (iset[8] > 1) {   /* continuation request: values in Rwork, Iwork*/
    ii = LENGTH(Rwork);
    for (j=0; j < ii; j++) fspace[j] = REAL(Rwork)[j];
    ii = LENGTH(Iwork);
    for (j=0; j < ii; j++) ispace[j] = INTEGER(Iwork)[j];
  }

  ii = LENGTH(Absent);
  absent = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) absent[j] = INTEGER(Absent)[j];

  ii = LENGTH(RRwork);
  rwork = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rwork[j] = REAL(RRwork)[j];

  epsmin = REAL(Eps)[0];
  epsini = REAL(Epsini)[0];

  ii =  LENGTH(Rpar);
  rpar = (double *) R_alloc(ii, sizeof(double));
    for (j = 0; j < ii;j++) rpar[j] = REAL(Rpar)[j];
  ii =  LENGTH(Ipar);
  ipar = (int *) R_alloc(ii, sizeof(int));
    for (j = 0; j < ii;j++) ipar[j] = INTEGER(Ipar)[j];
  icount = (int *)  R_alloc(7, sizeof(int));

/* initialise global R-variables...                                           */
  /* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  if (isDll == 0) {
    PROTECT(EPS = NEW_NUMERIC(1));              incr_N_Protect();
    PROTECT(Y = allocVector(REALSXP,mstar));    incr_N_Protect();
  }
  /* Initialization of Parameters and Forcings (DLL functions)  */
  isForcing = initForcings(flist);
  epsval = (double *) R_alloc(1, sizeof(double)); epsval[0] = 0.;

  initParms(Initfunc, Parms);

  R_envir = rho;

  /* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }
  /* pointers to functions passed to FORTRAN */
  if (isDll) {   /* DLL addresses passed to fortran */
     deriv_func    = (C_deriv_func2_type *)    dll_colmod_deriv_func;
     derfun        = (C_deriv_func_type *)     R_ExternalPtrAddrFn_(derivfunc);

     if (absent[0] == 0) { 
       jac_func      = (C_jac_func2_type *)      dll_colmod_jac_func;
       jacfun        = (C_jac_func_type *)       R_ExternalPtrAddrFn_(jacfunc);
     }
     if (absent[1] == 0) { 
       bound_func    = (C_bound_func2_type *)    dll_colmod_bound_func;
       boundfun      = (C_bound_func_type *)     R_ExternalPtrAddrFn_(boundfunc);
     }
     if (absent[2] == 0) { 
       jacbound_func = (C_jacbound_func2_type *) dll_colmod_jacbound_func;
       jacboundfun   = (C_jacbound_func_type *)  R_ExternalPtrAddrFn_(jacboundfunc);
     }

	  /* here overruling deriv_func if forcing */
      if (isForcing) {
        deriv_func = (C_deriv_func2_type *) dll_colmod_deriv_func_forc;
      }
  } else {      /* interface functions between fortran and R */

      deriv_func = (C_deriv_func2_type *) C_colmod_derivs;
      R_cont_deriv_func = derivfunc;

      if (absent[0] == 0) {
        jac_func = C_colmod_jac;
        R_cont_jac_func = jacfunc;
      }

      if (absent[1] == 0) {
        bound_func = C_colmod_bound;
        R_cont_bound_func = boundfunc;
      }
      
      if (absent[2] == 0) {
        jacbound_func = C_colmod_jacbound;
        R_cont_jacbound_func = jacboundfunc;
      }
   }

    if (absent[0] == 1) {
        dy     = (double *) R_alloc(ncomp, sizeof(double));
        dycopy = (double *) R_alloc(ncomp, sizeof(double));

        jac_func = (C_jac_func2_type *) C_num_epsjac_func;
        jepsderfun  = deriv_func;
    }
    
    if (absent[1] == 1) {
      bound_func = (C_bound_func2_type *) C_num_epsbound_func;
      iibb = (int *) R_alloc(mstar, sizeof(int));
      for (j = 0; j < mstar; j++)
        iibb[j] = absent[3 + j];
      bb = (double *) R_alloc(mstar, sizeof(double));
      for (j = 0; j < mstar; j++)
        bb[j] = rwork[j];
    }
    
    if (absent[2] == 1) {
      jacbound_func = (C_jacbound_func2_type *) C_num_epsjacbound_func;
      jepsbndfun = bound_func;
      g = (double *) R_alloc(1, sizeof(double));
      gcopy = (double *) R_alloc(1, sizeof(double));

    }
   ycopy  = (double *) R_alloc(mstar, sizeof(double)); 
   guess_func = (C_guess_func2_type *) C_colmod_guess;
   R_cont_guess_func = guessfunc;

/* Call the fortran function                                                  */

	  F77_CALL(colmod) (&ncomp, m, aleft, aright, zeta, iset, ltol,
        tol, fixpnt, ispace, fspace, &iflag, &epsini, &epsmin,
        deriv_func,jac_func,bound_func,jacbound_func,
        guess_func, rpar, ipar, icount);

/*             Call Appsln(Xx,Z,Fspace,Ispace)
C....   Iflag - The Mode Of Return From Colmod.
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
	  unprotect_all();
	  error("The collocation matrix is singular for the final continuation problem\n");
	}
  else if (iflag == -1)
	{
	  unprotect_all();
	  error("The Expected No. Of Subintervals Exceeds Storage Specifications.\n");
	}
  else if (iflag == -2)
	{
	  unprotect_all();
	  error("The Nonlinear Iteration Has Not Converged For The Final Continuation Problem.\n");
	}
  else  if (iflag == -3)
	{
	  unprotect_all();
	  error("Illegal input to colmod\n");
	}
  else
	{
    if (iflag == 2) {
 	    Rprintf("The Problem was **not** solved for the requested eps value\n");
      Rprintf("Results are for eps equal to %g.\n", epsmin);
	  }
    if (iflag == -4)
	   warning("The maximum number of continuation steps has been reached and the model did not yet converge - returning the last values \n");

    nx = LENGTH(Xout);
    z  =(double *) R_alloc(mstar, sizeof(double));

    PROTECT(yout = allocMatrix(REALSXP,mstar+1,nx));incr_N_Protect();
	  for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(mappsln)(&xout,z,fspace,ispace);
                 for (j=0;j<mstar;j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
      }  /* end main x loop */
  ii = ncomp+7;
  PROTECT(ICOUNT = allocVector(INTSXP, 7)); incr_N_Protect();
  PROTECT(ISTATE = allocVector(INTSXP, ii+6)); incr_N_Protect();
  INTEGER(ISTATE)[0] = iflag;
  for (k = 0; k < 7; k++) INTEGER(ICOUNT)[k] = icount[k];
  for (k = 0; k < 5; k++)  INTEGER(ISTATE)[k+1] = icount[k];
  for (k = 0; k < ii; k++) INTEGER(ISTATE)[k+6] = ispace[k];
  if (FullOut)
    ii = ispace[6];
  else 
    ii = 1;
  PROTECT(EPSS = allocVector(REALSXP, 2)); incr_N_Protect();
  REAL(EPSS)[0] = epsini;
  REAL(EPSS)[1] = epsmin;
  
  setAttrib(yout, install("eps"), EPSS);
  PROTECT(RWORK = allocVector(REALSXP, ii));incr_N_Protect();
  for (k = 0; k < ii; k++) REAL(RWORK)[k] = fspace[k];
  setAttrib(yout, install("istate"), ISTATE);
  setAttrib(yout, install("icount"), ICOUNT);
  setAttrib(yout, install("rstate"), RWORK);
  }
/*                       ####   termination   ####                            */    
  unprotect_all();
  return(yout);
}

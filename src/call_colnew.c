#include <time.h>
#include <string.h>
#include "bvpSolve.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   boundary value problem solvers.
   
   The C-wrappers that provide the interface between FORTRAN codes and R-code 
   are: C_bvp_deriv_func   : interface with R-code "derivfunc", passes derivatives  
        C_bvp_jac_func     : interface with R-code "jacfunc", passes jac_funcobian
        C_bvp_bound_func   : interface with R-code "bound_func", boundaries
        C_bvp_jacbound_func: interface with R-code "jacbound_func", jac_funcobian  of boundaries
        C_bvp_guess_func   : interface with R-code "guess_func", initial estimates
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
         void (*)(double *, double *, double *),                           /* guess_func */
         double *, int *, int*) ;

void F77_NAME(colsys)(int*, int*, double *, double *, double *, int *, int *,
         double *, double *, int *, double *, int *, 
         void (*)(int *, double *, double *, double *, double *, int *),   /* fsub  */
		     void (*)(int *, double *, double *, double *, double *, int *),   /* dfsub */
			   void (*)(int *, int *, double *, double *, double *, int *),      /* gsub  */
		     void (*)(int *, int *, double *, double *, double *, int *),      /* dgsub */
         void (*)(double *, double *, double *),                           /* guess_func */
         double *, int *, int*) ;

void F77_NAME(appsln)(double *, double *, double *, int *);
void F77_NAME(sysappsln)(double *, double *, double *, int *);


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
  SEXP R_fcall, ans, R_fcall2, ans2;

  REAL(X)[0]   = *x;

  PROTECT(R_fcall = lang2(R_bvp_guess_func, X));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));            incr_N_Protect();

  p = fmax(1e-7, *x*1e-7 );
  REAL(X)[0]   = *x+p;
  PROTECT(R_fcall2 = lang2(R_bvp_guess_func, X));   incr_N_Protect();
  PROTECT(ans2 = eval(R_fcall2, R_envir));          incr_N_Protect();

  /* both have the same dimensions... */
  for (i = 0; i < n_eq; i++) y[i] = REAL(ans)[i];
  for (i = 0; i < n_eq; i++) ydot[i] = (REAL(ans2)[i]-y[i])/p;

  my_unprotect(4);

}

/* derivative function                                                        */
static void C_bvp_deriv_func (int * n, double *x, double *y,
                              double *ydot, double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                                REAL(X)[0]   = *x;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang3(R_bvp_deriv_func,X,Y)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));     incr_N_Protect();

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);
}

/* jacobian                                                                   */
static void C_bvp_jac_func (int *n, double *x, double *y, double *pd,
                            double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                              REAL(X)[0]   = *x;
  for (i = 0; i < mstar; i++) REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang3(R_bvp_jac_func,X,Y)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));   incr_N_Protect();

  for (i = 0; i < n_eq * mstar; i++)  pd[i] = REAL(ans)[i];
  my_unprotect(2);
}

/*  boundary condition                                                        */
static void C_bvp_bound_func (int *ii, int * n, double *y, double *gout,
                              double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                                INTEGER(J)[0] = *ii;
  for (i = 0; i < mstar ; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(R_bvp_bound_func,J,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));            incr_N_Protect();
  /* only one element returned... */
  gout[0] = REAL(ans)[0];

  my_unprotect(2);
}

/* jacobian of boundary condition                                             */
static void C_bvp_jacbound_func (int *ii, int *n, double *y, double *dg,
                                 double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < mstar; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(R_bvp_jacbound_func,J,Y)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));             incr_N_Protect();

  for (i = 0; i < mstar ; i++)  dg[i] = REAL(ans)[i];
  my_unprotect(2);
}

/* -----------------------------------------------------------------------------
                  MAIN C-FUNCTION, CALLED FROM R-code
----------------------------------------------------------------------------- */

SEXP call_colnew(SEXP Ncomp, SEXP Xout, SEXP Aleft, SEXP Aright,
		SEXP Zeta, SEXP Mstar, SEXP M, SEXP Iset, SEXP Rwork, SEXP Iwork,
    SEXP Tol, SEXP Fixpnt, SEXP Rpar, SEXP Ipar,
		SEXP derivfunc, SEXP jacfunc, SEXP boundfunc,
    SEXP jacboundfunc, SEXP guessfunc, SEXP Initfunc, SEXP Parms, SEXP flist,
    SEXP Type, SEXP rho)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, ICOUNT, RWORK;

  int  j, ii, ncomp, k, nx, ntol, nfixpnt, isForcing, type;
  double aleft, aright, *zeta, *fspace, *tol, *fixpnt, *z, *rpar;
  double xout;
  int *m, *ispace, *iset, *icount, *ltol, *ipar, iflag, isDll, FullOut;

  C_deriv_func_type    *deriv_func;
  C_jac_func_type      *jac_func;
  C_bound_func_type    *bound_func;
  C_jacbound_func_type *jacbound_func;
  C_guess_func_type    *guess_func;
  
/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();

  aleft  =REAL(Aleft)[0];
  aright =REAL(Aright)[0];

  ncomp = INTEGER(Ncomp)[0];     /* number of equations -global variable */
  type  = INTEGER(Type)[0];      /* 2 = bvpcol */
  
  n_eq  = INTEGER(Ncomp)[0];     /* number of equations -global variable */
  mstar = INTEGER(Mstar)[0];     /* number of variables */

  /* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  m  = (int *) R_alloc(n_eq, sizeof(int));  /* order of diff eqns */
  for (j = 0; j < n_eq; j++) m[j] = INTEGER(M)[j];


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

  PROTECT(X  = NEW_NUMERIC(1));                 incr_N_Protect();
  if (isDll == 0) {
    PROTECT(J  = NEW_INTEGER(1));               incr_N_Protect();
    PROTECT(Y = allocVector(REALSXP,mstar));    incr_N_Protect();
  }
  /* Initialization of Parameters and Forcings (DLL functions)  */
  isForcing = initForcings(flist);
  initParms(Initfunc, Parms);

  R_envir = rho;

  /* pointers to functions passed to FORTRAN */
  if (isDll) {   /* DLL addresses passed to fortran */
      deriv_func    = (C_deriv_func_type *)    R_ExternalPtrAddr(derivfunc);
      jac_func      = (C_jac_func_type *)      R_ExternalPtrAddr(jacfunc);
      bound_func    = (C_bound_func_type *)    R_ExternalPtrAddr(boundfunc);
      jacbound_func = (C_jacbound_func_type *) R_ExternalPtrAddr(jacboundfunc);

	  /* here overruling deriv_func if forcing */
      if (isForcing) {
        derfun =     (C_deriv_func_type *) R_ExternalPtrAddr(derivfunc);
        deriv_func = (C_deriv_func_type *) dll_bvp_deriv_func_forc;
      }

  } else {      /* interface functions between fortran and R */
      deriv_func = C_bvp_deriv_func;
      R_bvp_deriv_func = derivfunc;

      jac_func = C_bvp_jac_func;
      R_bvp_jac_func = jacfunc;

      bound_func = C_bvp_bound_func;
      R_bvp_bound_func = boundfunc;

      jacbound_func = C_bvp_jacbound_func;
      R_bvp_jacbound_func = jacboundfunc;
    }

  guess_func = (C_guess_func_type *) C_bvp_guess_func;
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
   else
	  F77_CALL(colsys) (&ncomp, m, &aleft, &aright, zeta, iset, ltol,
        tol, fixpnt, ispace, fspace, &iflag, 
        deriv_func, jac_func, bound_func, jacbound_func, 
        guess_func, rpar, ipar, icount);

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
	  error("Illegal input to bvpcol\n");
	}
  else
	{
    nx = LENGTH(Xout);
    z  =(double *) R_alloc(mstar, sizeof(double));

    PROTECT(yout = allocMatrix(REALSXP,mstar+1,nx));incr_N_Protect();
   if (type == 0) 
	  for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(appsln)(&xout,z,fspace,ispace);
                 for (j=0;j<mstar;j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
      }  /* end main x loop */
    else
	   for (k = 0; k < nx; k++)
      {          xout = REAL(Xout)[k];
                 REAL(yout)[k*(mstar+1)] = xout;
                 F77_CALL(sysappsln)(&xout,z,fspace,ispace);
                 for (j=0;j<mstar;j++) REAL(yout)[k*(mstar+1) + j+1] = z[ j];
      }  /* end main x loop */

  ii = ncomp+7;
  PROTECT(ICOUNT = allocVector(INTSXP, 6));incr_N_Protect();
  PROTECT(ISTATE = allocVector(INTSXP, ii+6));incr_N_Protect();
  INTEGER(ISTATE)[0] = iflag;
  for (k = 0; k < 6; k++)  INTEGER(ICOUNT)[k] = icount[k];
  for (k = 0; k < 5; k++)  INTEGER(ISTATE)[k+1] = icount[k];
  for (k = 0; k < ii; k++) INTEGER(ISTATE)[k+6] = ispace[k];
  if (FullOut) 
    ii = ispace[6]; 
  else 
    ii = 1;

  PROTECT(RWORK = allocVector(REALSXP, ii));incr_N_Protect();
  for (k = 0; k<ii; k++) REAL(RWORK)[k] = fspace[k];
  setAttrib(yout, install("icount"), ICOUNT);
  setAttrib(yout, install("istate"), ISTATE);
  setAttrib(yout, install("rstate"), RWORK); 
  }
/*                       ####   termination   ####                            */    
  unprotect_all();
  return(yout);
}

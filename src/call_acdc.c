#include <time.h>
#include <string.h>
#include "bvpSolve.h"

/* -----------------------------------------------------------------------------
  definition of the calls to the fortran functions
      Subroutine acdc(Ncomp, Nlbc, Nucol, Aleft, Aright, Nfxpnt, Fixpnt,
     +            Ntol, Ltol, Tol, Linear, Givmsh, Giveu,
     +            Full,nmshguess, xguess, nugdim, uguess,Nmsh, Xx,
     +            Nudim, U, Nmax, Lwrkfl, Wrk, Lwrkin, Iwrk, Giveps,
     +            Eps, Epsmin, acfsub, acdfsub, acgsub, acdgsub,
     +            ckappa1,gamma1,sigma,ckappa,ckappa2,rpar,ipar,icount,
     +            precis, useC, Iflbvp)                                              */
void F77_NAME(acdc)(int*, int*, int*, double *, double *,  int *, double *,
 int *, int *,double *, int *, int *, int *, int *, int *, double *, int *,
 double *, int *, double *, int *, double *, int *, int *, double *, int *,
 int *, int *, double *, double *,
 void (*)(int *, double *, double *, double *, double *, double *, int *), /* fsub(n,x,u,f,eps,rp,ip)   */
 void (*)(int *, double *, double *, double *, double *, double *, int *), /* dfsub(n,x,u,f,eps,rp,ip)   */
 void (*)(int *, int *, double *, double *, double *, double *, int *),    /* gsub(i,n,u,g,eps,rp,ip)   */
 void (*)(int *, int *, double *, double *, double *, double *, int *),    /* dgsub(i,n,u,dg,eps,rp,ip) */
 double *, double *, double *, double *, double *, double *, int *,
          int *, double *, int *, int*);

/* -----------------------------------------------------------------------------
                        - when model in compiled code
----------------------------------------------------------------------------- */

/* wrapper above the derivate function that first estimates the
values of the forcing functions and puts new value of eps in par              */

static void dll_bvp_deriv_func_forc (int *n, double *x, double *y,
                         double *ydot, double *eps, double *rpar, int *ipar)
{
  updatedeforc(x);
  epsval[0] = eps[0];
  derfun(n, x, y, ydot, rpar, ipar);
}

/* wrapper above the functions that puts new value of eps in par
The interface is slightly different in colmod - overruled in C-code           */

static void dll_bvp_deriv_func (int *n, double *x, double *y,
                         double *ydot, double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  derfun(n, x, y, ydot, rpar, ipar);
}

static void dll_bvp_jac_func (int *n, double *x, double *y,
                         double *pd, double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  jacfun(n, x, y, pd, rpar, ipar);
}

static void dll_bvp_bound_func (int *ii, int *n, double *y, double *gout,
                        double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  boundfun(ii, n, y, gout, rpar, ipar);
}
static void dll_bvp_jacbound_func (int *ii, int *n, double *y, double *dg,
                        double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  jacboundfun(ii, n, y, dg, rpar, ipar);
}

/* -----------------------------------------------------------------------------
   interface between fortran function calls and R functions
   Note: passing of parameter values and "..." is done in R-function bvptwp
----------------------------------------------------------------------------- */

static void C_acdc_deriv_func (int *n, double *x, double *y,
                        double *ydot, double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, ans;
                                REAL(EPS)[0] = *eps;
                                REAL(X)[0]   = *x;
  for (i = 0; i < n_eq ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang4(R_cont_deriv_func,X,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));               incr_N_Protect();

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);
}

/* interface between fortran call to jacobian and R function                  */
static void C_acdc_jac_func (int *n, double *x, double *y, double *pd,
                        double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, ans;
                              REAL(EPS)[0] = *eps;
                              REAL(X)[0]   = *x;
  for (i = 0; i < n_eq; i++) REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang4(R_cont_jac_func,X,Y,EPS));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));               incr_N_Protect();

  for (i = 0; i < n_eq * n_eq; i++)  pd[i] = REAL(ans)[i];
  my_unprotect(2);
}

/* interface between fortran call to boundary condition and R function        */

static void C_acdc_bound_func (int *ii, int *n, double *y, double *gout,
                        double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, ans;
                             REAL(EPS)[0]  = *eps;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < n_eq ; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang4(R_cont_bound_func,J,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));               incr_N_Protect();
  /* only one element returned... */
  gout[0] = REAL(ans)[0];
  my_unprotect(2);
}
/*interface between fortran call to jacobian of boundary and R function      */

static void C_acdc_jacbound_func (int *ii, int *n, double *y, double *dg,
                           double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, ans;
                             REAL(EPS)[0]  = *eps;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < n_eq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang4(R_cont_jacbound_func,J,Y,EPS)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));                  incr_N_Protect();

  for (i = 0; i < n_eq ; i++)  dg[i] = REAL(ans)[i];
  my_unprotect(2);
}

/* -----------------------------------------------------------------------------
  give name to data types
----------------------------------------------------------------------------- */

typedef void C_acdc_deriv_func_type(int *, double *, double *,double *,
                                    double *, double *, int *);
typedef void C_acdc_bound_func_type(int *, int *, double *, double *,
                                    double *, double *, int *);
typedef void C_acdc_jac_func_type  (int *,  double *, double *, double *,
                                    double *, double *, int *);
typedef void C_acdc_jacbound_func_type(int *, int *, double *, double *,
                                    double *, double *, int *);

/* -----------------------------------------------------------------------------
                  MAIN C-FUNCTION, CALLED FROM R-code
----------------------------------------------------------------------------- */

SEXP call_acdc(SEXP Ncomp, SEXP Fixpnt, SEXP Aleft, SEXP Aright,
		SEXP Nlbc, SEXP Tol, SEXP Linear, SEXP Full, SEXP Givmesh, SEXP Givu,
    SEXP Nmesh, SEXP Nmax, SEXP Lwrkfl, SEXP Lwrkin, SEXP Xguess, SEXP Yguess,
    SEXP Rpar, SEXP Ipar, SEXP UseC, SEXP Epsini, SEXP Eps,
    SEXP derivfunc, SEXP jacfunc, SEXP boundfunc, SEXP jacboundfunc,
    SEXP Initfunc, SEXP Parms, SEXP flist, SEXP rho)
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, RSTATE, EPSS;

  int  j, ii, ncomp, nlbc, nmax, lwrkfl, lwrkin, nx, *ipar, isForcing;
  double *wrk, *tol, *fixpnt, *u, *xx, *rpar, *precis, *xguess, *yguess;
  double epsmin, epsini, aleft, aright, ckappa1, gamma1, sigma, ckappa, ckappa2;
  int *icount, *ltol, *iwrk, ntol, iflag, nfixpnt, linear, givmesh;
  int full, useC, givu, giveps, nmesh, isDll, nugdim, nmshguess;

/* pointers to functions passed to FORTRAN                                    */
  C_acdc_deriv_func_type    *deriv_func;
  C_acdc_jac_func_type      *jac_func;
  C_acdc_bound_func_type    *bound_func;
  C_acdc_jacbound_func_type *jacbound_func;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */
  init_N_Protect();

  aleft  =REAL(Aleft)[0];
  aright =REAL(Aright)[0];

  ncomp  = INTEGER(Ncomp)[0];    /* number of equations */
  n_eq   = ncomp;

  nlbc    = INTEGER(Nlbc)[0];    /* number of left boundary conditions */
  nmax    = INTEGER(Nmax)[0];    /* max number of mesh points */
  lwrkfl  = INTEGER(Lwrkfl)[0];  /* length of double workspace */
  lwrkin  = INTEGER(Lwrkin)[0];  /* length of integer workspace */
  linear  = INTEGER(Linear)[0];  /* true if linear problem */
  full    = INTEGER(Full)[0];    /* true if verbose output */
  givu    = INTEGER(Givu)[0];    /* true if initial trial solution given */
  givmesh = INTEGER(Givmesh)[0]; /* true if initial mesh given */
  nmesh   = INTEGER(Nmesh)[0];   /* size of mesh */
  useC    = INTEGER(UseC)[0];    /* conditioning or not */

  /* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  /* copies of variables that will be changed in the FORTRAN subroutine */
  ntol = LENGTH(Tol);
  tol   =(double *) R_alloc(ntol, sizeof(double));
  for (j = 0; j < ntol; j++) tol[j] = REAL(Tol)[j];

  ltol   =(int *) R_alloc(ntol, sizeof(int));
  for (j = 0; j < ntol; j++) ltol[j] = j+1;

  nfixpnt =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(nfixpnt, sizeof(double));
  for (j = 0; j < nfixpnt;j++) fixpnt[j] = REAL(Fixpnt)[j];

  xx   =(double *) R_alloc(nmax, sizeof(double));
  for (j = 0; j < nmesh; j++) xx[j] = REAL(Xguess)[j];
  for (j = nmesh; j < nmax; j++) xx[j] = 0;

  // to be used in fortran code initu
  if (givu) {
    xguess = (double *) R_alloc(nmesh, sizeof(double));
    for (j = 0; j < nmesh; j++) xguess[j] = REAL(Xguess)[j];

    yguess = (double *) R_alloc(nmesh*ncomp, sizeof(double));
    for (j = 0; j < nmesh*ncomp; j++) yguess[j] = REAL(Yguess)[j];

  } else { /* dummy variables */
    xguess = (double *) R_alloc(1, sizeof(double));
    xguess[0] = 0.;
    yguess = (double *) R_alloc(ncomp, sizeof(double));
    for (j = 0; j <  ncomp; j++) yguess[j] = 0.;

  }

  ii = nmax*ncomp;
  u   =(double *) R_alloc(ii, sizeof(double));
   for (j = 0; j < nmesh*ncomp; j++) u[j] = REAL(Yguess)[j];
   for (j = nmesh*ncomp; j < nmax*ncomp; j++) u[j] = 0;

  wrk = (double *) R_alloc(lwrkfl, sizeof(double));
     for (j = 0; j < lwrkfl; j++) wrk[j] = 0.;

  iwrk = (int *)   R_alloc(lwrkin, sizeof(int));
     for (j = 0; j < lwrkin; j++) iwrk[j] = 0;

  precis = (double *) R_alloc(3,sizeof(double));
  precis[0] = DBL_MIN;
  precis[1] = DBL_MAX;
  precis[2] = DBL_EPSILON/FLT_RADIX;

  epsmin = REAL(Eps)[0];
  epsini = REAL(Epsini)[0];
  giveps = 1;
  icount = (int *)  R_alloc(8, sizeof(int));

  ii = LENGTH(Ipar);
  ipar = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) ipar[j] = INTEGER(Ipar)[j];

  ii = LENGTH(Rpar);
  rpar = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rpar[j] = REAL(Rpar)[j];

  /* initialise global R-variables... */
  if (isDll == 0) {
    PROTECT(X  = NEW_NUMERIC(1));               incr_N_Protect();
    PROTECT(EPS = NEW_NUMERIC(1));              incr_N_Protect();
    PROTECT(J = NEW_INTEGER(1));                incr_N_Protect();
    PROTECT(Y = allocVector(REALSXP,ncomp));    incr_N_Protect();
  }

  /* Initialization of Parameters and Forcings (DLL functions)   */
  isForcing = initForcings(flist);
  initParms(Initfunc, Parms);

  R_envir = rho;

  /* pointers to functions passed to FORTRAN */
  if (isDll) {
      deriv_func    = (C_acdc_deriv_func_type *)    dll_bvp_deriv_func;
      jac_func      = (C_acdc_jac_func_type *)      dll_bvp_jac_func;
      bound_func    = (C_acdc_bound_func_type *)    dll_bvp_bound_func;
      jacbound_func = (C_acdc_jacbound_func_type *) dll_bvp_jacbound_func;

      derfun        = (C_deriv_func_type *)         R_ExternalPtrAddr(derivfunc);
      jacfun        = (C_jac_func_type *)           R_ExternalPtrAddr(jacfunc);
      boundfun      = (C_bound_func_type *)         R_ExternalPtrAddr(boundfunc);
      jacboundfun   = (C_jacbound_func_type *)      R_ExternalPtrAddr(jacboundfunc);

	  /* here overruling deriv_func if forcing  */
      if (isForcing) {
        deriv_func = (C_acdc_deriv_func_type *)     dll_bvp_deriv_func_forc;
      }
  } else {
      deriv_func = C_acdc_deriv_func;
      R_cont_deriv_func = derivfunc;

      jac_func = C_acdc_jac_func;
      R_cont_jac_func = jacfunc;

      bound_func = C_acdc_bound_func;
      R_cont_bound_func = boundfunc;

      jacbound_func = C_acdc_jacbound_func;
      R_cont_jacbound_func = jacboundfunc;
    }

/* Call the fortran function acdc               */
  nugdim = ncomp;
  nmshguess = nmesh;
/*  Rprintf("precis %g %g %g\n",precis[0],precis[1],precis[2]);*/
	F77_CALL(acdc) (&ncomp, &nlbc, &nmax, &aleft, &aright, &nfixpnt, fixpnt,
      &ntol, ltol, tol, &linear, &givmesh, &givu, &full, &nmshguess,
      xguess, &nugdim, yguess, &nmesh, xx, &ncomp, u, &nmax, &lwrkfl,
      wrk, &lwrkin, iwrk, &giveps, &epsini, &epsmin,
      deriv_func, jac_func, bound_func, jacbound_func,
      &ckappa1, &gamma1, &sigma, &ckappa, &ckappa2, rpar, ipar, icount,
      precis, &useC, &iflag);


/* error("Till here.\n"); iflag - The Mode Of Return From acdc                                      */
	if (iflag == 4)      {
	   unprotect_all();
     error("One of the input parameters is invalid.\n");
  } else if (iflag == 1) 	{
	  unprotect_all();
	  error("Terminated: final problem not solved.\n");
	} else if (iflag == 2) 	{
	  unprotect_all();
	  error("Terminated: too many continuation steps\n");
	} else if (iflag == 3)  	{
	  unprotect_all();
	  error("Terminated: ill conditioned problem.\n");
	} else	{
/*                   ####   returning output   ####                           */
    nx = nmesh;

    PROTECT(yout = allocVector(REALSXP,(ncomp+1)*(nx))); incr_N_Protect();
	  for (j = 0; j < nx; j++)       REAL(yout)[j]    = xx[j];
    for (j = 0; j < ncomp*nx; j++) REAL(yout)[nx+j] =  u[j];
  }

  PROTECT(ISTATE = allocVector(INTSXP, 13)); incr_N_Protect();
  INTEGER(ISTATE)[0] = iflag;
  for (j = 0; j < 7; j++)
    INTEGER(ISTATE)[1+j] = icount[j];
  INTEGER(ISTATE)[9] = nmax;
  INTEGER(ISTATE)[10] = nmesh;
  INTEGER(ISTATE)[11] = lwrkfl;
  INTEGER(ISTATE)[12] = lwrkin;

  setAttrib(yout, install("istate"), ISTATE);

  PROTECT(EPSS = allocVector(REALSXP, 2)); incr_N_Protect();
  REAL(EPSS)[0] = epsini;
  REAL(EPSS)[1] = epsmin;

  setAttrib(yout, install("eps"), EPSS);

  PROTECT(RSTATE = allocVector(REALSXP, 5)); incr_N_Protect();
  REAL(RSTATE)[0] = ckappa1;
  REAL(RSTATE)[1] = gamma1;
  REAL(RSTATE)[2] = sigma;
  REAL(RSTATE)[3] = ckappa;
  REAL(RSTATE)[4] = ckappa2;
  setAttrib(yout, install("rstate"), RSTATE);

/*               ####   termination   ####                                    */
  unprotect_all();
  return(yout);
}


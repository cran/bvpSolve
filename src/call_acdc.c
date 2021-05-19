#include <time.h>
#include <string.h>
#include "bvpSolve.h"
#include "externalptr.h"

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

static void C_num_acdcjac_func (int *n,  double *x, double *y, double *pd,
                            double *eps, double * rpar, int * ipar)
{
  int i, j;
  double perturb;

  for (i = 0; i < *n; i++) ycopy[i]   = y[i];

  jaderfun(n, x, y, dy, eps, rpar, ipar);
  for (i = 0; i < *n; i++) dycopy[i]   = dy[i];
  for (i = 0; i < *n * *n; i++) pd[i] = 0.;
  
  for (j = 0; j < *n; j++) {
     if (y[j] > 1.)
       perturb = y[j]*1e-8;
     else
       perturb = 1e-8 ;  
     ycopy[j] = y[j] + perturb;
     
     jaderfun(n, x, ycopy, dycopy, eps, rpar, ipar);
     
     ycopy[j] = y[j];
     
     for (i = 0; i < *n; i++) 
       pd[j* *n + i] = (dycopy[i] - dy[i])/perturb;

  }

}

static void C_num_acdcbound_func (int *ii, int *n, double *y, double *gout,
                              double *eps, double * rpar, int * ipar)
{
   int i, ib;
   i =  ii[0] - 1;        /*-1 to go from R to C indexing*/

   ib = iibb[i] - 1;

   gout[0] = y[ib] - bb[i];       
}

static void C_num_acdcjacbound_func (int *ii, int *n, double *y, double *dg,
                                 double *eps, double * rpar, int * ipar)
{
  int i;
  double perturb;
  double g, gcopy;
//  warning("entering numerical BOUNDARY jacobian %i %i %g %g\n", *ii, *n, y[0], dg[0]);

  for (i = 0; i < *n; i++) ycopy[i]  = y[i];
  for (i = 0; i < *n; i++) dg[i] = 0.;
  for (i = 0; i < *n; i++) {
    jabndfun(ii, n, y, &g, eps, rpar, ipar);

    if (y[i] > 1.)
       perturb = y[i]*1e-8;
    else
       perturb = 1e-8;  

    ycopy[i] = y[i] + perturb;
    jabndfun(ii, n, ycopy, &gcopy, eps, rpar, ipar);
    ycopy[i] = y[i];
    dg[i] = (gcopy - g)/perturb;
  }
}

/* wrapper above the derivate function that first estimates the
values of the forcing functions and puts new value of eps in rpar              */

static void dll_bvp_deriv_func_forc_eps (int *n, double *x, double *y,
                         double *ydot, double *eps, double *rpar, int *ipar)
{
  updatedeforc(x);
  epsval[0] = eps[0];   /* value of parameter */
  rpar[0] = eps[0];
  derfun(n, x, y, ydot, rpar, ipar);
}

/* wrapper above the functions that puts new value of eps in par
The interface is slightly different in colmod - overruled in C-code           */

static void dll_bvp_deriv_func (int *n, double *x, double *y,
                         double *ydot, double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  rpar[0] = eps[0];   /* value of parameter */
  derfun(n, x, y, ydot, rpar, ipar);
}

static void dll_bvp_jac_func (int *n, double *x, double *y,
                         double *pd, double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  rpar[0] = eps[0];   /* value of parameter */
  jacfun(n, x, y, pd, rpar, ipar);
}

static void dll_bvp_bound_func (int *ii, int *n, double *y, double *gout,
                        double *eps, double *rpar, int *ipar)
{
  epsval[0] = eps[0];   /* value of parameter */
  rpar[0] = eps[0];   /* value of parameter */
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
  SEXP R_fcall, X, ans;
                               REAL(EPS)[0] = *eps;
  for (i = 0; i < n_eq ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                         
  PROTECT(R_fcall = lang4(R_cont_deriv_func,X,Y,EPS)); 
  PROTECT(ans = eval(R_fcall, R_envir));               

  for (i = 0; i < n_eq ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  UNPROTECT(3);
}

/* interface between fortran call to jacobian and R function                  */
static void C_acdc_jac_func (int *n, double *x, double *y, double *pd,
                        double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, X, ans;
                             REAL(EPS)[0] = *eps;
  for (i = 0; i < n_eq; i++) REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                         
  PROTECT(R_fcall = lang4(R_cont_jac_func,X,Y,EPS));   
  PROTECT(ans = eval(R_fcall, R_envir));               

  for (i = 0; i < n_eq * n_eq; i++)  pd[i] = REAL(ans)[i];
  UNPROTECT(3);
}

/* interface between fortran call to boundary condition and R function        */

static void C_acdc_bound_func (int *ii, int *n, double *y, double *gout,
                        double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, J, ans;
                             REAL(EPS)[0]  = *eps;
  for (i = 0; i < n_eq ; i++)  REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                     
  PROTECT(R_fcall = lang4(R_cont_bound_func,J,Y,EPS)); 
  PROTECT(ans = eval(R_fcall, R_envir));               
  /* only one element returned... */
  gout[0] = REAL(ans)[0];
  UNPROTECT(3);
}
/*interface between fortran call to jacobian of boundary and R function      */

static void C_acdc_jacbound_func (int *ii, int *n, double *y, double *dg,
                           double *eps, double *rpar, int *ipar)
{
  int i;
  SEXP R_fcall, J, ans;
                             REAL(EPS)[0]  = *eps;

  for (i = 0; i < n_eq; i++) REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                        
  PROTECT(R_fcall = lang4(R_cont_jacbound_func,J,Y,EPS)); 
  PROTECT(ans = eval(R_fcall, R_envir));                  

  for (i = 0; i < n_eq ; i++)  dg[i] = REAL(ans)[i];
  UNPROTECT(3);
}

/* -----------------------------------------------------------------------------
  give name to data types
----------------------------------------------------------------------------- */


/* -----------------------------------------------------------------------------
                  MAIN C-FUNCTION, CALLED FROM R-code
----------------------------------------------------------------------------- */

SEXP call_acdc(SEXP Ncomp, SEXP Fixpnt, SEXP Aleft, SEXP Aright,
		SEXP Nlbc, SEXP Tol, SEXP Linear, SEXP Full, SEXP Givmesh, SEXP Givu,
    SEXP Nmesh, SEXP Nmax, SEXP Lwrkfl, SEXP Lwrkin, SEXP Xguess, SEXP Yguess,
    SEXP Rpar, SEXP Ipar, SEXP UseC, SEXP Epsini, SEXP Eps,
    SEXP derivfunc, SEXP jacfunc, SEXP boundfunc, SEXP jacboundfunc,
    SEXP Initfunc, SEXP Parms, SEXP flist, SEXP Absent, SEXP Rwork, SEXP rho)
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
  int *absent;
  double *rwork;

/* pointers to functions passed to FORTRAN                                    */
  C_acdc_deriv_func_type    *deriv_func;
  C_acdc_jac_func_type      *jac_func;
  C_acdc_bound_func_type    *bound_func;
  C_acdc_jacbound_func_type *jacbound_func = NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */
  int nprot = 0;

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

  ii = LENGTH(Absent);
  absent = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) absent[j] = INTEGER(Absent)[j];

  ii = LENGTH(Rwork);
  rwork = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rwork[j] = REAL(Rwork)[j];

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
    PROTECT(EPS = NEW_NUMERIC(1));              nprot++;
    PROTECT(Y = allocVector(REALSXP,ncomp));    nprot++;
  }

  /* Initialization of Parameters and Forcings (DLL functions)   */
  epsval = (double *) R_alloc(1, sizeof(double)); epsval[0] = 0.;
  isForcing = initForcings(flist);
  // initParms(Initfunc, Parms);
  if (Initfunc != NA_STRING) {
    if (inherits(Initfunc, "NativeSymbol"))  {
      init_func *initializer;
      
      PROTECT(bvp_gparms = Parms);     nprot++;
      initializer = (init_func *) R_ExternalPtrAddrFn_(Initfunc);
      initializer(Initbvpparms);
    }
  }

  R_envir = rho;

  /* pointers to functions passed to FORTRAN */
  if (isDll) {
      deriv_func      = (C_acdc_deriv_func_type *)    dll_bvp_deriv_func;
      if (absent[0] == 0)  
        jac_func      = (C_acdc_jac_func_type *)      dll_bvp_jac_func;

      if (absent[1] == 0)
        bound_func    = (C_acdc_bound_func_type *)    dll_bvp_bound_func;

      if (absent[2] == 0)
        jacbound_func = (C_acdc_jacbound_func_type *) dll_bvp_jacbound_func;

      derfun        = (C_deriv_func_type *)         R_ExternalPtrAddrFn_(derivfunc);
      jacfun        = (C_jac_func_type *)           R_ExternalPtrAddrFn_(jacfunc);
      boundfun      = (C_bound_func_type *)         R_ExternalPtrAddrFn_(boundfunc);
      jacboundfun   = (C_jacbound_func_type *)      R_ExternalPtrAddrFn_(jacboundfunc);

	  /* here overruling deriv_func if forcing  */
      if (isForcing) {
        deriv_func = (C_acdc_deriv_func_type *) dll_bvp_deriv_func_forc_eps;
      }

  } else {
      deriv_func = (C_acdc_deriv_func_type *)    C_acdc_deriv_func;
      R_cont_deriv_func = derivfunc;

      if (absent[0] == 0) {
        jac_func = C_acdc_jac_func;
        R_cont_jac_func = jacfunc;
      }

      if (absent[1] == 0) {
        bound_func = C_acdc_bound_func;
        R_cont_bound_func = boundfunc;
      }
      
      if (absent[2] == 0) {
        jacbound_func = C_acdc_jacbound_func;
        R_cont_jacbound_func = jacboundfunc;
      } 
    }

/* if numerical approximates should be used */    
    if (absent[0] == 1) {
        dy     = (double *) R_alloc(ncomp, sizeof(double));
        dycopy = (double *) R_alloc(ncomp, sizeof(double));
        ycopy  = (double *) R_alloc(ncomp, sizeof(double)); 
        jac_func = (C_acdc_jac_func_type *)      C_num_acdcjac_func;
        jaderfun  = (C_acdc_deriv_func_type *)    deriv_func;
    }
    if (absent[1] == 1) {
      bound_func = (C_acdc_bound_func_type *) C_num_acdcbound_func;
      iibb = (int *) R_alloc(ncomp, sizeof(int));
      for (j = 0; j < ncomp; j++)
        iibb[j] = absent[3 + j];
      bb = (double *) R_alloc(ncomp, sizeof(double));
      for (j = 0; j < ncomp; j++)
        bb[j] = rwork[j];
    }
    if (absent[2] == 1) {
        jacbound_func = (C_acdc_jacbound_func_type *) C_num_acdcjacbound_func;
        jabndfun = (C_acdc_bound_func_type *) bound_func;
        if (absent[0] != 1)
          ycopy  = (double *) R_alloc(ncomp, sizeof(double)); 
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
	   UNPROTECT(nprot);
     error("One of the input parameters is invalid.\n");
  } else if (iflag == 1) 	{
    UNPROTECT(nprot);
    error("Terminated: final problem not solved.\n");
	} else if (iflag == 2) 	{
	  UNPROTECT(nprot);
	  error("Terminated: too many continuation steps\n");
	} else if (iflag == 3)  	{
	  UNPROTECT(nprot);
	  error("Terminated: ill conditioned problem.\n");
	} else	{
/*                   ####   returning output   ####                           */
    nx = nmesh;

    PROTECT(yout = allocVector(REALSXP,(ncomp+1)*(nx))); nprot++;
	  for (j = 0; j < nx; j++)       REAL(yout)[j]    = xx[j];
    for (j = 0; j < ncomp*nx; j++) REAL(yout)[nx+j] =  u[j];
  }

  PROTECT(ISTATE = allocVector(INTSXP, 13));             nprot++;
  for (j = 0; j < 13; j++)
    INTEGER(ISTATE)[j] = 0;
  INTEGER(ISTATE)[0] = iflag;
  for (j = 0; j < 7; j++)
    INTEGER(ISTATE)[1+j] = icount[j];
  INTEGER(ISTATE)[9] = nmax;
  INTEGER(ISTATE)[10] = nmesh;
  INTEGER(ISTATE)[11] = lwrkfl;
  INTEGER(ISTATE)[12] = lwrkin;

  setAttrib(yout, install("istate"), ISTATE);

  PROTECT(EPSS = allocVector(REALSXP, 2));             nprot++;
  REAL(EPSS)[0] = epsini;
  REAL(EPSS)[1] = epsmin;

  setAttrib(yout, install("eps"), EPSS);

  PROTECT(RSTATE = allocVector(REALSXP, 5));          nprot++;
  REAL(RSTATE)[0] = ckappa1;
  REAL(RSTATE)[1] = gamma1;
  REAL(RSTATE)[2] = sigma;
  REAL(RSTATE)[3] = ckappa;
  REAL(RSTATE)[4] = ckappa2;
  setAttrib(yout, install("rstate"), RSTATE);

/*               ####   termination   ####                                    */
  UNPROTECT(nprot);
  return(yout);
}

#include <time.h>
#include <string.h>
#include "twpbvpc.h"

/* definition of the calls to the fortran functions - in file twpbvpc.f
*/

void F77_NAME(twpbvpc)(int*, int*, double *, double *,
         int *, double *, int *, int *, double *,
         int *, int *, int *, int *, int *,
         double *, int *, double *, int *,
         int *, double *, int*, int*, double *,
         void (*)(int *, double *, double *, double *, double *, int *), /* fsub(n,x,u,f,rp,ip)   */
		     void (*)(int *, double *, double *, double *, double *, int *), /* dfsub(n,x,u,df,rp,ip) */
			   void (*)(int *, int *, double *, double *, double *, int *),    /* gsub(i,n,u,g,rp,ip)   */
		     void (*)(int *, int *, double *, double *, double *, int *),    /* dgsub(i,n,u,dg,rp,ip) */
         double *, double *, double *, double *, double *, 
         double *, int *, int *, int *, int *, int *, int *, int *);

/* interface between fortran function calls and R functions */

static void C_bvpderiv_func (int *n,  double *x, double *y, double *ydot,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                             REAL(X)[0]   = *x;
  for (i = 0; i < *n ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang3(R_bvpderiv_func,X,Y)); incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));      incr_N_Protect();

  for (i = 0; i < *n ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
                                                my_unprotect(2);
}

/* wrapper above the derivate function that first estimates the
values of the forcing functions */

static void C_bvpderiv_func_forc (int *neq, double *x, double *y,
                         double *ydot, double *rpar, int *ipar)
{
  updatedeforc(x);
  derfun(neq, x, y, ydot, rpar, ipar);
}

/* interface between fortran call to jacobian and R function */
static void C_bvpjac_func (int *n,  double *x, double *y, double *pd,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                           REAL(X)[0]   = *x;
  for (i = 0; i < *n; i++) REAL(Y)[i]   = y[i];

  PROTECT(R_fcall = lang3(R_bvpjac_func,X,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));      incr_N_Protect();

  for (i = 0; i < *n * *n; i++)  pd[i] = REAL(ans)[i];
                                                my_unprotect(2);
}

/* interface between fortran call to boundary condition and corresponding R function */

static void C_bvpbound_func (int *ii, int *n, double *y, double *gout,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                             INTEGER(J)[0] = *ii;
  for (i = 0; i < *n ; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(R_bvpbound_func,J,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));        incr_N_Protect();
  /* only one element returned... */
  gout[0] = REAL(ans)[0];
                                                  my_unprotect(2);
}
/*interface between fortran call to jacobian of boundary and corresponding R function */

static void C_bvpjacbound_func (int *ii, int *n, double *y, double *dg,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, ans;
                           INTEGER(J)[0] = *ii;
  for (i = 0; i < *n; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(R_bvpjacbound_func,J,Y));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));          incr_N_Protect();

  for (i = 0; i < *n ; i++)  dg[i] = REAL(ans)[i];
                                                    my_unprotect(2);
}

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_bvptwp(SEXP Ncomp, SEXP Nlbc, SEXP Fixpnt, SEXP Aleft, SEXP Aright,
		SEXP Tol, SEXP Linear, SEXP Full, SEXP Givmesh, SEXP Givu, SEXP Nmesh,
		SEXP Nmax, SEXP Lwrkfl, SEXP Lwrkin, SEXP Xguess, SEXP Yguess,
    SEXP Rpar, SEXP Ipar, SEXP UseC, SEXP derivfunc, SEXP jacfunc, SEXP boundfunc,
    SEXP jacboundfunc, SEXP Initfunc, SEXP Parms, SEXP flist, SEXP rho)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, RSTATE;

  int  j, ii, ncomp, nlbc, nmax, lwrkfl, lwrkin, nx, *ipar, isForcing;
  double aleft, aright, *wrk, *tol, *fixpnt, *u, *xx, *rpar, *precis;
  double ckappa1, gamma1, sigma, ckappa, ckappa2; 
  int  liseries, *iseries, indnms, nxdim;
  int *ltol, *iwrk, ntol, iflag, nfixpnt, linear, givmesh, givu, nmesh, isDll;
  int full, useC;
  
  /* pointers to functions passed to FORTRAN */
  deriv_func_type    *deriv_func;
  jac_func_type      *jac_func;
  jacbound_func_type *jacbound_func;
  bound_func_type    *bound_func;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();
  useC   = INTEGER(UseC)[0];     /* conditioning or not */
  ncomp  = INTEGER(Ncomp)[0];    /* number of equations */
  nlbc   = INTEGER(Nlbc)[0];     /* number of left boundary conditions */
  nmax   = INTEGER(Nmax)[0];     /* max number of mesh points */
  lwrkfl = INTEGER(Lwrkfl)[0];   /* length of double workspace */
  lwrkin = INTEGER(Lwrkin)[0];   /* length of integer workspace */
  linear = INTEGER(Linear)[0];   /* true if linear problem */
  full = INTEGER(Full)[0];       /* true if full output */
  givu   = INTEGER(Givu)[0];     /* true if initial trial solution given */
  givmesh = INTEGER(Givmesh)[0]; /* true if initial mesh given */
  nmesh  = INTEGER(Nmesh)[0];    /* size of mesh */

  aleft  =REAL(Aleft)[0];
  aright =REAL(Aright)[0];

  /* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  /* copies of variables that will be changed in the FORTRAN subroutine */
  ntol = LENGTH(Tol);
  tol   =(double *) R_alloc(ntol, sizeof(double));
    for (j = 0; j < ntol;j++) tol[j] = REAL(Tol)[j];

  ltol   =(int *) R_alloc(ntol, sizeof(int));
    for (j = 0; j < ntol;j++) ltol[j] = j+1;

  nfixpnt =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(nfixpnt, sizeof(double));
   for (j = 0; j < nfixpnt;j++) fixpnt[j] = REAL(Fixpnt)[j];
  // check this - Francesca/Jeff:
  nxdim = nmax;
  xx   =(double *) R_alloc(nmax, sizeof(double));
   for (j = 0; j < nmesh; j++) xx[j] = REAL(Xguess)[j];
   for (j = nmesh; j < nmax; j++) xx[j] = 0;

  ii = nmax*ncomp;
  u   =(double *) R_alloc(ii, sizeof(double));
   for (j = 0; j < nmesh*ncomp; j++) u[j] = REAL(Yguess)[j];
   for (j = nmesh*ncomp; j < nmax*ncomp; j++) u[j] = 0;
      
  wrk = (double *) R_alloc(lwrkfl, sizeof(double));
     for (j = 0; j < lwrkfl; j++) wrk[j] = 0.;

  iwrk= (int *)    R_alloc(lwrkin, sizeof(int));
     for (j = 0; j < lwrkin; j++) iwrk[j] = 0;

  precis = (double *) R_alloc(3,sizeof(double));
  precis[0] = DBL_MIN;
  precis[1] = DBL_MAX;
  precis[2] = DBL_EPSILON/FLT_RADIX;
  
  ii = LENGTH(Ipar);
  ipar = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) ipar[j] = INTEGER(Ipar)[j];

  ii = LENGTH(Rpar);
  rpar = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rpar[j] = REAL(Rpar)[j];

  /* initialise global R-variables... */
  if (isDll == 0) {
    PROTECT(X  = NEW_NUMERIC(1));               incr_N_Protect();
    PROTECT(J = NEW_INTEGER(1));                incr_N_Protect();
    PROTECT(Y = allocVector(REALSXP,ncomp));    incr_N_Protect();
  }

  /* Initialization of Parameters and Forcings (DLL functions)  */
  isForcing = initForcings(flist);
  initParms(Initfunc, Parms);

  R_envir = rho;

  /* pointers to functions passed to FORTRAN */
  if (isDll) {   /* DLL addresses passed to fortran */
      deriv_func = (deriv_func_type *) R_ExternalPtrAddr(derivfunc);
      jac_func = (jac_func_type *) R_ExternalPtrAddr(jacfunc);
      bound_func = (bound_func_type *) R_ExternalPtrAddr(boundfunc);
      jacbound_func = (jacbound_func_type *) R_ExternalPtrAddr(jacboundfunc);

	  /* here overruling deriv_func if forcing */
      if (isForcing) {
        derfun = (deriv_func_type *) R_ExternalPtrAddr(derivfunc);
        deriv_func = (deriv_func_type *) C_bvpderiv_func_forc;
      }

  } else {      /* interface functions between fortran and R */
      deriv_func = (deriv_func_type *) C_bvpderiv_func;
      R_bvpderiv_func = derivfunc;

      jac_func = C_bvpjac_func;
      R_bvpjac_func = jacfunc;

      bound_func = C_bvpbound_func;
      R_bvpbound_func = boundfunc;

      jacbound_func = C_bvpjacbound_func;
      R_bvpjacbound_func = jacboundfunc;
    }

/* Call the fortran function twpbvpc
// CHECK liseries with jeff/francesca!

*/
  liseries = nmax;
  iseries = (int *)    R_alloc(liseries, sizeof(int));
    
	F77_CALL(twpbvpc) (&ncomp, &nlbc, &aleft, &aright, &nfixpnt, fixpnt,
        &ntol, ltol, tol, &linear, &givmesh, &givu, &nmesh, &nxdim, xx,
        &ncomp, u, &nmax, &lwrkfl, wrk, &lwrkin, iwrk, precis,
        deriv_func, jac_func, bound_func, jacbound_func, 
        &ckappa1, &gamma1, &sigma, &ckappa, &ckappa2, 
        rpar, ipar,
        &iflag, &liseries, iseries, &indnms, &full, &useC);

/*  iflag - The Mode Of Return From twpbvp  */
	if (iflag == 4)      {
	   unprotect_all();
     error("One of the input parameters is invalid.\n");
  }  else if (iflag == 1) 	{
	  unprotect_all();
	  error("The Expected No. Of mesh points Exceeds Storage Specifications.\n");
	}     else if (iflag == 2) 	{
	  unprotect_all();
	  error("The Expected No. Of meshes Exceeds Storage Specifications. Increase liseries\n");
	}     else if (iflag == 3)  	{
	  unprotect_all();
	  error("Terminated: ill conditioned problem.\n");
	}   else	{
  /*                   ####   returning output   ####                           */    
    nx = nmesh;

    PROTECT(yout = allocVector(REALSXP,(ncomp+1)*(nx)));incr_N_Protect();
	  for (j = 0; j < nx; j++)       REAL(yout)[j]    = xx[j];
    for (j = 0; j < ncomp*nx; j++) REAL(yout)[nx+j] =  u[j];
  }
 
  PROTECT(ISTATE = allocVector(INTSXP, 5));incr_N_Protect();
  INTEGER(ISTATE)[0] = iflag;
  INTEGER(ISTATE)[1] = nmax;
  INTEGER(ISTATE)[2] = nmesh;
  INTEGER(ISTATE)[3] = lwrkfl;
  INTEGER(ISTATE)[4] = lwrkin;
  setAttrib(yout, install("istate"), ISTATE);

  PROTECT(RSTATE = allocVector(REALSXP, 5));incr_N_Protect();
  REAL(RSTATE)[0] = ckappa1;
  REAL(RSTATE)[1] = gamma1; 
  REAL(RSTATE)[2] = sigma; 
  REAL(RSTATE)[3] = ckappa;
  REAL(RSTATE)[4] = ckappa2;
  setAttrib(yout, install("rstate"), RSTATE);
  
/*               ####   termination   ####                            */
  unprotect_all();
  return(yout);
}

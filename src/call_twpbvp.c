#include <time.h>
#include <string.h>
#include "bvpSolve.h"
#include "externalptr.h"

/* definition of the calls to the fortran functions - in files twpbvpc.f,
twpbvplc.f, twpbvpa.f */

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
       double *, int *, int *, int *, int *, int *, int *, int *,
       int *, double *, int *, double *, int*);

void F77_NAME(twpbvplc)(int*, int*, double *, double *,
       int *, double *, int *, int *, double *,
       int *, int *, int *, int *, int *,
       double *, int *, double *, int *,
       int *, double *, int*, int*, double *,
       void (*)(int *, double *, double *, double *, double *, int *), /* fsub(n,x,u,f,rp,ip)   */
		   void (*)(int *, double *, double *, double *, double *, int *), /* dfsub(n,x,u,df,rp,ip) */
		   void (*)(int *, int *, double *, double *, double *, int *),    /* gsub(i,n,u,g,rp,ip)   */
		   void (*)(int *, int *, double *, double *, double *, int *),    /* dgsub(i,n,u,dg,rp,ip) */
       double *, double *, double *, double *, double *,
       double *, int *, int *, int *, int *, int *, int *, int *,
       int *, double *, int *, double *, int*);

/* wrapper above the derivate function that first estimates the
values of the forcing functions - when model in compiled code */

static void dll_bvp_deriv_func_forc (int *neq, double *x, double *y,
                         double *ydot, double *rpar, int *ipar)
{
  updatedeforc(x);
  derfun(neq, x, y, ydot, rpar, ipar);
}

/* KARLINE -> FRANCESCA: NUMERICAL FUNCTIONS, IF NOT GIVEN */

static void C_num_jac_func (int *n,  double *x, double *y, double *pd,
                            double * rpar, int * ipar)
{
  int i, j;
  double perturb;

  for (i = 0; i < *n; i++) ycopy[i]   = y[i];

  jderfun(n, x, y, dy, rpar, ipar);
  for (i = 0; i < *n; i++) dycopy[i]   = dy[i];
  for (i = 0; i < *n * *n; i++) pd[i] = 0.;
  
  for (j = 0; j < *n; j++) {
     if (y[j] > 1.)
       perturb = y[j]*1e-8;
     else
       perturb = 1e-8 ;  
     ycopy[j] = y[j] + perturb;
     
     jderfun(n, x, ycopy, dycopy, rpar, ipar);
     
     ycopy[j] = y[j];
     
     for (i = 0; i < *n; i++) 
       pd[j* *n + i] = (dycopy[i] - dy[i])/perturb;

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
  double g, gcopy;
//  warning("entering numerical BOUNDARY jacobian %i %i %g %g\n", *ii, *n, y[0], dg[0]);

  for (i = 0; i < *n; i++) ycopy[i]  = y[i];
  for (i = 0; i < *n; i++) dg[i] = 0.;
  for (i = 0; i < *n; i++) {
    jbndfun(ii, n, y, &g, rpar, ipar);

    if (y[i] > 1.)
       perturb = y[i]*1e-8;
    else
       perturb = 1e-8;  

    ycopy[i] = y[i] + perturb;
    jbndfun(ii, n, ycopy, &gcopy, rpar, ipar);
    ycopy[i] = y[i];
    dg[i] = (gcopy - g)/perturb;
  }
}


/* interface between fortran function calls and R functions */

static void C_bvp_deriv_func (int *n,  double *x, double *y, double *ydot,
  double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, X, ans;

  for (i = 0; i < *n ; i++)  REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                    
  PROTECT(R_fcall = lang3(R_bvp_deriv_func,X,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));          

  for (i = 0; i < *n ; i++) ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
  
  UNPROTECT(3);
}

/* interface between fortran call to jacobian and R function */
static void C_bvp_jac_func (int *n,  double *x, double *y, double *pd,
                            double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, X, ans;

  for (i = 0; i < *n; i++) REAL(Y)[i]   = y[i];

  PROTECT(X = ScalarReal(*x));                 
  PROTECT(R_fcall = lang3(R_bvp_jac_func,X,Y));
  PROTECT(ans = eval(R_fcall, R_envir));       

  for (i = 0; i < *n * *n; i++)  pd[i] = REAL(ans)[i];
  
  UNPROTECT(3);
}


/* interface between fortran call to boundary condition and corresponding R function */

static void C_bvp_bound_func (int *ii, int *n, double *y, double *gout,
                              double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, J, ans;

  for (i = 0; i < *n ; i++)  REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                
  PROTECT(R_fcall = lang3(R_bvp_bound_func,J,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));          
  
  gout[0] = REAL(ans)[0];       /* only one element returned... */
  
  UNPROTECT(3);
}

/*interface between fortran call to jacobian of boundary and corresponding R function */

static void C_bvp_jacbound_func (int *ii, int *n, double *y, double *dg,
                                 double * rpar, int * ipar)
{
  int i;
  SEXP R_fcall, J, ans;

  for (i = 0; i < *n; i++) REAL(Y)[i] = y[i];

  PROTECT(J = ScalarInteger(*ii));                   
  PROTECT(R_fcall = lang3(R_bvp_jacbound_func,J,Y)); 
  PROTECT(ans = eval(R_fcall, R_envir));             

  for (i = 0; i < *n ; i++)  dg[i] = REAL(ans)[i];
  UNPROTECT(3);
}


/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_bvptwp(SEXP Ncomp, SEXP Fixpnt, SEXP Aleft, SEXP Aright,
		SEXP Nlbc, SEXP Tol, SEXP Linear, SEXP Full, 
    SEXP Givmesh, SEXP Givu, SEXP Nmesh,
		SEXP Nmax, SEXP Lwrkfl, SEXP Lwrkin, SEXP Xguess, SEXP Yguess,
    SEXP Rpar, SEXP Ipar, SEXP UseC, 
    SEXP derivfunc, SEXP jacfunc, SEXP boundfunc,
    SEXP jacboundfunc, SEXP Initfunc, SEXP Parms, SEXP flist, 
    SEXP Lobatto, SEXP Absent, SEXP Rwork, SEXP rho)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout=NULL, ISTATE, RSTATE;

  int  j, ii, ncomp, nlbc, nmax, lwrkfl, lwrkin, nx, *ipar, isForcing;
  double *wrk, *tol, *fixpnt, *u, *xx, *rpar, *precis, *xguess, *yguess;
  double aleft, aright, ckappa1, gamma1, sigma, ckappa, ckappa2; 
  int  liseries, *iseries, indnms, nxdim; // type;
  int *ltol, *iwrk, *iset, ntol, iflag, nfixpnt, linear, givmesh, 
      givu, nmesh, isDll, nugdim, nmshguess, lobatto;
  int *absent;
  double *rwork;
  int full, useC;
  
  /* pointers to functions passed to FORTRAN */
  C_deriv_func_type    *deriv_func;
  C_jac_func_type      *jac_func = NULL;
  C_bound_func_type    *bound_func = NULL;
  C_jacbound_func_type *jacbound_func = NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  int nprot = 0;

  aleft  =REAL(Aleft)[0];
  aright =REAL(Aright)[0];

  ncomp   = INTEGER(Ncomp)[0];    /* number of equations */
  lobatto = INTEGER(Lobatto)[0];  /* 0 = bvptwp , 1 = bvptwpl*/

  ii = LENGTH(Absent);
  absent = (int *) R_alloc(ii, sizeof(int));
     for (j=0; j<ii; j++) absent[j] = INTEGER(Absent)[j];

  ii = LENGTH(Rwork);
  rwork = (double *) R_alloc(ii, sizeof(double));
     for (j=0; j<ii; j++) rwork[j] = REAL(Rwork)[j];

  nlbc   = INTEGER(Nlbc)[0];     /* number of left boundary conditions */
  nmax   = INTEGER(Nmax)[0];     /* max number of mesh points */
  lwrkfl = INTEGER(Lwrkfl)[0];   /* length of double workspace */
  lwrkin = INTEGER(Lwrkin)[0];   /* length of integer workspace */
  linear = INTEGER(Linear)[0];   /* true if linear problem */
  full   = INTEGER(Full)[0];     /* true if full output */
  givu   = INTEGER(Givu)[0];     /* true if initial trial solution given */
  givmesh = INTEGER(Givmesh)[0]; /* true if initial mesh given */
  nmesh  = INTEGER(Nmesh)[0];    /* size of mesh */
  useC   = INTEGER(UseC)[0];     /* conditioning or not */

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
  for (j = 0; j < ntol; j++) ltol[j] = j+1;

  nfixpnt =  LENGTH(Fixpnt);
  fixpnt   =(double *) R_alloc(nfixpnt, sizeof(double));
  for (j = 0; j < nfixpnt;j++) fixpnt[j] = REAL(Fixpnt)[j];

  // check this - Francesca/Jeff:
  nxdim = nmax;
  xx   =(double *) R_alloc(nmax, sizeof(double));
  if (givmesh) {
  for (j = 0; j < nmesh; j++) xx[j] = REAL(Xguess)[j];
  } else {
  for (j = 0; j < nmesh; j++) xx[j] = 0;
  }
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
   if (givu) {
   for (j = 0; j < nmesh*ncomp; j++) u[j] = REAL(Yguess)[j];
   } else { 
   for (j = 0; j < nmesh*ncomp; j++) u[j] = 0;
   }
   for (j = nmesh*ncomp; j < nmax*ncomp; j++) u[j] = 0;
      
  wrk = (double *) R_alloc(lwrkfl, sizeof(double));
     for (j = 0; j < lwrkfl; j++) wrk[j] = 0.;

  iwrk = (int *)   R_alloc(lwrkin, sizeof(int));
     for (j = 0; j < lwrkin; j++) iwrk[j] = 0;

  iset = (int*)    R_alloc(6, sizeof(int));
     for (j = 0; j < 6; j++) iset[j] = 0;
  
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
    PROTECT(Y = allocVector(REALSXP,ncomp));    nprot++;
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

  /* pointers to functions passed to FORTRAN */
  if (isDll) {   /* DLL addresses passed to fortran */
      deriv_func    = (C_deriv_func_type *)    R_ExternalPtrAddrFn_(derivfunc);
      
      if (absent[0] == 0)  
        jac_func      = (C_jac_func_type *)   R_ExternalPtrAddrFn_(jacfunc);
      
      if (absent[1] == 0)
        bound_func    = (C_bound_func_type *)  R_ExternalPtrAddrFn_(boundfunc);
      
      if (absent[2] == 0)   /* not given*/
        jacbound_func = (C_jacbound_func_type *) R_ExternalPtrAddrFn_(jacboundfunc);

	  /* here overruling deriv_func if forcing */
      if (isForcing) {
        derfun =     (C_deriv_func_type *) R_ExternalPtrAddrFn_(derivfunc);
        deriv_func = (C_deriv_func_type *) dll_bvp_deriv_func_forc;
      }
      
  } else {      /* interface functions between fortran and R */
      deriv_func = C_bvp_deriv_func;
      R_bvp_deriv_func = derivfunc;

      if (absent[0] == 0) {
        jac_func = C_bvp_jac_func;
        R_bvp_jac_func = jacfunc;
      } 

      if (absent[1] == 0) {
        bound_func = C_bvp_bound_func;
        R_bvp_bound_func = boundfunc;
      }
      
      if (absent[2] == 0) {
        jacbound_func = C_bvp_jacbound_func;
        R_bvp_jacbound_func = jacboundfunc;
      }
    }
/* if numerical approximates should be used */    
    if (absent[0] == 1) {
        dy     = (double *) R_alloc(ncomp, sizeof(double));
        dycopy = (double *) R_alloc(ncomp, sizeof(double));
        ycopy  = (double *) R_alloc(ncomp, sizeof(double)); 
        jac_func = (C_jac_func_type *)      C_num_jac_func;
        jderfun  = deriv_func;
    }
    if (absent[1] == 1) {
      bound_func = (C_bound_func_type *) C_num_bound_func;
      iibb = (int *) R_alloc(ncomp, sizeof(int));
      for (j = 0; j < ncomp; j++)
        iibb[j] = absent[3 + j];
      bb = (double *) R_alloc(ncomp, sizeof(double));
      for (j = 0; j < ncomp; j++)
        bb[j] = rwork[j];
    }
    if (absent[2] == 1) {
        jacbound_func = (C_jacbound_func_type *) C_num_jacbound_func;
        jbndfun = bound_func;
        if (absent[0] != 1)
          ycopy  = (double *) R_alloc(ncomp, sizeof(double)); 
    }
    
    
/* Call the fortran function twpbvpc
// CHECK liseries with jeff/francesca!
// francesca liseries is the maximum possible different mesh, 
// it is not related to nmax, usually 500 is ok.

*/
  liseries = 500;
  iseries = (int *)    R_alloc(liseries, sizeof(int));
  nugdim = ncomp;
  if (givu){
  nmshguess = nmesh;
  } else
  { 
  nmshguess = 1;
  }

  if (lobatto == 1)
	F77_CALL(twpbvplc) (&ncomp, &nlbc, &aleft, &aright, &nfixpnt, fixpnt,
        &ntol, ltol, tol, &linear, &givmesh, &givu, &nmesh, &nxdim, xx,
        &ncomp, u, &nmax, &lwrkfl, wrk, &lwrkin, iwrk, precis,
        deriv_func, jac_func, bound_func, jacbound_func,
        &ckappa1, &gamma1, &sigma, &ckappa, &ckappa2,
        rpar, ipar, &iflag, &liseries, iseries, &indnms,
        &full, &useC, &nmshguess, xguess, &nugdim,
        yguess, iset);
  else
	F77_CALL(twpbvpc) (&ncomp, &nlbc, &aleft, &aright, &nfixpnt, fixpnt,
        &ntol, ltol, tol, &linear, &givmesh, &givu, &nmesh, &nxdim, xx,
        &ncomp, u, &nmax, &lwrkfl, wrk, &lwrkin, iwrk, precis,
        deriv_func, jac_func, bound_func, jacbound_func, 
        &ckappa1, &gamma1, &sigma, &ckappa, &ckappa2, 
        rpar, ipar, &iflag, &liseries, iseries, &indnms, 
        &full, &useC, &nmshguess, xguess, &nugdim,
        yguess, iset);

/*  iflag - The Mode Of Return From twpbvp  */
	if (iflag == 4)      {
	  UNPROTECT(nprot);
	  error("One of the input parameters is invalid.\n");
  }  else if (iflag == 1) 	{
    UNPROTECT(nprot);
    error("The Expected No. Of mesh points Exceeds Storage Specifications.\n");
	}     else if (iflag == 2) 	{
	  UNPROTECT(nprot);
	  error("The Expected No. Of meshes Exceeds Storage Specifications. Increase liseries\n");
	}     else if (iflag == 3)  	{
	  UNPROTECT(nprot);
	  error("Terminated: ill conditioned problem.\n");
	}   else	{
  /*                   ####   returning output   ####                           */    
    nx = nmesh;

    PROTECT(yout = allocVector(REALSXP,(ncomp+1)*(nx)));   nprot++;
	  for (j = 0; j < nx; j++)       REAL(yout)[j]    = xx[j];
    for (j = 0; j < ncomp*nx; j++) REAL(yout)[nx+j] =  u[j];
  }
 
  PROTECT(ISTATE = allocVector(INTSXP, 13));   nprot++;
  INTEGER(ISTATE)[0] = iflag;
  for (j = 0; j < 6; j++)
    INTEGER(ISTATE)[1+j] = iset[j];
  INTEGER(ISTATE)[7]=0;
  INTEGER(ISTATE)[8]=0;  
  INTEGER(ISTATE)[9] = nmax;
  INTEGER(ISTATE)[10] = nmesh;
  INTEGER(ISTATE)[11] = lwrkfl;
  INTEGER(ISTATE)[12] = lwrkin;
    
  setAttrib(yout, install("istate"), ISTATE);

  PROTECT(RSTATE = allocVector(REALSXP, 5));   nprot++;
  REAL(RSTATE)[0] = ckappa1;
  REAL(RSTATE)[1] = gamma1; 
  REAL(RSTATE)[2] = sigma; 
  REAL(RSTATE)[3] = ckappa;
  REAL(RSTATE)[4] = ckappa2;
  setAttrib(yout, install("rstate"), RSTATE);
  
/*               ####   termination   ####                            */
  UNPROTECT(nprot);
  return(yout);
}


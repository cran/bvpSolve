/* Define some global variables and functions that operate on some of them */
/* Patterned on code odesolve_utils.c from package odesolve */
#include <R.h>
#include <Rdefines.h>
#include "bvpSolve.h"
#include "externalptr.h"

/* some functions for keeping track of how many SEXPs 
 * 	are PROTECTed, and UNPROTECTing them in the case of a fortran stop.
 */
 
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}
/* */
int Leftbc;
double A;
double B;
double yguess;

/* Globals : the R-functions, R-environment, R-parameters */

SEXP R_bvp_deriv_func;
SEXP R_bvp_jac_func;
SEXP R_bvp_bound_func;
SEXP R_bvp_jacbound_func;
SEXP R_bvp_guess_func;
SEXP R_cont_deriv_func;
SEXP R_cont_jac_func;
SEXP R_cont_bound_func;
SEXP R_cont_jacbound_func;
SEXP R_cont_guess_func;

SEXP R_envir;
SEXP bvp_gparms;

/*==================================================
Parameter initialisation functions
note: forcing initialisation ftion is in forcings.c
===================================================*/

void initParms(SEXP Initfunc, SEXP Parms) {
  // ks: added this to prevent entering this if initfunc does not exist
  if (Initfunc == NA_STRING) return;
  if (inherits(Initfunc, "NativeSymbol"))  {
    init_func *initializer;

    PROTECT(bvp_gparms = Parms);     incr_N_Protect();
    initializer = (init_func *) R_ExternalPtrAddrFn_(Initfunc);
    initializer(Initbvpparms);
  }
}

void Initbvpparms(int *N, double *parms) {
  int i, Nparms;

  Nparms = LENGTH(bvp_gparms);
  if ((*N) != Nparms)  {
      warning("Number of parameters passed to solver, %i; number in DLL, %i\n",Nparms, *N);
      PROBLEM "Confusion over the length of parms"
      ERROR;
  }  else  {
      for (i = 0; i < *N; i++) parms[i] = REAL(bvp_gparms)[i];
      epsval = parms;     /* set pointer to c globals or fortran common block */
  }
}
  
SEXP get_bvpSolve_gparms(void) {
  return bvp_gparms;
}

/*==================================================
 extracting elements from a list
===================================================*/

SEXP getListElement(SEXP list, const char *str)  {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
	 if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	 }
   return elmt;
}

/*==================================================
 Slightly different printing, based on R-core dblep0 and intp0
===================================================*/


void F77_NAME(dblep0k) (const char *label, int *nchar, double *data, int *ndata)
{
    int k, nc = *nchar;

    if(nc < 0) nc = strlen(label);
    if(nc > 255) {
	warning("invalid character length in dblepr");
	nc = 0;
    } else if(nc > 0) {
	for (k = 0; k < nc; k++)
	    Rprintf("%c", label[k]);
    }
/*    if(*ndata > 0) printRealVector(data, *ndata, 1); */
    
    if(*ndata > 0) 
      for (k = 0; k < *ndata; k++)  
        Rprintf("%g",data[k]);
	Rprintf("\n");/* put that at end, so that vector is on same line as string*/
/*    return(0);  and removed this one */
}

void F77_NAME(intp0k) (const char *label, int *nchar, int *data, int *ndata)
{
    int k, nc = *nchar;

    if(nc < 0) nc = strlen(label);
    if(nc > 255) {
	warning("invalid character length in intpr");
	nc = 0;
    } else if(nc > 0) {
	for (k = 0; k < nc; k++)
	    Rprintf("%c", label[k]);
    }
/*    if(*ndata > 0) printIntegerVector(data, *ndata, 1); */
      for (k = 0; k < *ndata; k++)  
        Rprintf("%i",data[k]);
    	Rprintf("\n");
/*    return(0);  and removed this one */
}

/* Define some global variables and functions that operate on some of them */
/* Patterned on code desolve_utils.c from package desolve */
#include <R.h>
#include <Rdefines.h>
#ifndef  R_INTERNALS_H_
#include <Rinternals.h>
#endif

#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>

#include "bvpSolve.h"
#include "externalptr.h"

/* functions for keeping track of how many SEXPs not used anymore.
long int N_Protected;
void init_N_Protect(void) { N_Protected = 0; }
void incr_N_Protect(void) { N_Protected++; }
void unprotect_all(void) { UNPROTECT((int) N_Protected); }
void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}  */
/* */
int Leftbc;
double A;
double B;
double yguess;

/* Globals : the R-functions, R-environment, R-parameters */

//SEXP R_bvp_deriv_func;
//SEXP R_bvp_jac_func;
//SEXP R_bvp_bound_func;
//SEXP R_bvp_jacbound_func;
//SEXP R_bvp_guess_func;
//SEXP R_cont_deriv_func;
//SEXP R_cont_jac_func;
//SEXP R_cont_bound_func;
//SEXP R_cont_jacbound_func;
//SEXP R_cont_guess_func;

//SEXP R_envir;
//SEXP bvp_gparms;

/*==================================================
Parameter initialisation functions - removed
note: forcing initialisation ftion is in forcings.c
===================================================*/


void Initbvpparms(int *N, double *parms) {
  int i, Nparms;

  Nparms = LENGTH(bvp_gparms);
  if ((*N) != Nparms)  {
      warning("Number of parameters passed to solver, %i; number in DLL, %i\n",Nparms, *N);
      Rf_error("Confusion over the length of parms.");
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


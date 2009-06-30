/* Define some global variables and functions that operate on some of them */
/* Patterned on code odesolve_utils.c from package odesolve */
#include <R.h>
#include <Rdefines.h>

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

SEXP colsys_deriv_func;
SEXP colsys_jac_func;
SEXP colsys_bound_func;
SEXP colsys_jacbound_func;

SEXP col_deriv_func;
SEXP col_jac_func;
SEXP col_bound_func;
SEXP col_jacbound_func;

SEXP bvp_envir;
SEXP bvpcolmod_envir;
SEXP bvp_gparms;

void Initbvpparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(bvp_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(bvp_gparms)[i];
    }
}
  
SEXP get_bvpSolve_gparms(void)
{
  return bvp_gparms;
}

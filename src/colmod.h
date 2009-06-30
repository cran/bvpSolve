#include <R.h>
#include <Rdefines.h>
/* global variables */
SEXP X, Y, J, EPS;

/* bvp globals */
extern SEXP colsys_deriv_func;
extern SEXP colsys_jac_func;
extern SEXP colsys_bound_func;
extern SEXP colsys_jacbound_func;
extern SEXP bvpcolmod_envir;

/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

/* declarations for initibvpparms;*/
void Initbvpparms(int *, double *);

/* use in colmod */
long int n_eq;   /* number of equations */
long int ng;     /* number of boundary conditions */
long int mstar;  /* number of derivates (including higher order ones) */
long int ml;
long int nrowpd;

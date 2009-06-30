#include <R.h>
#include <Rdefines.h>
/* global variables */
SEXP X, Y, J;
extern SEXP bvp_gparms;

/* bvp globals */
extern SEXP col_deriv_func;
extern SEXP col_jac_func;
extern SEXP col_bound_func;
extern SEXP col_jacbound_func;
extern SEXP bvp_envir;

extern SEXP bvp_envir;
/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

/* declarations for initibvpparms;*/
void Initbvpparms(int *, double *);


#include <R.h>
#include <Rdefines.h>


/*============================================================================
  global R variables 
============================================================================*/
SEXP X, Y, J, EPS;

/*============================================================================
  global C variables 
============================================================================*/
long int nforc;  /* the number of forcings */
/* use in colnew */
int n_eq;   /* number of equations */
int ng;     /* number of boundary conditions */
int mstar;  /* number of derivates (including higher order ones) */
int ml;
int nrowpd;

/* Input data. three vectors:
  tmat, fmat: time, forcing function data value
  imat: index to start of each forcing function in tmat, fmat*/
double * tvec;
double * fvec;
int    * ivec;
int    fmethod;

/* for each forcing function: index to current position in tmat, fmat,
 current value, interpolation factor, current forcing time, next forcing time,
 max time (to be removed).....
*/
int    * findex;
double * intpol;
int    * maxindex;

double * forcings;
double * epsval;    /* when eps and model in compiled code */

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type    (int *, double *,double *, double *, double *, int *);
typedef void C_bound_func_type    (int *, int *, double *,double *, double *, int *);
typedef void C_jac_func_type      (int *, double *, double *,double *, double *, int *);
typedef void C_jacbound_func_type (int *, int *, double *, double *, double *, int *);
typedef void C_guess_func_type    (double *, double *, double *);

C_deriv_func_type    *derfun;    /* if DLL */
C_bound_func_type    *boundfun;
C_jac_func_type      *jacfun;
C_jacbound_func_type *jacboundfun;

/*============================================================================
  solver R- global functions 
============================================================================*/

extern SEXP bvp_gparms;

/* bvp globals */
extern SEXP R_bvp_deriv_func;
extern SEXP R_bvp_jac_func;
extern SEXP R_bvp_bound_func;
extern SEXP R_bvp_jacbound_func;
extern SEXP R_bvp_guess_func;
extern SEXP R_cont_deriv_func;
extern SEXP R_cont_jac_func;
extern SEXP R_cont_bound_func;
extern SEXP R_cont_jacbound_func;
extern SEXP R_cont_guess_func;

extern SEXP R_envir;


/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

/* declarations for initibvpparms;*/
void Initbvpparms(int *, double *);
void initParms(SEXP Initfunc, SEXP Parms);
typedef void init_func (void (*)(int *, double *));

/* forcings */
void updatedeforc(double *);
void Initdeforc(int *, double *);

int initForcings(SEXP list);

SEXP getListElement(SEXP list, const char *str);


#ifndef R_R_H
#  include <R.h>
#endif

#ifndef R_DEFINES_H
#  include <Rdefines.h>
#endif

#ifndef R_INTERNALS_H_
#  include <Rinternals.h>
#endif

#ifndef EXTERN
# define EXTERN extern
#endif 

/*============================================================================
  global R variables 
============================================================================*/
EXTERN  SEXP Y, EPS;
EXTERN  SEXP yout, ISTATE, RSTATE;

/*============================================================================
  global C variables 
============================================================================*/
EXTERN  long int nforc;  /* the number of forcings */
/* use in colnew */
EXTERN  int n_eq;   /* number of equations */
EXTERN  int ng;     /* number of boundary conditions */
EXTERN  int mstar;  /* number of derivates (including higher order ones) */
EXTERN  int ml;
EXTERN  int nrowpd;

/* Input data. three vectors:
  tmat, fmat: time, forcing function data value
  imat: index to start of each forcing function in tmat, fmat*/
EXTERN  double * tvec;
EXTERN  double * fvec;
EXTERN  int    * ivec;
EXTERN  int    fmethod;

/* for each forcing function: index to current position in tmat, fmat,
 current value, interpolation factor, current forcing time, next forcing time,
 max time (to be removed).....
*/
EXTERN  int    * findex;
EXTERN  double * intpol;
EXTERN  int    * maxindex;

EXTERN  double * forcings;
EXTERN  double * epsval;    /* when eps and model in compiled code */

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type    (int *, double *,double *, double *, double *, int *);
typedef void C_bound_func_type    (int *, int *, double *,double *, double *, int *);
typedef void C_jac_func_type      (int *, double *, double *,double *, double *, int *);
typedef void C_jacbound_func_type (int *, int *, double *, double *, double *, int *);
typedef void C_guess_func_type    (double *, double *, double *);

EXTERN  C_deriv_func_type    *jderfun;    /* if DLL */
EXTERN  C_deriv_func_type    *derfun;    /* if DLL */
EXTERN  C_bound_func_type    *boundfun;
EXTERN  C_bound_func_type    *jbndfun;
EXTERN  C_jac_func_type      *jacfun;
EXTERN  C_jacbound_func_type *jacboundfun;

typedef void C_acdc_deriv_func_type(int *, double *, double *,double *,
                                    double *, double *, int *);
typedef void C_acdc_bound_func_type(int *, int *, double *, double *,
                                    double *, double *, int *);
typedef void C_acdc_jac_func_type  (int *,  double *, double *, double *,
                                    double *, double *, int *);
typedef void C_acdc_jacbound_func_type(int *, int *, double *, double *,
                                    double *, double *, int *);
EXTERN  C_acdc_deriv_func_type    *jaderfun;    /* if DLL */
EXTERN  C_acdc_bound_func_type    *jabndfun;

typedef void C_deriv_func2_type(double *, double *,double *, double *,
                                double *, int *);
typedef void C_bound_func2_type(int *, double *, double *,double *,
                                double *, int *);
typedef void C_jac_func2_type  (double *, double *, double *, int *, double *,
                                double *, int *);
typedef void C_jacbound_func2_type(int *, double *, double *, double *,
                                double *, int *);
typedef void C_guess_func2_type(double *, double *, double *, double *, 
                                double *, int *);
EXTERN  C_deriv_func2_type     *jepsderfun;
EXTERN  C_bound_func2_type     *jepsbndfun;
EXTERN  double *dy, *dycopy, *ycopy, *ycopy2, *bb, *g, *gcopy;
EXTERN  int * iibb;
EXTERN  int  isDll;


/*============================================================================
  solver R- global functions 
============================================================================*/

EXTERN  SEXP bvp_gparms;
EXTERN  SEXP getListElement(SEXP list, const char *str);

/* bvp globals */
EXTERN  SEXP R_bvp_deriv_func;
EXTERN  SEXP R_bvp_jac_func;
EXTERN  SEXP R_bvp_bound_func;
EXTERN  SEXP R_bvp_jacbound_func;
EXTERN  SEXP R_bvp_guess_func;
EXTERN  SEXP R_cont_deriv_func;
EXTERN  SEXP R_cont_jac_func;
EXTERN  SEXP R_cont_bound_func;
EXTERN  SEXP R_cont_jacbound_func;
EXTERN  SEXP R_cont_guess_func;

EXTERN  SEXP R_envir;


/* utilities -not used anymore
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);
void initParms(SEXP Initfunc, SEXP Parms);
*/

/* declarations for initibvpparms;*/
void Initbvpparms(int *, double *);
typedef void init_func (void (*)(int *, double *));

/* forcings */
void updatedeforc(double *);
void Initdeforc(int *, double *);

int initForcings(SEXP list);


#undef EXTERN

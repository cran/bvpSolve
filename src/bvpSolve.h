#ifndef R_R_H
#  include <R.h>
#endif

#ifndef R_DEFINES_H
#  include <Rdefines.h>
#endif

#ifndef R_INTERNALS_H_
#  include <Rinternals.h>
#endif


/*============================================================================
  global R variables 
============================================================================*/
extern SEXP Y, EPS;

/*============================================================================
  global C variables 
============================================================================*/
extern long int nforc;  /* the number of forcings */
/* use in colnew */
extern int n_eq;   /* number of equations */
extern int ng;     /* number of boundary conditions */
extern int mstar;  /* number of derivates (including higher order ones) */
extern int ml;
extern int nrowpd;

/* Input data. three vectors:
  tmat, fmat: time, forcing function data value
  imat: index to start of each forcing function in tmat, fmat*/
extern double * tvec;
extern double * fvec;
extern int    * ivec;
extern int    fmethod;

/* for each forcing function: index to current position in tmat, fmat,
 current value, interpolation factor, current forcing time, next forcing time,
 max time (to be removed).....
*/
extern int    * findex;
extern double * intpol;
extern int    * maxindex;

extern double * forcings;
extern double * epsval;    /* when eps and model in compiled code */

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type    (int *, double *,double *, double *, double *, int *);
typedef void C_bound_func_type    (int *, int *, double *,double *, double *, int *);
typedef void C_jac_func_type      (int *, double *, double *,double *, double *, int *);
typedef void C_jacbound_func_type (int *, int *, double *, double *, double *, int *);
typedef void C_guess_func_type    (double *, double *, double *);

extern C_deriv_func_type    *jderfun;    /* if DLL */
extern C_deriv_func_type    *derfun;    /* if DLL */
extern C_bound_func_type    *boundfun;
extern C_bound_func_type    *jbndfun;
extern C_jac_func_type      *jacfun;
extern C_jacbound_func_type *jacboundfun;

typedef void C_acdc_deriv_func_type(int *, double *, double *,double *,
                                    double *, double *, int *);
typedef void C_acdc_bound_func_type(int *, int *, double *, double *,
                                    double *, double *, int *);
typedef void C_acdc_jac_func_type  (int *,  double *, double *, double *,
                                    double *, double *, int *);
typedef void C_acdc_jacbound_func_type(int *, int *, double *, double *,
                                    double *, double *, int *);
extern C_acdc_deriv_func_type    *jaderfun;    /* if DLL */
extern C_acdc_bound_func_type    *jabndfun;

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
extern C_deriv_func2_type     *jepsderfun;
extern C_bound_func2_type     *jepsbndfun;
extern double *dy, *dycopy, *ycopy, *ycopy2, *bb, *g, *gcopy;
extern int * iibb;

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

SEXP getListElement(SEXP list, const char *str);


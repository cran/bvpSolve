\name{bvptwp}
\alias{bvptwp}
\title{
  Solves two-point boundary value problems of ordinary differential
  equations, using a mono-implicit Runge-Kutta formula
}

\description{
  Solves a boundary value problem of a system of ordinary differential
  equations. This is an implementation of the fortran
  code twpbvpc, based on mono-implicit Runge-Kutta formulae
  of orders 4, 6 and 8 in a deferred correction framework and 
  that uses conditioning in the mesh selection.

  written by J.R. Cash, F. Mazzia and M.H. Wright.

  Rather than MIRK, it is also possible to select a lobatto method. This
  then uses the code 'twpbvplc', written  by Cash and Mazzia.

  It is possible to solve stiff systems, by using an automatic continuation
  strategy. This then uses the code 'acdc'.
}

\usage{
bvptwp(yini = NULL, x, func, yend = NULL, parms = NULL, 
       order = NULL, ynames = NULL, xguess = NULL, yguess = NULL, 
       jacfunc = NULL, bound = NULL, jacbound = NULL, 
       leftbc = NULL, posbound = NULL, islin = FALSE, nmax = 1000, 
       ncomp = NULL, atol = 1e-8, cond = FALSE, lobatto = FALSE,
       allpoints = TRUE, dllname = NULL, initfunc = dllname,
       rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, 
       forcings = NULL, initforc = NULL, fcontrol = NULL, 
       verbose = FALSE, epsini = NULL, eps = epsini, ...)
}
\arguments{
  \item{yini }{either a vector with the initial (state) variable values for
    the ODE system, or \code{NULL}.

    If \code{yini} is a vector, use \code{NA} for an
    initial value which is not specified. 
    
    If \code{yini} has a
    \code{names} attribute, the names will be available within \code{func} and
    used to label the output matrix.
    
    If \code{yini = NULL}, then the boundary conditions must be specified
    via function \code{bound};
    if not \code{NULL} then \code{yend} should also be not \code{NULL}.
  }
  \item{x }{sequence of the independent variable for which output is wanted;
    the first value of \code{x} must be the initial value (at which
    \code{yini} is defined), the final value the end condition (at which
    \code{yend} is defined).
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the model definition) at point \code{x},
    or a character string giving the name of a compiled function in a
    dynamically loaded shared library.
    
    If \code{func} is an \R-function, it must be defined as:
    \code{func = function(x, y, parms,...)}. \code{x} is the current point of
    the independent variable in the integration, \code{y} is the current
    estimate of the (state) variables in the ODE system. If the initial
    values \code{yini} has a names attribute, the names will be available
    inside \code{func}.  \code{parms} is a vector or list of parameters;
    ... (optional) are any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{x}, and whose next elements are global values that are required at
    each point in \code{x}.
    
    If \code{func} is a string, then \code{dllname} must give the name of the
    shared library (without extension) which must be loaded before
    \code{bvptwp} is called. See package vignette "bvpSolve" for more details.
  }
  \item{yend }{either a vector with the final (state) variable values for the
    ODE system, or \code{NULL};

    if \code{yend} is a vector, use \code{NA}
    for a final value which is not specified.
    
    If \code{yend} has a  \code{names} attribute, and \code{yini} does not,
    the names will be available within the functions and 
    used to label the output matrix.

    If \code{yend = NULL}, then the boundary conditions must be specified
    via function \code{bound}; 
    if not \code{NULL} then \code{yini} should also be not \code{NULL}.
  }
  \item{parms }{vector or a list with parameters passed to \code{func},
    \code{jacfunc}, \code{bound} and \code{jacbound} (if present).

    If \code{eps} is given a value then it should be the **first** element in
    \code{parms}.
  }
  \item{epsini }{the initial value of the continuation parameter. If
    \code{NULL} and \code{eps} is given a value, then \code{epsini} takes
    the default starting value of 0.5. For many singular perturbation type
    problems, the choice of 0.1 < \code{eps} < 1 represents a (fairly) easy
    problem. The user should attempt to specify an initial
    problem that is not `too' challenging.
    \code{epsini} must be initialised strictly less than 1 and greater than 0.
  }
  \item{eps }{the desired value of precision for which the user would like
    to solve the problem. \code{eps} must be less than or equal to
    \code{epsini}. If this is given a value, it must be the first value in
    \code{parms}.
  }
  \item{ynames }{The names of the variables; used to label the output, and
    avaliable within the functions.
    
    If \code{ynames} is \code{NULL}, names can also be passed via \code{yini},  
    \code{yend} or \code{yguess}.
      }
  \item{xguess }{Initial grid \code{x}, a vector. \code{bvptwp} requires the length
    of \code{xguess} to be at least equal to the length of \code{x}. If this
    is not the case, then \code{xguess} and \code{yguess} will be interpolated
    to \code{x} and a warning printed. If \code{xguess} is given, so should
    \code{yguess} be.
    
    Supplying \code{xguess} and \code{yguess}, based on results from a previous
    (simpler) BVP-ODE can be used for model continuation, see example 2.
  }
  \item{yguess }{First guess values of \code{y}, corresponding to initial grid
    \code{xguess}; a matrix with number of rows equal to the number of
    equations, and whose number of columns equals the length of \code{xguess}.
    
    if the rows of \code{yguess} have a names attribute, the names will be 
    available within the functions and used to label the output matrix.
  }
  \item{jacfunc }{jacobian (optional) - either an \R-function that evaluates the
    jacobian of \code{func} at point \code{x}, or a string with the name of a
    function or subroutine in \code{dllname} that computes the Jacobian
    (see vignette \code{"bvpSolve"} for more about this option).
    
    If \code{jacfunc} is an \R-function, it must be defined as:
    \code{jacfunc = function(x, y, parms,...)}.
    It should return the partial derivatives of \code{func} with respect to
    \code{y}, i.e. df(i,j) = dfi/dyj. See last example.
    
    If \code{jacfunc} is \code{NULL}, then a numerical approximation using
    differences is used. This is the default.
  }
  \item{bound }{boundary function (optional) - only if \code{yini} and
    \code{yend} are not available. Either an \R function that evaluates
    the i-th boundary element at point \code{x}, or a string with the name of a
    function or subroutine in \code{dllname} that computes the boundaries
    (see vignette \code{"bvpSolve"} for more about this option).

    If \code{bound} is an \R-function, it should be defined as:
    \code{bound = function(i, y, parms, ...)}. It should return the i-th
    boundary condition. See last example.
  }
  \item{jacbound }{jacobian of the boundary function (optional) - only if
    \code{bound} is defined. Either an \R function that evaluates
    the gradient of the i-th boundary element with respect to the state
    variables, at point \code{x}, or a string with the name of a
    function or subroutine in \code{dllname} that computes the boundary jacobian
    (see vignette \code{"bvpSolve"} for more about this option).

    If \code{jacbound} is an \R-function, it should be defined as:
    \code{jacbound = function(i, y, parms, ...)}.
    It should return the gradient of the i-th boundary condition.
    See last example.

    If \code{jacbound} is \code{NULL}, then a numerical approximation using
    differences is used. This is the default.
  }
  \item{leftbc }{only if \code{yini} and \code{yend} are not available and
    \code{posbound} is not specified: the number of left boundary conditions.
  }
  \item{posbound }{only used if \code{bound} is given: a vector with the position 
    (in the mesh) of the boundary conditions - only the boundary points are 
    allowed. Note that it is simpler to use \code{leftbc}. 
  }
  \item{islin }{set to \code{TRUE} if the problem is linear - this will
    speed up the simulation.
  }
  \item{nmax }{maximal number of subintervals during the calculation.
  }
  \item{order }{the order of each derivative in \code{func}. The default is that
    all derivatives are 1-st order, in which case \code{order} can be set = 
    \code{NULL}.
    
    If \code{order} is not \code{NULL}, the number of equations in \code{func} 
    must equal the length of \code{order}; the summed values of \code{order}
    must equal the number of variables (ncomp). 
    The jacobian as specified in \code{jacfunc}
    must have number of rows = number of equations and number of columns = number 
    of variables.                            
    \code{bound} and \code{jacbound} remain defined in the number of variables.
    See example 4 and 4b. 
  }   
  \item{ncomp }{used if the model is specified by compiled code, the number of
    components. See package vignette \code{"bvpSolve"}.
    
    Also to be used if the boundary conditions are specified by \code{bound}, and 
    there is no \code{yguess}
  }
  \item{atol }{error tolerance, a scalar.
  }
  \item{cond }{if \code{TRUE}, uses conditioning in the mesh selection
  }
  \item{lobatto }{if \code{TRUE}, selects a lobatto method.
  }
  \item{allpoints }{sometimes the solver estimates the solution in a number
    of extra points, and by default the solutions at these extra points will
    also be returned.

    By setting  \code{allpoints} equal to \code{FALSE},
    only output corresponding to the elements in \code{x} will be returned.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions referred to in \code{func}, \code{jacfunc},
    \code{bound} and  \code{jacbound}.  Note that ALL these subroutines
    must be defined in the shared library; it is not allowed to merge R-functions
    with compiled functions.
    
    See package vignette \code{"bvpSolve"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"bvpSolve"}.
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculated in the dll - you have
    to perform this check in the code. See deSolve's package vignette
    \code{"compiledCode"}.
  }
  \item{outnames }{only used if function is specified in compiled code and 
    \code{nout} > 0: the names of output variables calculated in the compiled 
    function. These names will be used to label the output matrix. 
    The length of \code{outnames} should be = \code{nout}.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme. This feature is included for consistency with
    the initial value problem solvers from package \code{deSolve}.

    See package vignette \code{"compiledCode"} from package \code{deSolve}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.

    See package vignette \code{"compiledCode"} from package \code{deSolve}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.

    See package vignette \code{"compiledCode"} from package \code{deSolve}.
  }
  \item{verbose }{if \code{TRUE} then more verbose output will be generated
    as "warnings".
  }
   \item{... }{additional arguments passed to the model functions.
  }
}

\value{
  A matrix of class \code{bvpSolve}, with up to as many rows as elements 
  in \code{x} and as many columns
  as elements in \code{yini} or \code{ncomp} plus the number of "global"
  values returned from \code{func}, plus an additional
  column (the first) for the \code{x}-value.

  Typically, there will be one row for each element in \code{x} unless the
  solver returns with an unrecoverable error. In certain cases, \code{bvptwp}
  will return the solution also in a number of extra points. This will occur
  if the number of points as in \code{x} was not sufficient. In order to
  not return these extra points, set \code{allpoints} equal to \code{FALSE}.

  If \code{ynames} is given, or \code{yini}, \code{yend} has a names attribute, 
  or \code{yguess} has named rows, the names will be used to label the
  columns of the output value.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
  
  Jeff Cash <j.cash@imperial.ac.uk>
  
  Francesca Mazzia <mazzia@dm.uniba.it>
}
\details{
  This is an implementation of the method \emph{twpbvpC}, written by
  Cash, Mazzia and Wright, to solve two-point boundary value problems of
  ordinary differential equations.

  A boundary value problem does not have all initial values of
  the state variable specified. Rather some conditions are specified at
  the end of the integration interval. The number of unknown boundary
  conditions must be equal to the number of equations (or the number
  of dependent variables \code{y}).

  The ODEs and boundary conditions are made available through the
  user-provided routines, \code{func} and vectors \code{yini} and \code{yend}
  or (optionally) \code{bound}.

  The corresponding partial derivatives for \code{func} and \code{bound}
  are optionally available through the
  user-provided routines, \code{jacfunc} and \code{jacbound}. Default is that
  they are automatically generated by \code{bvptwp}, using numerical differences.

  The user-requested tolerance is provided through \code{tol}. The
  default is \eqn{1e^-6}

  If the function terminates because the maximum
  number of subintervals was exceeded, then it is recommended that
  'the program be run again with a larger value for this maximum.'
  It may also help to start with a simple version of the model, and use
  its result as initial guess to solve the more complex problem
  (continuation strategy, see example 2, and package vignette "bvpSolve").

  Models may be defined in \bold{compiled C or Fortran} code, as well as
  in an  \bold{R}-function.

  This is similar to the initial value problem solvers from package
  \code{deSolve}. See package vignette \code{"bvpSolve"} for details about
  writing compiled code.

  The \bold{fcontrol} argument is a list that can supply any of the following
  components (conform the definitions in the \link{approx} function):
  \describe{
    \item{method }{specifies the interpolation method to be used.
      Choices are "linear" or "constant",}
    \item{rule }{an integer describing how interpolation is to take place
      outside the interval [min(times), max(times)].
      If \code{rule} is \code{1} then an error will be triggered and the
      calculation will stop if \code{times} extends the interval of the
      forcing function data set. If it is \code{2}, the *default*, the
      value at the closest data extreme is used, a warning will be printed if
      \code{verbose} is TRUE,

      Note that the default differs from the \code{approx} default}
    \item{f }{For method=\code{"constant"} a number between \code{0} and
      \code{1} inclusive, indicating a compromise between left- and
      right-continuous step functions. If \code{y0} and \code{y1} are the
      values to the left and right of the point then the value is
      \code{y0*(1-f)+y1*f} so that \code{f=0} is right-continuous and
      \code{f=1} is left-continuous,
      }
    \item{ties }{Handling of tied \code{times} values. Either a function
      with a single vector argument returning a single number result or the
      string "ordered".

      Note that the default is "ordered", hence the existence of ties will
      NOT be investigated; in the \code{C} code this will mean that -if ties
      exist, the first value will be used; if the dataset is not ordered,
      then nonsense will be produced.

      Alternative values for \code{ties} are \code{mean}, \code{min} etc
      }
   }
   The defaults are:

   \code{fcontrol=list(method="linear", rule = 2,  f = 0, ties = "ordered")}

   Note that only ONE specification is allowed, even if there is more than
   one forcing function data set.
   
   This -advanced- feature is explained in \code{deSolve}'s package 
   vignette "compiledCode".
}
\seealso{
  \code{\link{bvpshoot}} for the shooting method

  \code{\link{bvpcol}} for the collocation method

  \code{\link{diagnostics.bvpSolve}}, for a description of diagnostic messages
  
  \code{\link{plot.bvpSolve}}, for a description of plotting the output of the 
  BVP solvers.

}
\note{
 When \code{order} is not \code{NULL}, then it should contain the order of all
 \emph{equations} in \code{func}. If the order of some equations > 1, then 
 there will be less equations than variables. The number of equations should
 be equal to the length of \code{order}, while the number of variables will be
 the sum of order.

 For instance, if \code{order = c(1,2,3,4)}, then the first equation will be
 of order 1, the second one of order 2, ...and the last of order 4. 

 There will be 1+2+3+4 = 10 variables. For instance, if the derivative 
 function defines (A', B'', C''', D'''') respectively, then the variable vector 
 will contain values for A, B, B', C, C', C'', D, D', D'', D'''; in that 
 order. This is also the order in which the initial and end conditions of 
 all variables must be specified.

 If neq are the number of equations, and ncomp the number of variables, then 
 the Jacobian of the derivative function as specified in \code{jacfunc}
 must be of dimension (neq, ncomp).
 
 The jacobian of the boundaries, as specified in \code{jacbound} should return
 a vector of length = ncomp
}

\keyword{math}

\references{
  J.R. Cash and M.H. Wright, A deferred correction method for nonlinear
  two-point boundary value problems: implementation and numerical evaluation,
  SIAM J. Sci. Stat. Comput., 12 (1991) 971 989.

  Cash, J. R. and F. Mazzia, A new mesh selection algorithm, based on 
   conditioning, for two-point boundary value codes. J. Comput. Appl. Math. 184  
   (2005), no. 2, 362--381.

  Cash, J. R. and F. Mazzia, in press.
      Hybrid Mesh Selection Algorithms Based on Conditioning for
      Two-Point Boundary Value Problems, Journal of Numerical Analysis,
     Industrial and Applied Mathematics.
}

\examples{
## =============================================================================
## Example 1: simple standard problem
## solve the BVP ODE:
## d2y/dt^2=-3py/(p+t^2)^2
## y(t= -0.1)=-0.1/sqrt(p+0.01)
## y(t=  0.1)= 0.1/sqrt(p+0.01)
## where p = 1e-5
##
## analytical solution y(t) = t/sqrt(p + t^2).
##
## The problem is rewritten as a system of 2 ODEs:
## dy=y2
## dy2=-3p*y/(p+t^2)^2
## =============================================================================

#--------------------------------
# Derivative function
#--------------------------------
fun <- function(t, y, pars) {
  dy1 <- y[2]
  dy2 <- - 3*p*y[1] / (p+t*t)^2
  return(list(c(dy1,
                dy2))) }


# parameter value
p    <- 1e-5

# initial and final condition; second conditions unknown
init <- c(y = -0.1 / sqrt(p+0.01), dy=NA)
end  <- c(     0.1 / sqrt(p+0.01), NA)

# Solve bvp
sol  <- as.data.frame(bvptwp(yini = init, x = seq(-0.1, 0.1, by = 0.001),
        func = fun, yend = end))
plot(sol$x, sol$y, type = "l")

# add analytical solution
curve(x/sqrt(p+x*x), add = TRUE, type = "p")


## =============================================================================
## Example 1b: 
## Same problem, now solved as a second-order equation.
## =============================================================================

fun2 <- function(t, y, pars)  {
  dy <- - 3 * p * y[1] / (p+t*t)^2
  list(dy)
}
sol2  <- bvptwp(yini = init, yend = end, order = 2, 
                x = seq(-0.1, 0.1, by = 0.001), func = fun2)

## =============================================================================
## Example 2: simple
## solve d2y/dx2 + 1/x*dy/dx + (1-1/(4x^2)y = sqrt(x)*cos(x),
## on the interval [1,6] and with boundary conditions:
## y(1)=1, y(6)=-0.5
##
## Write as set of 2 odes
## dy/dx = y2
## dy2/dx  = - 1/x*dy/dx - (1-1/(4x^2)y + sqrt(x)*cos(x)
## =============================================================================

f2 <- function(x, y, parms) {
  dy  <- y[2]
  dy2 <- -1/x*y[2] - (1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
  list(c(dy, dy2))
}

x    <- seq(1, 6, 0.1)
sol  <- bvptwp(yini = c(y = 1, dy = NA),
               yend = c(-0.5, NA), x = x, func = f2)
plot(sol, which = "y")

# add the analytic solution
curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x), add = TRUE, type = "l")

## =============================================================================
## Example 3  - solved with automatic continuation
## d2y/dx2 = y/xi
## =============================================================================

Prob1 <- function(t, y, xi)
   list(c( y[2] , y[1]/xi ))

x <-  seq(0, 1, by = 0.01)
xi <- 0.1
twp <- bvptwp(yini = c(1, NA), yend = c(0, NA), x = x,
             islin = TRUE, func = Prob1, parms = xi, eps = xi)

xi <-0.00001
twp2 <- bvptwp(yini = c(1, NA), yend = c(0, NA), x = x,
               islin = TRUE, func = Prob1, parms = xi, eps = xi)

plot(twp, twp2, which = 1, main = "test problem 1")

# exact solution
curve(exp(-x/sqrt(xi))-exp((x-2)/sqrt(xi))/(1-exp(-2/sqrt(xi))),
      0, 1, add = TRUE, type = "p")

curve(exp(-x/sqrt(0.1))-exp((x-2)/sqrt(0.1))/(1-exp(-2/sqrt(0.1))),
      0, 1, add = TRUE, type = "p")

## =============================================================================
## Example 4  - solved with specification of boundary, and jacobians
## d4y/dx4 =R(dy/dx*d2y/dx2 -y*dy3/dx3)
## y(0)=y'(0)=0
## y(1)=1, y'(1)=0
##
## dy/dx  = y2
## dy2/dx = y3    (=d2y/dx2)
## dy3/dx = y4    (=d3y/dx3)
## dy4/dx = R*(y2*y3 -y*y4)
## =============================================================================

f2<- function(x, y, parms, R) {
  list(c(y[2], y[3], y[4], R*(y[2]*y[3] - y[1]*y[4]) ))
}

df2 <- function(x, y, parms, R) {
  matrix(nrow = 4, ncol = 4, byrow = TRUE, data = c(
             0,        1,     0,     0,
             0,        0,     1,     0,
             0,        0,     0,     1,
             -1*R*y[4],R*y[3],R*y[2],-R*y[1]))
}

g2 <- function(i, y, parms, R) {
  if (i == 1) return(y[1])
  if (i == 2) return(y[2])
  if (i == 3) return(y[1]-1)
  if (i == 4) return(y[2])
}

dg2 <- function(i, y, parms, R) {
  if (i == 1) return(c(1, 0, 0, 0))
  if (i == 2) return(c(0, 1, 0, 0))
  if (i == 3) return(c(1, 0, 0, 0))
  if (i == 4) return(c(0, 1, 0, 0))
}

init <- c(1, NA)
R    <- 100
sol  <- bvptwp(x = seq(0, 1, by = 0.01), leftbc = 2,
          func = f2, R = R, ncomp = 4,
          bound = g2, jacfunc = df2, jacbound = dg2)
plot(sol[,1:2])  # columns do not have names

mf <- par ("mfrow")
sol  <- bvptwp(x = seq(0, 1, by = 0.01), leftbc = 2,
          func = f2, ynames = c("y", "dy", "d2y", "d3y"), R=R,
          bound = g2, jacfunc = df2, jacbound = dg2)
plot(sol)        # here they do
par(mfrow = mf)

## =============================================================================
## Example 4b - solved with specification of boundary, and jacobians
## and as a higher-order derivative
## d4y/dx4 =R(dy/dx*d2y/dx2 -y*dy3/dx3)
## y(0)=y'(0)=0
## y(1)=1, y'(1)=0
## =============================================================================

# derivative function: one fourth-order derivative
f4th <- function(x, y, parms, R) {
  list(R * (y[2]*y[3] - y[1]*y[4]))
}

# jacobian of derivative function
df4th <- function(x, y, parms, R)  {
 df <- matrix(nrow = 1, ncol = 4, byrow = TRUE, data = c(
             -1*R*y[4], R*y[3], R*y[2], -R*y[1]))
}

# boundary function - same as previous example

# jacobian of boundary - same as previous

# order = 4 specifies the equation to be 4th order
sol2 <- bvptwp(x = seq(0, 1, by = 0.01), 
          ynames = c("y", "dy", "d2y", "d3y"),
          posbound = c(0, 0, 1, 1), func = f4th, R = R, order = 4, 
          bound = g2, jacfunc = df4th, jacbound = dg2)

max(abs(sol2-sol))

}
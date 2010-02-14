\name{bvptwp}
\alias{bvptwp}
\title{
  Solver for two-point boundary value problems of ordinary differential
  equations, using a mono-implicit Runge-Kutta formula
}

\description{
  Solves a boundary value problem of a system of ordinary differential
  equations. This is an implementation of the fortran
  code twpbvpc, based on mono-implicit Runge-Kutta formulae
  of orders 4, 6 and 8 in a deferred correction framework and 
  that uses conditioning in the mesh selection.

  written by J.R. Cash, F. Mazzia and M.H. Wright.
}

\usage{
bvptwp(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL,
  xguess=NULL, yguess=NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
  leftbc=NULL, islin=FALSE, nmax=1000, atol=1e-8, cond = FALSE,
  allpoints=TRUE, dllname=NULL, initfunc=dllname, ncomp=NULL,
  forcings=NULL, initforc = NULL, fcontrol=NULL, verbose=FALSE, ...)
}
\arguments{
  \item{yini }{either a vector with the initial (state) variable values for
    the ODE system, or \code{NULL}.

    If \code{yini} is a vector, use \code{NA} for an
    initial value which is not available. 
    
    If \code{yini} has a
    \code{names} attribute, the names will be used to label the output matrix.
    
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
    for a final value which is not available.
    
    If \code{yend} has a  \code{names} attribute, and \code{yini} does not,
    the names will be available within the functions and 
    used to label the output matrix.

    If \code{yend = NULL}, then the boundary conditions must be specified
    via function \code{bound}; 
    if not \code{NULL} then \code{yini} should also be not \code{NULL}.
  }
  \item{parms }{vector or a list with parameters passed to \code{func},
    \code{jacfunc}, \code{bound} and \code{jacbound} (if present).
  }
  \item{guess }{guess for the value(s) of the unknown \emph{initial} conditions;
    i.e. one value for each \code{NA} in \code{yini}. Ignored if \code{xguess}
    and \code{yguess} are given.
    
    The length of \code{guess} should equal the number of \code{NA}s
    in \code{yini}. If not provided, a value = 0 is assumed for each
    \code{NA} in \code{yini} and a warning printed.

    if \code{guess} has a names attributes, the names will be 
    available within the functions and used to label the output matrix.
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
    
    if the rows of \code{guess} have a names attribute, the names will be 
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
  \item{leftbc }{only if \code{yini} and \code{yend} are not available: the
    number of left boundary conditions.
  }
  \item{islin }{set to \code{TRUE} if the problem is linear - this will
    speed up the simulation.
  }
  \item{nmax }{maximal number of subintervals during the calculation.
  }
  \item{atol }{error tolerance, a scalar.
  }
  \item{cond }{if \code{TRUE}, uses conditioning in the mesh selection
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
  \item{ncomp }{only used if the model is specified by compiled code, the number of
    components (or equations). See package vignette \code{"bvpSolve"}.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme. This feature is included for consistency with
    the intial value problem solvers from packatge \code{deSolve}.

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
  A matrix with up to as many rows as elements in \code{x} and as many columns
  as elements in \code{yini} or \code{ncomp} plus the number of "global"
  values returned from \code{func}, plus an additional
  column (the first) for the \code{x}-value.

  Typically, there will be one row for each element in \code{x} unless the
  solver returns with an unrecoverable error. In certain cases, \code{bvptwp}
  will return the solution also in a number of extra points. This will occur
  if the number of points as in \code{x} was not sufficient. In order to
  not return these extra points, set \code{allpoints} equal to \code{FALSE}.

  If \code{yini} has a names attribute, it will be used to label the columns
  of the output value. If \code{yini} is not named, the solver will try to find
  the names in \code{yend}, \code{yguess} and ultimately \code{guess}.
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
  
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
}
\seealso{
  \code{\link{bvpshoot}} for the shooting method
}

\keyword{math}

\references{
  J.R. Cash and M.H. Wright, A deferred correction method for nonlinear
  two-point boundary value problems: implementation and numerical evaluation,
  SIAM J. Sci. Stat. Comput., 12 (1991) 971 989.

  Cash, J. R. and F. Mazzia, A new mesh selection algorithm, based on 
   conditioning, for two-point boundary value codes. J. Comput. Appl. Math. 184  
   (2005), no. 2, 362--381.

  \url{http://www.ma.ic.ac.uk/~jcash/BVP\_software}
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
fun <- function(t,y,pars)
{ dy1 <-y[2]
  dy2 <- - 3*p*y[1]/(p+t*t)^2
  return(list(c(dy1,
             dy2))) }


# parameter value
p    <-1e-5

# initial and final condition; second conditions unknown
init <- c(y=-0.1/sqrt(p+0.01), dy=NA)
end  <- c(0.1/sqrt(p+0.01), NA)

# Solve bvp
sol  <- as.data.frame(bvptwp(yini=init,x=seq(-0.1,0.1,by=0.001),
        func=fun, yend=end, guess=1))
plot(sol$x,sol$y,type="l")

# add analytical solution
curve(x/sqrt(p+x*x),add=TRUE,type="p")

## =============================================================================
## Example 1b: simple
## solve d2y/dx2 + 1/x*dy/dx + (1-1/(4x^2)y = sqrt(x)*cos(x),
## on the interval [1,6] and with boundary conditions:
## y(1)=1, y(6)=-0.5
##
## Write as set of 2 odes
## dy/dx = y2
## dy2/dx  = - 1/x*dy/dx - (1-1/(4x^2)y + sqrt(x)*cos(x)
## =============================================================================

f2 <- function(x,y,parms)
{
 dy  <- y[2]
 dy2 <- -1/x*y[2]-(1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy,dy2))
}

x    <- seq(1,6,0.1)
sol  <- bvptwp(yini=c(y=1,dy=NA),yend=c(-0.5,NA),x=x,func=f2,guess=1)
plot(sol,which="y")

# add the analytic solution
curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x),add=TRUE,type="l")

## =============================================================================
## Example 2  - using continuation
## Problem nr 24 from http://www.ma.ic.ac.uk/~jcash/BVP_software/PROBLEMS.PDF
## xi*A*y*y'' -((1+ga)/2-xi*A'*y*y'+y'/y+(A'/A*(1-(ga-1)/2*y^2)=0
## where A = 1 + x2, y(x=0) = 0:9129, y(x=1) = 0:375
## =============================================================================

Prob24 <- function(x, y, pars, xi) {
  A <- 1+x*x
  AA <- 2*x
  ga <- 1.4
  list(c(y[2],
       (((1+ga)/2 -xi*AA)*y[1]*y[2]-y[2]/y[1]-
       (AA/A)*(1-(ga-1)*y[1]^2/2))/(xi*A*y[1])  ))
}

ini <- c(0.9129,NA)
end <- c(0.375,NA)

# This problem needs an initial guess for y; (as it divides by y,
# the default, 0, does not work)

# guess of y; set equal to 1
xguess <- c(0,1)
yguess <- matrix(nrow=2,ncol=2,data=1)

# start with relatively large value of xi
xi <-0.01
print(system.time(
mod1  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24,xi=xi, xguess=xguess,yguess=yguess)
))        
plot(mod1[,1:2], ylim = c(0.4,1.4), type="l", lwd=2)

# use previous solution to initialise solver for smaller xi
xi <-0.005
mod1b  <- bvptwp(yini=ini,yend=end,x=seq(0,1,by=0.01),
        func=Prob24,xi=xi, xguess=mod1[,1],yguess=t(mod1[,2:3]))
lines(mod1b, lwd=2, col="red")


## =============================================================================
## Example 3  - solved with specification of boundary, and jacobians
## d4y/dx4 =R(dy/dx*d2y/dx2 -y*dy3/dx3)
## y(0)=y'(0)=0
## y(1)=1, y'(1)=0
##
## dy/dx  = y2
## dy2/dx = y3    (=d2y/dx2)
## dy3/dx = y4    (=d3y/dx3)
## dy4/dx = R*(y2*y3 -y*y4)
## =============================================================================

f2<- function(x,y,parms,R)
{
  list(c(y[2],y[3],y[4],R*(y[2]*y[3]-y[1]*y[4]) ))
}

df2 <- function(x,y,parms,R)
{
 df <- matrix(nr=4,nc=4,byrow=TRUE,data=c(
             0,1,0,0,
             0,0,1,0,
             0,0,0,1,
             -1*R*y[4],R*y[3],R*y[2],-R*y[1]))
}

g2 <- function(i,y,parms,R)
{
  if ( i ==1) return(y[1])
  if ( i ==2) return(y[2])
  if ( i ==3) return(y[1]-1)
  if ( i ==4) return(y[2])
}

dg2 <- function(i,y,parms,R)
{
  if ( i ==1) return(c(1,0,0,0))
  if ( i ==2) return(c(0,1,0,0))
  if ( i ==3) return(c(1,0,0,0))
  if ( i ==4) return(c(0,1,0,0))
}

require(bvpSolve)
init <- c(1,NA)
R    <- 0.01
sol  <- bvptwp(x=seq(0,1,by=0.01), leftbc=2,
        func=f2, guess=1,R=R,
        bound=g2, jacfunc=df2, jacbound=dg2)
plot(sol[,1:2])  # columns do not have names

mf <- par ("mfrow")
sol  <- bvptwp(x=seq(0,1,by=0.01), leftbc=2,
        func=f2, guess=c(y=1,dy=1,dy2=1,dy3=1),R=R,
        bound=g2, jacfunc=df2, jacbound=dg2)
plot(sol)        # here they do
par(mfrow=mf)

}


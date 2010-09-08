\name{bvpshoot}
\alias{bvpshoot}
\title{
  Solver for two-point boundary value problems of ordinary differential
  equations, using the single shooting method
}

\description{
  Solves a boundary value problem of a system of ordinary differential
  equations using the single shooting method. 
  This combines the integration routines from package 
  \code{deSolve} with root-finding methods from package \code{rootSolve}.
  
  Preferentially \code{bvptwp} or \code{bvpcol} should be used rather than
  \code{bvpshoot}, as they give more precise output.
}

\usage{
bvpshoot(yini = NULL, x, func, yend = NULL, parms = NULL, 
    order = NULL, guess = NULL,
    jacfunc = NULL, bound = NULL, jacbound = NULL, 
    leftbc = NULL, posbound = NULL, ncomp = NULL, 
    atol = 1e-8, rtol = 1e-8, extra = NULL, 
    maxiter = 100, positive = FALSE, method = "lsoda",...)
}

\arguments{
  \item{yini }{either a \emph{vector} with the initial (state) variable
    values for the ODE system, or a \emph{function} that calculates the
    initial condition, or \code{NULL}.

    If \code{yini} is a function, it should be defined as:
    \code{yini <- function(y, parms,...)}; where \code{y} are the initial values,
    and \code{parms} the parameters.

    if \code{yini} is a vector then use \code{NA}
    for an initial value which is not available. 
    
    If \code{yini} has a  \code{names} attribute, the names will be 
    available within the functions and used to label the output matrix.
    
    if \code{yini = NULL} then \code{bound} should be specified;
    if not \code{NULL} then \code{yend} should also be not \code{NULL}
  }
  \item{x }{ sequence of the independent variable for which output is wanted;
    the first value of \code{x} must be the initial value (at which
    \code{yini} is defined), the final value the end condition (at which
    \code{yend} is defined).
  }
  \item{func }{an \R-function that computes the values of the derivatives in
    the ODE system (the model definition) at x. \code{func} must be defined as:
    \code{func  = function(x, y, parms, ...)}.  \code{x} is the current point of
    the independent variable in the integration, \code{y} is the current
    estimate of the (state) variables in the ODE system.  If the initial
    values \code{yini} or \code{yend} has a names attribute, the names will 
    be available inside \code{func}.  \code{parms} is a vector or list of parameters;
    ... (optional) are any other arguments passed to the function.

    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{x}, and whose next elements are global values that are required at
    each point in \code{x}.
    
    Note that it is not possible to use \code{bvpshoot} with functions
    defined in compiled code. Use \link{bvptwp} instead.
  }
  \item{yend }{either a vector with the final (state) variable values for the
    ODE system, a \emph{function} that calculates the final condition 
    or \code{NULL}; 
    
    if \code{yend} is a vector use \code{NA}
    for a final value which is not available.

    If \code{yend} is a function, it should be defined as:
    \code{yend <- function (y, yini, parms, ...)}; where \code{y} are the 
    final values, \code{yini} the initial values and \code{parms} the 
    parameters.
    
    If \code{yend} has a  \code{names} attribute, and \code{yini} does not,
    the names will be available within the functions and 
    used to label the output matrix.

    if \code{yend = NULL} then \code{bound} should be specified; 
    if not \code{NULL} then \code{yini} should also be not \code{NULL}. 
  }
  \item{parms }{vector or a list with parameters passed to \code{func},
    \code{jacfunc}, \code{bound} and \code{jacbound} (if present).
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
    See examples. 
  }   
  \item{guess }{guess for the value(s) of the unknown \emph{initial} conditions;
    
    if initial and final conditions are specified by \code{yini} and 
    \code{yend}, then \code{guess} should contain one 
    value for each \code{NA} in \code{yini}.
    The length of \code{guess} should thus equal the number of \emph{unknown initial 
    conditions} (=\code{NA}s in \code{yini}). If \code{guess} is not provided, 
    a value = 0 is assumed for each \code{NA} in \code{yini} and a warning 
    printed.
    
    If initial and final conditions are specified by the boundary function 
    \code{bound}, then \code{guess} should contain the initial guess for 
    \emph{all initial conditions}, i.e. its length should equal the number of state
    variables in the ODE system;
    if in this case \code{guess} has a names attribute, the names will be 
    available within the functions and used to label the output matrix.
    If \code{guess} is not provided, then \code{ncomp} should specify the 
    total number of variables, a value = 0 will be assumed for the initial
    conditions and a warning printed.
  }
  \item{jacfunc }{jacobian (optional) - an \R-function that evaluates the
    jacobian of \code{func} at point \code{x}.
    
    \code{jacfunc} must be defined as \code{jacfunc = function(x, y, parms,...)}.
    It should return the partial derivatives of \code{func} with respect to
    \code{y}, i.e. df(i,j) = dfi/dyj.  
    
    If \code{jacfunc} is \code{NULL}, then a numerical approximation using
    differences is used. This is the default.
    
    \code{jacfunc} is passed to the initial value problem solver.
  }
  \item{bound }{boundary function (optional) - only if \code{yini} and
    \code{yend} are not available. An \R function that evaluates
    the i-th boundary element at point \code{x}.

    \code{bound} should be defined as:
    \code{bound = function(i, y, parms, ...)}. It should return the i-th
    boundary condition.  
    
    if not \code{NULL}, \code{bound} defines the root to be solved by the 
    root solving algorithm.
  }
  \item{jacbound }{jacobian of the boundary function (optional) - only if
    \code{bound} is defined. An \R function that evaluates
    the gradient of the i-th boundary element with respect to the state
    variables, at point \code{x}.

    \code{jacbound} should be defined as:
    \code{jacbound = function(i, y, parms, ...)}.
    It should return the gradient of the i-th boundary condition.
    See last example.

    \code{jacbound} is passed to the root solver.
  }
  \item{leftbc }{only if \code{yini} and \code{yend} are not available: the
    number of left boundary conditions.
  }
  \item{posbound }{only used if \code{bound} is given: a vector with the position 
    (in the mesh) of the boundary conditions - only points that are in \code{x}
    are allowed. Note that, if the boundary conditions are at the ends of the 
    integration interval, it is simpler to use \code{leftbc}. 
  }
  \item{ncomp }{only used if the boundaries are specified via the boundary 
    function \code{bound} and \code{guess} is not specified. 
    The number of components.  
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one
    value for each unknown element - passed to function
    \link[rootSolve]{multiroot} - see help of this function.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector, one
    value for each unknown element - passed to function
    \link[rootSolve]{multiroot} - see help of this function.
  }
  \item{extra }{if too many boundary conditions are given, then it is assumed
    that an extra parameter has to be estimated.

    \code{extra} should contain the initial guess of this extra parameter.
  }
  \item{maxiter }{the maximal number of iterations allowed in the root solver.
  }
  \item{positive }{set to \code{TRUE} if dependent variables (y) have to be
    positive numbers.
  }
  \item{method }{the integration method used, one of ("lsoda", "lsode",
    "lsodes", "vode", "euler", "rk4", "ode23" or "ode45").
  }
  \item{... }{additional arguments passed to the integrator and (possibly)
    the model functions.
  }
}
\value{
  A matrix with up to as many rows as elements in \code{x} and as many columns
  as the number of state variables in the ODE system plus the number of 
  "global" values returned in the next elements of the return from \code{func}, 
  plus an additional column (the first) for the x-value.

  There will be one row for each element in \code{x} unless the solver
  returns with an unrecoverable error.

  If \code{yini} has a names attribute, it will be used to label the columns
  of the output value. If \code{yini} is not named, the solver will try to find
  the names in \code{yend}.  If the boundaries are specified by \code{bound} 
  then the names from \code{guess} will be used.

  The output will have the attribute \code{roots}, which returns the value(s)
  of the root(s) solved for (\code{root}), the function value (\code{f.root}),
  and the number of iterations (\code{iter}) required to find the root.
}

\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\details{
  This is a simple implementation of the shooting method to solve boundary
  value problems of ordinary differential equations.

  A boundary value problem does not have all initial values of
  the state variable specified. Rather some conditions are specified at
  the end of the integration interval.

  The shooting method, is a root-solving method.
  There are two strategies:
  
  \describe{
    \item{yini and yend specified }{
     If initial and end conditions are specified with \code{yini} and \code{yend}
     then the (unspecified) initial conditions are the unknown values to be solved 
     for; the function value whose root has to be found are the deviations from 
     the specified conditions at the end of the integration interval. 

     Thus, starting with an initial guess of the initial conditions (as
     provided in \code{guess}), the ODE model is solved as an initial value
     problem, and after termination, the discrepancy of the modeled final
     conditions with the known final condition is assessed (the cost function).
     The root of this cost function is to be found.
    }
    \item{bound specified }{If initial and end conditions are specified with 
      \code{bound}, then the 
      unknowns are all initial conditions; the function whose root is to be 
      found is given by \code{bound}.
    }
  }
  
  Starting from a guess of the initial values, one of the integrators from 
  package \code{deSolve} (as specified with \code{method}) is used to solve the 
  resulting initial value problem.

  Function \code{multiroot} from package \code{rootSolve} is used to
  retrieve the root.

  For this method to work, the model should be even determined, i.e. the number
  of equations should equal the number of unknowns.

  \code{bvpshoot} distinguishes two cases:
  
  1. the total number of specified boundary conditions (on both the start and
  end of the integration interval) equals the number of boundary value
  problem equations (or the number of dependent variables \code{y}).
  
  2. The number of boundary conditions specified \emph{exceeds}
  the number of equations. In this case, \code{extra} parameters have to be solved
  for to make the model even determined.

  See example nr 4.
}

\seealso{
  \code{\link{bvptwp}} for the MIRK method

  \code{\link[deSolve]{lsoda}}, \code{\link[deSolve]{lsode}},
    \code{\link[deSolve]{lsodes}}, \code{\link[deSolve]{vode}},

  \code{\link[deSolve]{rk}}, \code{\link[deSolve]{rkMethod}}
    for details about the integration method
    
  \code{\link[rootSolve]{multiroot}}, the root-solving method used

  \code{\link{diagnostics.bvpSolve}}, for a description of diagnostic messages
  
  \code{\link{plot.bvpSolve}}, for a description of plotting the output of the 
  BVP solvers.
}
\keyword{math}

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
fun <- function(t, y, pars)
{ dy1 <- y[2]
  dy2 <- - 3*p*y[1] / (p+t*t)^2
  return(list(c(dy1,
               dy2))) }


# parameter value
p    <- 1e-5

# initial and final condition; second conditions unknown
init <- c(y = -0.1 / sqrt(p+0.01), dy = NA)
end  <- c(     0.1 / sqrt(p+0.01), NA)

# Solve bvp
sol  <- bvpshoot(yini = init, x = seq(-0.1, 0.1, by = 0.001),
          func = fun, yend = end, guess = 1)
          
plot(sol, which = "y", type = "l")

# add analytical solution
curve(x/sqrt(p+x*x), add = TRUE, type = "p")

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

f2 <- function(x, y, parms)
{
 dy  <- y[2]
 dy2 <- -1/x*y[2] - (1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy, dy2))
}

x    <- seq(1, 6, 0.1)
sol  <- bvpshoot(yini = c(y = 1, dy = NA), yend = c(-0.5, NA), 
          x = x, func = f2, guess = 1)

# plot both state variables
plot(sol, type = "l", lwd = 2)

# plot only y and add the analytic solution
plot(sol, which = "y")

curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x), add = TRUE, type = "l")


## =============================================================================
## Example 2 - initial condition is a function of the unknown x
## tubular reactor with axial dispersion
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## on the interval [0,1] and with initial conditions:
## y'0=Pe(y(0)-1), y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
## =============================================================================

Reactor<-function(x, y, parms)
{
  list(c(y[2],
         Pe * (y[2]+R*(y[1]^n))))
}

Pe <- 1
R  <- 2
n  <- 2
x  <- seq(0, 1, by = 0.01)

# 1. yini is a function here
yini <- function (x, parms) c(x, Pe*(x-1))

system.time(
  sol <- bvpshoot(func = Reactor, yend = c(y = NA, dy = 0), 
          yini = yini, x = x, extra = 1)
)
plot(sol, which = "y", main = "Reactor", type = "l", lwd = 2)
attributes(sol)$roots

# 2. using boundary function rather than yini...
bound <- function(i, y, p) {
  if (i == 1) return(y[2] - Pe*(y[1]-1))
  if (i == 2) return(y[2])
}

# need to give number of left boundary conditions + guess of all initial 
# conditions (+ names)
system.time(
Sol2<- bvpshoot(func = Reactor, x = x, leftbc = 1,
            bound = bound, guess = c(y = 1, dy = 1) )
)
attributes(Sol2)$roots


# 3. boundary function with jacobian of boundary function
jacbound <- function(i, y, p) {
  if (i == 1) return(c(-Pe*y[1], 1))
  if (i == 2) return(c(0, 1))
}
system.time(
Sol3<-bvpshoot(func = Reactor, x = x, leftbc = 1, bound = bound, 
           jacbound = jacbound, guess = c(y = 1, dy = 1) )
)
attributes(Sol3)$roots

## =============================================================================
## Example 2b - same as 2 but written as higher-order equation
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## on the interval [0,1] and with initial conditions:
## y'0=Pe(y(0)-1), y'(1)=0
## =============================================================================

Reactor2<-function(x, y, parms) 
  list (Pe * (y[2]+R*(y[1]^n)))

Pe <- 1
R  <- 2
n  <- 2
x  <- seq(0, 1, by = 0.01)

# 1. yini is a function here
yini <- function (x, parms) c(x, Pe*(x-1))

# need to specify that order = 2
system.time(
  sol2 <- bvpshoot(func = Reactor2, yend = c(y = NA, dy = 0), order=2,
          yini = yini, x = x, extra = 1)
)
max(abs(sol2-sol))

## =============================================================================
## Example 3 - final condition is a residual function
## The example MUSN from Ascher et al., 1995.
## Numerical Solution of Boundary Value Problems for Ordinary Differential
## Equations, SIAM, Philadelphia, PA, 1995.
##
## The problem is
##      u' =  0.5*u*(w - u)/v
##      v' = -0.5*(w - u)
##      w' = (0.9 - 1000*(w - y) - 0.5*w*(w - u))/z
##      z' =  0.5*(w - u)
##      y' = -100*(y - w)
##   on the interval [0 1] and with boundary conditions:
##      u(0) = v(0) = w(0) = 1,  z(0) = -10,  w(1) = y(1)
## =============================================================================

musn <- function(t, Y, pars)
{
  with (as.list(Y),
  {
   du <- 0.5 * u * (w-u)/v
   dv <- -0.5 * (w-u)
   dw <- (0.9 - 1000 * (w-y) - 0.5 * w * (w-u))/z
   dz <- 0.5 * (w-u)
   dy <- -100 * (y-w)
   return(list(c(du, dv, dw, dy, dz)))
  })
}

#--------------------------------
# Residuals of end conditions
#--------------------------------
res  <- function (Y, yini, pars)  with (as.list(Y), w-y)

#--------------------------------
# Initial values; y= NA= not available
#--------------------------------

init <- c(u = 1, v = 1, w = 1, y = NA, z = -10)
sol  <-bvpshoot(y = init, x = seq(0, 1, by = 0.05), func = musn,
           yend = res, guess = 1, atol = 1e-10, rtol = 0)
pairs(sol, main = "MUSN")

## =============================================================================
## Example 4 - solve also for unknown parameter
## Find the 4th eigenvalue of Mathieu's equation:
## y''+(lam-10cos2t)y=0   on the interval [0,pi]
## y(0)=1, y'(0)=0  and y'(pi)=0
## The 2nd order problem is rewritten as 2 first-order problems:
## dy=y2
## dy2= -(lam-10cos(2t))*y
## =============================================================================

mathieu<- function(t, y, lam)
{
 list(c(y[2], -(lam - 10 * cos(2*t)) * y[1]))
}

yini <- c(1, 0)   # initial condition(y1=1,dy=y2=0)
yend <- c(NA, 0)  # final condition at pi (y1=NA, dy=0)

# there is one extra parameter to be fitted: "lam"; its initial guess = 15
Sol  <- bvpshoot(yini = yini, yend = yend, x = seq(0, pi, by = 0.01),
          func = mathieu, guess = NULL, extra = 15)
plot(Sol)
attr(Sol, "roots")  # root gives the value of "lam" (17.1068)
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

 Do not specify the jacobian if the maximal order>1.
}


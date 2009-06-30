\name{bvpshoot}
\alias{bvpshoot}
\title{
  Solver for boundary value problems of ordinary differential equations
}

\description{
  Solves a boundary value problem of a system of ordinary differential
  equations using the shooting method.
}

\usage{
bvpshoot(yini,x,func,yend,parms=NULL,
  guess=NULL, extra=NULL,atol=1e-8,rtol=1e-8,
  maxiter=100, positive = FALSE, method="lsoda",...)
}

\arguments{
  \item{yini }{either a vector with the initial (state) values for the ODE
    system, or a function that calculates the initial condition.

    If \code{yini} has a name attribute, the names will be used to label the
    output matrix. If \code{yini} is a function, it will be called as:
    \code{yini(y,parms,...)};

    if \code{yini} is a vector then use \code{NA}
    for an initial value which is not available.
  }
  \item{x }{ sequence of the independent variable for which output is wanted;
    the first value of \code{x} must be the initial value (at which
    \code{yini} is defined), the final value the end condition (at which
    \code{yend} is defined).
  }
  \item{func }{an \R-function that computes the values of the derivatives in
    the ODE system (the model definition) at x. \code{func} must be defined as:
    \code{yprime = func(x, y, parms,...)}.  \code{x} is the current point of
    the independent variable in the integration, \code{y} is the current
    estimate of the (state) variables in the ODE system.  If the initial
    values \code{yini} has a names attribute, the names will be available
    inside \code{func}.  \code{parms} is a vector or list of parameters;
    ... (optional) are any other arguments passed to the function.

    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{x}, and whose next elements are global values that are required at
    each point in \code{x}.
  }
  \item{yend }{either a vector with the final (state) values for the
    ODE system, or \code{NULL}; if \code{yend} is a vector use \code{NA}
    for a final value which is not available.
  }
  \item{parms }{parameters passed to \code{func}.
  }
  \item{guess }{guess for the value(s) of the unknown initial conditions;
    i.e. one value for each \code{NA} in \code{yini}.
    The length of \code{guess} should thus equal the number of NAs
    in \code{yini}. If not provided, a value = 0 is assumed for each
    \code{NA} and a warning printed.
  }
  \item{extra }{if too many boundary conditions are given, then an extra
    parameter can be estimated. \code{extra} should contain the initial guess.
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one
    value for each unknown element - passed to function
    \link[rootSolve]{multiroot} - see help of this function.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector, one
    value for each unknown element - passed to function
    \link[rootSolve]{multiroot} - see help of this function.
  }
  \item{maxiter }{the maximal number of iterations allowed in the root solver.
  }
  \item{positive }{set to \code{TRUE} if state variables have to be
    positive numbers.
  }
  \item{method }{the integration method used, one of ("lsoda", "lsode",
    "lsodes", "vode", "euler", "rk4", "ode23" or "ode45").
  }
  \item{... }{additional arguments passed to the integrator (and possibly
    the model functions).
  }
}
\value{
  A matrix with up to as many rows as elements in times and as many columns
  as elements in \code{yini} plus the number of "global" values returned
  in the second element of the return from \code{func}, plus an additional
  column (the first) for the x-value.

  There will be one row for each element in \code{x} unless the solver
  returns with an unrecoverable error.

  If \code{y} has a names attribute, it will be used to label the columns
  of the output value.

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

  The shooting method, is a root-solving method, where the unknown initial
  conditions are the unknown values to be solved for; the function value
  whose root has to be found are the deviations from the final conditions.

  Thus, starting with an initial guess of the initial conditions (as
  provided in \code{guess}, the ODE model is solved (as an initial value
  problem), and after termination, the discrepancy of the modeled final
  conditions with the known final condition is assessed (the cost function).
  The root of this cost function is to be found.

  In \code{bvpshoot} one of the integrators from package \code{deSolve} (as
  specified with \code{method}) are used to solve the resulting initial
  value problem.

  Function \code{multiroot} from package \code{rootSolve} is used to
  retrieve the root.

  For this method to work, the model should be even determined, i.e.
  the total number of specified boundary conditions (on both the start and
  end of the integration interval) should equal the number of boundary value
  problem equations.
  An exception is when the number of boundary conditions specified exceeds
  the number of equations. In this case, extra parameters have to be solved
  for to make the model even determined.

  See examples
}

\seealso{
  \code{\link{bvptwp}} for the MIRK method

  \code{\link[deSolve]{lsoda}}, \code{\link[deSolve]{lsode}},
    \code{\link[deSolve]{lsodes}}, \code{\link[deSolve]{vode}},

  \code{\link[deSolve]{rk}}, \code{\link[deSolve]{rkMethod}}
    for details about the integration method}
\keyword{math}

\examples{
#################################################################
# Example 1: simple standard problem
# solve the BVP ODE:
# d2y/dt^2=-3py/(p+t^2)^2
# y(t= -0.1)=-0.1/sqrt(p+0.01)
# y(t=  0.1)= 0.1/sqrt(p+0.01)
# where p = 1e-5
#
# analytical solution y(t) = t/sqrt(p + t^2).
#
# The problem is rewritten as a system of 2 ODEs:
# dy=y2
# dy2=-3p*y/(p+t^2)^2
################################################################################

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
init <- c(-0.1/sqrt(p+0.01), NA)
end  <- c(0.1/sqrt(p+0.01), NA)

# Solve bvp
sol  <- as.data.frame(bvpshoot(yini=init, x=seq(-0.1,0.1,by=0.001),
        func=fun, yend=end, guess=1))
plot(sol$time,sol[,2],type="l")

# add analytical solution
curve(x/sqrt(p+x*x),add=TRUE,type="p")

#################################################################
# Example 1b: simple
# solve d2y/dx2 + 1/x*dy/dx + (1-1/(4x^2)y = sqrt(x)*cos(x),
# on the interval [1,6] and with boundary conditions:
# y(1)=1, y(6)=-0.5
#
# Write as set of 2 odes
# dy/dx = y2
# dy2/dx  = - 1/x*dy/dx - (1-1/(4x^2)y + sqrt(x)*cos(x)
#################################################################

f2 <- function(x,y,parms)
{
 dy  <- y[2]
 dy2 <- -1/x*y[2]-(1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy,dy2))
}

x    <- seq(1,6,0.1)
sol  <- bvpshoot(yini=c(1,NA), yend=c(-0.5,NA), x=x, func=f2, guess=1)
plot(sol)

# add the analytic solution
curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x),add=TRUE,type="l")


#################################################################
# Example 2 - initial condition is a function of the unknown x
# tubular reactor with axial dispersion
# y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
# on the interval [0,1] and with initial conditions:
# y'0=Pe(y(0)-1), y'(1)=0
#
# dy=y2
# dy2=Pe(dy-Ry^n)
#################################################################

Reactor<-function(x,y,parms)
{
  list(c(y[2],
         Pe*(y[2]+R*(y[1]^n))))
}

Pe <- 1
R  <- 2
n  <- 2

yini <- function (x,parms) c(x,Pe*(x-1))
x    <- seq(0,1,by=0.01)
sol<-bvpshoot(func=Reactor, yend=c(NA,0), y=yini, x=x, extra=1)
plot(sol,main="Reactor")

#################################################################
# Example 3 - final condition is a residual function
# The example for MUSN in Ascher et al., 1995.
# Numerical Solution of Boundary Value Problems for Ordinary Differential
# Equations, SIAM, Philadelphia, PA, 1995.
# MUSN is a multiple shooting code for nonlinear BVPs.
#
# The problem is
#      u' =  0.5*u*(w - u)/v
#      v' = -0.5*(w - u)
#      w' = (0.9 - 1000*(w - y) - 0.5*w*(w - u))/z
#      z' =  0.5*(w - u)
#      y' = -100*(y - w)
#   on the interval [0 1] and with boundary conditions:
#      u(0) = v(0) = w(0) = 1,  z(0) = -10,  w(1) = y(1)
#################################################################

musn <- function(t,Y,pars)
{
  with (as.list(Y),
  {
   du=0.5*u*(w-u)/v
   dv=-0.5*(w-u)
   dw=(0.9-1000*(w-y)-0.5*w*(w-u))/z
   dz=0.5*(w-u)
   dy=-100*(y-w)
   return(list(c(du,dv,dw,dy,dz)))
  })
}

#--------------------------------
# Residuals
#--------------------------------
res  <- function (Y,yini,pars)  with (as.list(Y), w-y)

#--------------------------------
# Initial values; y= NA= not available
#--------------------------------

init <- c(u=1,v=1,w=1,y=NA,z=-10)
sol  <-bvpshoot(y= init, x=seq(0,1,by=0.05), func=musn,
           yend=res, guess=1, atol=1e-10, rtol=0)
pairs(sol, main="MUSN")

#################################################################
# Example 4 - solve also for unknown parameter
# Find the 4th eigenvalue of Mathieu's equation:
# y''+(lam-10cos2t)y=0   on the interval [0,pi]
# y(0)=1, y'(0)=0  and y'(pi)=0
# The 2nd order problem is rewritten as 2 first-order problems:
# dy=y2
# dy2= -(lam-10cos(2t))*y
#################################################################

mathieu<- function(t,y,lam)
{
 list(c(y[2],-(lam-10*cos(2*t))*y[1]))
}

yini <- c(1,0)   # initial condition(y1=1,dy=y2=0)
yend <- c(NA,0)  # final condition at pi (y1=NA, dy=0)

# there is one extra parameter to be fitted: "lam"; its initial guess = 15
Sol  <- bvpshoot(yini=yini, yend=yend, x=seq(0,pi,by=0.01),
        func=mathieu, guess=NULL, extra=15)
plot(Sol)
attr(Sol,"roots")  # root gives the value of "lam" (17.10684)
}


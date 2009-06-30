\name{bvptwp}
\alias{bvptwp}
\title{
  Solver for two-point boundary value problems of ordinary differential
  equations, a mono-implicit Runge-Kutta (MIRK) formula
}

\description{
  Solves a boundary value problem of a system of ordinary differential
  equations. This is an implementation of the fortran
 code twpbvp, based on mono-implicit Runge-Kutta formulae (MIRK)
  of orders 4, 6 and 8 in a deferred correction framework

  written by J.R. Cash and M.H. Wright.
}

\usage{
bvptwp(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL,
  xguess=NULL, yguess=NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
  leftbc=NULL, islin=FALSE, nmax=1000, colp= NULL, atol=1e-8, ...)
}
\arguments{
  \item{yini }{either a vector with the initial (state) values for the ODE
    system, or \code{NULL}. If \code{yini} is a vector, use \code{NA} for an
    initial value which is not available.
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
  \item{xguess }{Initial grid x, a vector.
  }
  \item{yguess }{First guess values y, corresponding to initial grid
    \code{xguess}; a matrix with number of rows=number of equations, and
    whose number of columns = length of \code{xguess}.
  }
  \item{jacfunc }{jacobian (optional) - an \R-function that evaluates the
    jacobian of \code{func} at point \code{x},
    \code{jacfunc} must be defined as: \code{jac = func(x, y, parms,...)}.
    It should return the partial derivatives of func with respect to y,
    i.e. df(i,j) = dfi/dyj. See last example.
    
    If \code{jacfunc} is \code{NULL}, then a numerical approximation using
    differences is used.
  }
  \item{bound }{boundary function (optional) - only if \code{yini} and
    \code{yend} are not available. An \R function that evaluates
    the i-th boundary element at point \code{x}. It should be defined as:
    \code{bound = func(i, y, parms, ...)}. It should return the i-th
    boundary condition. See last example.
  }
  \item{jacbound }{jacobian of the boundary function (optional) - only if
    \code{bound} is defined. An \R function that evaluates
    the gradient of the i-th boundary element with respect to the state
    variables, at point \code{x}.
    It should be defined as: \code{jacbound = func(i, x, y, parms, ...)}.
    It should return the gradient of the i-th boundary condition.
    See last example.

    If \code{jacbound} is \code{NULL}, then a numerical approximation using
    differences is used.
  }
  \item{leftbc }{only if \code{yini} and \code{yend} are not available: the
    number of left boundary conditions.
  }
  \item{islin }{set to \code{TRUE} if the problem is linear - this will
    speed up the simulation.
  }
  \item{nmax }{maximal number of subintervals.
  }
  \item{colp }{number of points per subinterval.
  }
  \item{atol }{error tolerance, a scalar.
  }
  \item{... }{additional arguments passed to the model functions.
  }
}

\value{
  A matrix with up to as many rows as elements in times and as many columns
  as elements in \code{yini} plus the number of "global" values returned
  in the second element of the return from \code{func}, plus an additional
  column (the first) for the x-value.

  There will be one row for each element in \code{x} unless the solver
  returns with an unrecoverable error.

  If \code{yini} has a names attribute, it will be used to label the columns
  of the output value.

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\details{
  This is an implementation of the method twpbvp, to solve
  two-point boundary value problems of ordinary differential equations.

  A boundary value problem does not have all initial values of
  the state variable specified. Rather some conditions are specified at
  the end of the integration interval.

  The ODEs and boundary conditions are made available through the
  user-provided routines, \code{func} and vectors \code{yini} and \code{yend}
  or (optionally) \code{bound}.

  The corresponding partial derivatives are optionally available through the
  user-provided routines, \code{jacfunc} and \code{jacbound}. Default is that
  they are automatically generated by \R, using numerical differences.

  The user-requested tolerance is provided through \code{tol}. The
  default is $1e^-6$.


  If the function terminates because the maximum
  number of subintervals was exceeded, then it is recommended that'
  the program be run again with a larger value for this maximum.'

  For this method to work, there should be ONLY one solution.
  This means that the number of specified boundary value problems must
  equal the number of unspecified initial conditions (so that there is
  exactly ONE root).
}
\seealso{
  \code{\link{bvpshoot}} for the shooting method
}

\keyword{math}

\references{
  J.R. Cash and M.H. Wright, A deferred correction method for nonlinear
  two-point boundary value problems: implementation and numerical evaluation,
  SIAM J. Sci. Stat. Comput., 12 (1991) 971 989.

  \url{http://www.ma.ic.ac.uk/~jcash/BVP\_software}
}

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
sol  <- as.data.frame(bvptwp(yini=init,x=seq(-0.1,0.1,by=0.001),
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
sol  <- bvptwp(yini=c(1,NA),yend=c(-0.5,NA),x=x,func=f2,guess=1)
plot(sol)

# add the analytic solution
curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x),add=TRUE,type="l")


#################################################################
# Problem 2  - solved with specification of boundary, and jacobians
# d4y/dx4 =R(dy/dx*d2y/dx2 -y*dy3/dx3)
# y(0)=y'(0)=0
# y(1)=1, y'(1)=0
#
# dy/dx  = y2
# dy2/dx = y3    (=d2y/dx2)
# dy3/dx = y4    (=d3y/dx3)
# dy4/dx = R*(y2*y3 -y*y4)
#################################################################

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
sol  <- as.data.frame(bvptwp(x=seq(0,1,by=0.01), leftbc=2,
        func=f2, guess=1,R=R,
        bound=g2, jacfunc=df2, jacbound=dg2))
plot(sol[,1:2])

}


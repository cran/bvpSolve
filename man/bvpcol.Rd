\name{bvpcol}
\alias{bvpcol}
\title{
  Solves multipoint boundary value problems of ordinary differential equations,
  using a collocation method.
}
\description{
  Solves Boundary Value Problems For Ordinary
  Differential Equations By Collocation.

  It is possible to solve stiff systems, by using an automatic continuation strategy

  This is an implementation of the fortran codes colsys.f and colnew.f written by
  respectively U. Ascher, J. christiansen  and R.D. Russell (colsys), 
  and U. Ascher and G. Bader (colnew).

  The continuation strategy is an implementation of the fortran code colmod
  written by J.R. Cash, M.H. Wright and F. Mazzia.
}

\note{
  colnew.f (Bader and Ascher, 1987), is a modification of the code
  colsys.f (Ascher, Christiansen and Russell, 1981), which incorporates 
  a new basis representation replacing b-splines, and improvements for
  the linear and nonlinear algebraic equation solvers. To toggle on/off
  colsys, set \code{bspline} = \code{TRUE}/\code{FALSE}

  colmod is a revised version of the package colnew by
  Bader and Ascher (1987), which in turn is a modification of the package
  colsys by Ascher, Christiansen and Russell (1981). Colmod has been
  adapted to allow an automatic continuation strategy to
  be used (Cash et al., 1995).

  The mesh selection algorithm used in colmod differs from that used in colnew
}

\usage{
bvpcol (yini = NULL, x, func, yend = NULL, parms = NULL, 
        order = NULL, ynames = NULL, xguess = NULL, yguess = NULL, 
        jacfunc = NULL, bound = NULL, jacbound = NULL, 
        leftbc = NULL, posbound = NULL, islin = FALSE, nmax = 1000, 
        ncomp = NULL, atol = 1e-8, colp = NULL, bspline = FALSE,
        fullOut = TRUE, dllname = NULL, initfunc = dllname, 
        rpar = NULL, ipar = NULL, nout = 0, forcings = NULL, 
        initforc = NULL, fcontrol = NULL, verbose = FALSE,
        epsini = NULL, eps = epsini,...)
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
    vector containing the values of the equations. In case where the equations
    are first-order, this will be the derivatives of \code{y} with respect to
    \code{x}. After this can come global values that are required at
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
 \item{xguess }{Initial grid \code{x}, a vector. 
    If \code{xguess} is given, so should \code{yguess} be.
    
    Supplying \code{xguess} and \code{yguess}, based on results from a previous
    (simpler) BVP-ODE can be used for model continuation, see example 2 of
    \code{bvptwp}.
  }
  \item{yguess }{First guess values of \code{y}, corresponding to initial grid
    \code{xguess}; a matrix with number of rows equal to the number of
    variables, and whose number of columns equals the length of \code{xguess}.
    
    if the rows of \code{yguess} have a names attribute, the names will be 
    available within the functions and used to label the output matrix.
    
    It is also allowed to pass the output of a previous run for continuation.
    This will use the information that is stored in the attributes 
    \code{istate} and \code{rstate}. It will only work when for the previous 
    run, \code{fullOut} was set equal to \code{TRUE} (the default).
    In this case, \code{xguess} need not be provided.
    
    See example 3b.
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
    See examples.

    If \code{jacbound} is \code{NULL}, then a numerical approximation using
    differences is used. This is the default.
  }
  \item{leftbc }{only if \code{yini} and \code{yend} are not available and
    \code{posbound} is not specified: the number of left boundary conditions.
  }
  \item{posbound }{only used if \code{bound} is given: a vector with the position 
    (in the mesh) of the boundary conditions - its values should be sorted - and 
    it should be within the range of \code{x};
    (posbound corresponds to fortran input "Zeta" in the colnew/colsys 
     FORTRAN codes.  )
    See last example.
    Note that two-point boundary value problems can also be specified via 
    \code{leftbc} (which is simpler). 
  }
  \item{islin }{set to \code{TRUE} if the problem is linear - this will
    speed up the simulation.
  }
  \item{nmax }{maximal number of subintervals during the calculation.
  }
  \item{order }{the order of each derivative in \code{func}. The default is that
    all derivatives are 1-st order, in which case \code{order} can be set = 
    \code{NULL}.
    
    For higher-order derivatives, specifying the order can improve computational
    efficiency, but this interface is more complex. 
    
    If \code{order} is not \code{NULL}, the number of equations in \code{func} 
    must equal the length of \code{order}; the summed values of \code{order}
    must equal the number of variables (ncomp). 
    The jacobian as specified in \code{jacfunc}
    must have number of rows = number of equations and number of columns = number 
    of variables.                            
    \code{bound} and \code{jacbound} remain defined in the number of variables.
    See example 3 and 3b. 
  }   
  \item{ncomp }{used if the model is specified by compiled code, the number of
    components (or equations). See package vignette \code{"bvpSolve"}.
    
    Also to be used if the boundary conditions are specified by \code{bound}, and 
    there is no \code{yguess}
  }
  \item{atol }{ error tolerance, a scalar.
  }
  \item{colp }{number of collocation points per subinterval.
  }
  \item{bspline }{if \code{FALSE}, then code colnew is used the default, 
    if \code{TRUE}, then fortran code colsys is used. Code colnew 
    incorporates a new basis representation, while colsys uses b-splines.
  }
  \item{fullOut }{if set to \code{TRUE}, then the collocation output required
    e.g. for continuation will be returned in attributes \code{rwork} and 
    \code{iwork}. Use \code{attributes(out)\$rwork, attributes(out)\$rwork} 
    to see their contents
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions referred to in \code{func}, \code{jacfunc},
    \code{bound} and  \code{jacbound}.  Note that ALL these subroutines
    must be defined in the shared library; it is not allowed to merge R-functions
    with compiled functions.
    
    See package vignette \code{"bvpSolve"} or deSolve's package vignette
    \code{"compiledCode"}.
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
  \item{verbose }{when \code{TRUE}, writes messages (warnings) to the screen.
  }
  \item{... }{additional arguments passed to the model functions.
  }
}
\value{
  A matrix of class \code{bvpSolve}, with up to as many rows as elements in 
  \code{x} and as many columns
  as elements in \code{yini} plus the number of "global" values returned
  in the second element of the return from \code{func}, plus an additional
  column (the first) for the \code{x}-value.

  There will be one row for each element in \code{x} unless the solver returns
  with an unrecoverable error.

  If \code{ynames} is given, or \code{yini}, \code{yend} has a names attribute, 
  or \code{yguess} has named rows, the names will be used to label the
  columns of the output value.

  The output will also have \code{attributes} \code{istate} and \code{rstate}
  which contain the collocation output required e.g. for continuation of a 
  problem, unless \code{fullOutput} is \code{FALSE}
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\details{
  If \code{eps} does not have a value, then the method is based on an implementation
  of the Collocation methods called "colnew" and
  "colsys" to solve multi-point boundary value problems of ordinary
  differential equations.
  
  The ODEs and boundary conditions are made available through the
  user-provided routines, \code{func} and vectors \code{yini} and \code{yend}
  or (optionally) \code{bound}. \code{bvpcol} can also solve multipoint 
  boundary value problems (see last example).

  The corresponding partial derivatives are optionally available through the
  user-provided routines, \code{jacfunc} and \code{jacbound}. Default is that
  they are automatically generated by \R, using numerical differences.

  The user-requested tolerance is provided through \code{atol}.

  If the function terminates because the maximum
  number of subintervals was exceeded, then it is recommended that
  'the program be run again with a larger value for this maximum.'

  If \code{eps} does have a value, then the method is based on an implementation
  of the Collocation methods called "colmod".
  The type of problems which this is designed to solve typically
  involve a small positive parameter 0 < eps << 1.
  As eps becomes progressively smaller, the problem normally becomes
  increasingly difficult to approximate numerically (for example, due
  to the appearance of narrow transition layers in the profile of the
  analytic solution).

  The idea of continuation is to solve a chain of problems in which the
  parameter eps decreases monotonically towards some desired value.
  That is, a sequence of problems is attempted to be solved:

  epsini  > eps1  > eps2  > eps3  >  .....  > eps  > 0

  where \code{epsini} is a user provided starting value and \code{eps} is a
  user desired final value for the parameter.
}

\seealso{
\code{\link{bvpshoot}} for the shooting method

\link{bvptwp} for a MIRK formula

\code{\link{diagnostics.bvpSolve}}, for a description of diagnostic messages
  
\code{\link{approx.bvpSolve}}, for approximating solution in new values
  
\code{\link{plot.bvpSolve}}, for a description of plotting the output of the 
  BVP solvers.
}

\references{
  U. Ascher, J. Christiansen and R. D. Russell,  (1981)
  collocation software for boundary-value odes,
  acm trans. math software 7, 209-222.

  G. Bader and U. Ascher, (1987)
  a new basis implementation for a mixed order
  boundary value ode solver,
  siam j. scient. stat. comput. 8, 487-483.

  U. Ascher, J. Christiansen and R.D. Russell, (1979)
  a collocation solver for mixed order
  systems of boundary value problems,
  math. comp. 33, 659-679.

  U. Ascher, J. Christiansen and R.D. Russell, (1979)
  colsys - a collocation code for boundary value problems,
  lecture notes comp.sc. 76, springer verlag,
  B. Childs et. al. (eds.), 164-185.

  J. R. Cash, G. Moore and R. W. Wright, (1995)
  an automatic continuation strategy for the solution of
  singularly perturbed linear two-point boundary value problems,
  j. comp. phys. 122, 266-279.
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
fun <- function(t, y, pars) { 
  dy1 <- y[2]
  dy2 <- - 3 * p * y[1] / (p+t*t)^2
  return(list(c(dy1,
                dy2))) }

# parameter value
p    <- 1e-5

# initial and final condition; second conditions unknown
init <- c(-0.1 / sqrt(p+0.01), NA)
end  <- c( 0.1 / sqrt(p+0.01), NA)

# Solve bvp
sol  <- bvpcol(yini = init, yend = end, 
               x = seq(-0.1, 0.1, by = 0.001), func = fun)
plot(sol, which = 1)

# add analytical solution
curve(x/sqrt(p+x*x), add = TRUE, type = "p")

diagnostics(sol)

zoom <- approx(sol, xout = seq(-0.005, 0.005, by  = 0.0001))
plot(zoom, which = 1, main = "zoom in on [-0.0005,0.0005]")


## =============================================================================
## Example 1b: 
## Same problem, now solved as a second-order equation 
## and with different value of "p".
## =============================================================================

fun2 <- function(t, y, pars)
{ dy <- - 3 * p * y[1] / (p+t*t)^2
  list(dy)
}

p <- 1e-4
sol2  <- bvpcol(yini = init, yend = end, order = 2, 
                x = seq(-0.1, 0.1, by = 0.001), func = fun2)

# plot both runs at once:
plot(sol, sol2, which = 1)

## =============================================================================
## Example 1c: simple
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
 dy2 <- -1/x * y[2]- (1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy, dy2))
}

x    <- seq(1, 6, 0.1)
sol  <- bvpcol(yini = c(1, NA), yend = c(-0.5, NA), bspline = TRUE,
               x = x, func = f2)
plot(sol, which = 1)

# add the analytic solution
curve(0.0588713*cos(x)/sqrt(x) + 1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x), add = TRUE, type = "l")


## =============================================================================
## Example 2. Uses continuation
## Test problem 24
## =============================================================================

Prob24<- function(t, y, ks) {     #eps is called ks here
  A <- 1+t*t
  AA <- 2*t
  ga <- 1.4
  list(c(y[2],(((1+ga)/2 -ks*AA)*y[1]*y[2]-y[2]/y[1]-
               (AA/A)*(1-(ga-1)*y[1]^2/2))/(ks*A*y[1])))
}

ini <- c(0.9129, NA)
end <- c(0.375, NA)
xguess <- c(0, 1)
yguess <- matrix(nrow = 2, ncol = 2, 0.9 )

# bvpcol works with eps NOT too small, and good initial condition ...
sol <- bvpcol(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
              xguess = xguess, yguess = yguess,
              parms = 0.1, func = Prob24, verbose = FALSE)

# when continuation is used: does not need a good initial condition
sol2 <- bvpcol(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                  parms = 0.05, func = Prob24,
                  eps = 0.05)
                  
#zoom <- approx(sol2, xout = seq(0.01, 0.02, by  = 0.0001))
#plot(zoom, which = 1, main = "zoom in on [0.01, 0.02]")

sol3 <- bvpcol(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                  parms = 0.01, func = Prob24 , eps = 0.01)

sol4 <- bvpcol(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                  parms = 0.001, func = Prob24, eps = 0.001)

print(system.time(
sol5 <- bvpcol(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
                  parms = 1e-4, func = Prob24, eps = 1e-4)
))

plot(sol, sol2, sol3, sol4, sol5, which = 1, main = "test problem 24",
     lwd = 2)

legend("topright", col = 1:5, lty = 1:5, lwd = 2,
       legend = c("0.05", "0.01", "0.001", "0.0001"), title = "eps")

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

# derivative function: 4 first-order derivatives
f1st<- function(x, y, S) {
  list(c(y[2],
         y[3],
         y[4],
         1/S*(y[2]*y[3] - y[1]*y[4]) ))
}

# jacobian of derivative function
df1st <- function(x, y, S) {
 matrix(nrow = 4, ncol = 4, byrow = TRUE, data = c(
             0,         1,      0,       0,
             0,         0,      1,       0,
             0,         0,      0,       1,
             -1*y[4]/S, y[3]/S, y[2]/S, -y[1]/S))
}

# boundary
g2 <- function(i, y, S)  {
  if (i == 1) return (y[1])
  if (i == 2) return (y[2])
  if (i == 3) return (y[1] - 1)
  if (i == 4) return (y[2])
}

# jacobian of boundary
dg2 <- function(i, y, S)  {
  if (i == 1) return(c(1, 0, 0, 0))
  if (i == 2) return(c(0, 1, 0, 0))
  if (i == 3) return(c(1, 0, 0, 0))
  if (i == 4) return(c(0, 1, 0, 0))
}

# we use posbound to specify the position of boundary conditions
# we can also use leftbc = 2 rather than posbound = c(0,0,1,1)
S    <- 1/100
sol  <- bvpcol(x = seq(0, 1, by = 0.01),
          ynames = c("y", "dy", "d2y", "d3y"),
          posbound = c(0, 0, 1, 1), func = f1st, parms = S, eps = S,
          bound = g2, jacfunc = df1st, jacbound = dg2)

plot(sol)

## =============================================================================
## Example 3b - solved with specification of boundary, and jacobians
## and as a higher-order derivative
## d4y/dx4 =R(dy/dx*d2y/dx2 -y*dy3/dx3)
## y(0)=y'(0)=0
## y(1)=1, y'(1)=0
## =============================================================================

# derivative function: one fourth-order derivative
f4th <- function(x, y, S) {
  list(1/S * (y[2]*y[3] - y[1]*y[4]))
}

# jacobian of derivative function
df4th <- function(x, y, S)  {
  matrix(nrow = 1, ncol = 4, byrow = TRUE, data = c(
             -1*y[4]/S, y[3]/S, y[2]/S, -y[1]/S))
}

# boundary function - same as previous example

# jacobian of boundary - same as previous

# order = 4 specifies the equation to be 4th order
# solve with bspline false
S    <- 1/100
sol  <- bvpcol (x = seq(0, 1, by = 0.01),
          ynames = c("y", "dy", "d2y", "d3y"),
          posbound = c(0, 0, 1, 1), func = f4th, order = 4,
          parms = S, eps = S, bound = g2, jacfunc = df4th,
          jacbound = dg2 )

plot(sol)

# Use (manual) continuation to find solution of a more difficult example
# Previous solution collocation from sol passed ("guess = sol")

sol2  <- bvpcol(x = seq(0, 1, by = 0.01),
          ynames = c("y", "dy", "d2y", "d3y"),
          posbound = c(0, 0, 1, 1), func = f4th,
          parms = 1e-6, order = 4, eps = 1e-6,
          bound = g2, jacfunc = df4th, jacbound = dg2 )

# plot both at same time
plot(sol, sol2, lwd = 2)

legend("bottomright", leg = c(100, 10000), title = "R = ",
         col = 1:2, lty = 1:2, lwd = 2)


## =============================================================================
## Example 4  - a multipoint bvp
## dy1 = (y2 - 1)/2
## dy2 = (y1*y2 - x)/mu
## over interval [0,1]
## y1(1) = 0; y2(0.5) = 1
## =============================================================================

multip <- function (x, y, p) {
  list(c((y[2] - 1)/2, 
         (y[1]*y[2] - x)/mu))
}

bound <- function (i, y, p) {
  if (i == 1) y[2] -1    # at x=0.5: y2=1
  else y[1]              # at x=  1: y1=0
}

mu  <- 0.1
sol <- bvpcol(func = multip, bound = bound, 
              x = seq(0, 1, 0.01), posbound = c(0.5, 1))

plot(sol)

# check boundary value
sol[sol[,1] == 0.5,]
}
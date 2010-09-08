\name{diagnostics.bvpSolve}
\alias{diagnostics.bvpSolve}
\alias{approx}
\alias{approx.bvpSolve}
\title{Prints Diagnostic Characteristics of BVP Solvers}
\description{
  Prints several diagnostics of the simulation to the screen, e.g.
  conditioning parameters
}
\usage{
\S3method{diagnostics}{bvpSolve}(obj, ...)
\S3method{approx}{bvpSolve}(x, xout=NULL, ...)
}
\arguments{
  \item{obj }{the output as produced by \code{bvptwp}, \code{bvpcol} 
    or \code{bvpshoot}.    
  }
  \item{x }{the output as produced by \code{bvpcol}    
  }
  \item{xout }{points \code{x} for which new variable values should 
    be generated. 
  }
  \item{...}{optional arguments to the generic function.
  }
}
\value{
  S3 method \code{diagnostics} prints diagnostic features of the simulation. 
  
  It will depend on the solution method.
  
  In case the problem was solved with \code{bvpshoot}, the diagnostics
  of the initial value problem solver will also be written to screen.
  
  The diagnostics of all solvers include the 
  number of function evaluations, the number of jacobian evaluations, and
  the number of steps
  
  Note that the number of function evaluations are *without* the extra 
  calls performed to generate the ordinary output variables (if present). 

  The diagnostics of both \code{bvptwp} and \code{bvpcol} also include the 
  the number of boundary evaluations and the number of boundary jacobian 
  evaluations.

  In case the method used was
  \code{bvptwp}, will also return the \emph{conditioning parameters}. They are:
  ckappa, ckappa1, ckappa2, sigma and gamma.
  
  See \url{http://pitagora.dm.uniba.it/~mazzia/bvp/tom/lucidinorvegiad.pdf}
  
  The kappa's are based on the Inf-norm, gamma is based on the 1-norm,
  If ckappa and gamma are both of moderate size, the problem is well conditioned.
  If both are large, the problem is ill-conditioned.
  If ckappa is large and gamma is small, the problem is ill-conditioned in the 
  maximum and well conditioned in the 1-norm. This is typical for problems that
  involve different time scales ("stiff" problems).

  S3 method \code{approx} calculates an approximate solution vector at points
  inbetween the original \code{x}-values. If beyond the integration interval,
  it will not extrapolate, but just return the values at the edges.
  This works only when the solution was generated with \code{bvpcol}, and 
  usses information in the arrays rwork and iwork, stored as attributes.
  The returned matrix will be of class "bvpSolve" 

}
\details{
  When the integration output is saved as a data.frame, then the required
  attributes are lost and method \code{diagnostics} will not work anymore.
}
\examples{
## =============================================================================
## Diagnostic messages
## =============================================================================
f2 <- function(x, y, parms) {
 dy  <- y[2]
 dy2 <- -1/x*y[2] - (1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy, dy2))
}

x    <- seq(1, 6, 0.1)
yini <- c(y = 1, dy = NA)
yend <- c(-0.5, NA)

sol   <- bvptwp(yini = yini, yend = yend, x = x, func = f2)
sol2  <- bvpcol(yini = yini, yend = yend, x = x, func = f2)
sol3  <- bvpshoot(yini = yini, yend = yend, x = x, func = f2, guess = 0)

plot(sol, which = "y")
diagnostics(sol)
diagnostics(sol2)
diagnostics(sol3)

## =============================================================================
## approx
## =============================================================================

soldetail <- approx(sol2, xout = seq(2,4,0.01))
plot(soldetail)

# beyond the interval
approx(sol2, xout = c(0,1,2))
approx(sol2, xout = c(6,100))

}
\seealso{
  \link{diagnostics.deSolve} for a description of diagnostic messages of the 
  initial value problem solver as used by \code{bvpshoot}

  \code{\link{plot.bvpSolve}}, for a description of plotting the output of the 
  BVP solvers.
}
\keyword{ utilities }
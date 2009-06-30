\name{bvpSolve-package}
\alias{bvpSolve-package}
\alias{bvpSolve}
\docType{package}
\title{
  Solvers for Boundary Value Problems (BVP) of Ordinary Differential Equations
}
\description{
  Functions that solve boundary value problems of a system of stiff
  ordinary differential equations (ODE)

  The functions provide an interface to the FORTRAN code
  twpbvp written by J.R. Cash and M.H. Wright.

  and also implements a shooting method
}
\details{
  \tabular{ll}{
    Package: \tab bvpSolve\cr
    Type:    \tab Package\cr
    Version: \tab 1.0\cr
    Date:    \tab 2009-06-12\cr
    License: \tab GNU Public License 2 or above\cr
  }

  The system of ODE's is written as an \R function, similar as the initial
  value problems that are solved by integration routines from package
  \code{deSolve}.
  
}
\author{
  Karline Soetaert (Maintainer)
}
\references{
  J.R. Cash and M.H. Wright, (1991) A deferred correction method for nonlinear
  two-point boundary value problems: implementation and numerical evaluation,
  SIAM J. Sci. Stat. Comput. 12,  971-989.
}
\seealso{
  \code{\link{bvptwp}}, a MIRK method to solve
    two-point boundary value problems  (Cash and Wright, 1991).
    
  \code{\link{bvpshoot}}, a shooting method, using solvers from packages
    deSolve and rootSolve.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(bvptwp)
example(bvpshoot)

## open the directory with examples
browseURL(paste(system.file(package = "bvpSolve"), "/examples", sep = ""))

## show package vignette with how to use bvpSolve
## + source code of the vignette
vignette("bvpSolve")
edit(vignette("bvpSolve"))

## show directory with source code of the vignette
browseURL(paste(system.file(package = "bvpSolve"), "/doc", sep = ""))
}
}

\keyword{ package }
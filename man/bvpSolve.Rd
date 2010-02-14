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

  The functions provide (1) an interface to the FORTRAN code
  twpbvpC written by J.R. Cash, F. Mazzia and M.H. Wright,
  and (2) also implements a shooting method.
}
\details{
  \tabular{ll}{
    Package: \tab bvpSolve\cr
    Type:    \tab Package\cr
    Version: \tab 1.1\cr
    License: \tab GNU Public License 2 or above\cr
  }

  The system of ODE's can be written as an \R function, or in compiled code
  (FORTRAN, C), similar as the initial value problems that are solved by
  integration routines from package \code{deSolve}.
  
  A large number of examples have been implemented to show the functionalities
  of the package.
  
  \itemize{
    \item All test problems from the website of J.R. Cash
      \url{http://www.ma.ic.ac.uk/~jcash/BVP_software/PROBLEMS.PDF} are
      implemented in package vignette "bvpTests"
    \item Other test problems, in \R code are in the packages \url{doc/example}
      subdirectory.
    \item Test problems implemented in compiled code can be found in the
      packages \url{doc/dynload} subdirectory.
    \item Still more examples, both in R and compiled code are in the package
      vignetter "bvpSolve".
  }
}
\author{
  Karline Soetaert (Maintainer)
  
  Jeff Cash
  
  Francesca Mazzia
}
\references{
  J.R. Cash and M.H. Wright, (1991) A deferred correction method for nonlinear
  two-point boundary value problems: implementation and numerical evaluation,
  SIAM J. Sci. Stat. Comput. 12,  971-989.

  Cash, J. R. and F. Mazzia, (2005) A new mesh selection algorithm, based on 
   conditioning, for two-point boundary value codes. J. Comput. Appl. Math. 184  
   no. 2, 362--381.

  \url{http://www.ma.ic.ac.uk/~jcash/BVP\_software}
}
\seealso{
  \code{\link{bvptwp}}, a deferred correction method based on mono-
    implicit Runge-Kutta formulas and adaptive mesh refinement, based on 
    conditioning to solve two-point boundary value problems (Cash and Mazzia, 2005).
    
  \code{\link{bvpshoot}}, a shooting method, using solvers from packages
    deSolve and rootSolve.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(bvptwp)
example(bvpshoot)

## open the directory with R- examples
browseURL(paste(system.file(package = "bvpSolve"), "doc/examples", sep = ""))

## open the directory with examples in compiled code
browseURL(paste(system.file(package = "bvpSolve"), "doc/dynload", sep = ""))

## show package vignette with how to use bvpSolve
## + source code of the vignette
vignette("bvpSolve")
edit(vignette("bvpSolve"))

## package vignette with the test problems from J.Cash
## + source code of the vignette
vignette("bvpTests")
edit(vignette("bvpTests"))


## show directory with source code of the vignettes
browseURL(paste(system.file(package = "bvpSolve"), "/doc", sep = ""))
}
}

\keyword{ package }
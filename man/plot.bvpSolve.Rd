\name{plot.bvpSolve}
\alias{plot.bvpSolve}
\alias{print.bvpSolve}

\title{
 Plot and Print Methods for Output of bvp solvers
}
\description{
Plot the output of boundary value solver routines.
}
\usage{
\method{plot}{bvpSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
\method{print}{bvpSolve}(x, ...)
}

\arguments{
  \item{x }{the output of \code{bvpSolve}, as returned by the
    boundary value solvers, and to be plotted.
  }
  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables.
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.
  }
    \item{\dots}{additional graphics arguments passed to
    \code{\link{plot.default}}}
}

\details{
  \code{print.bvpSolve} prints the matrix without the attributes.

  \code{plot.bvpSolve} plots multiple figures on a page.

  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow=c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.

  Other graphical parameters can be passed as well. Parameters
  \code{xlab} and \code{ylab} are vectorized, so it is possible to
  assign specific axis labels to individual plots.
}

\examples{
## =============================================================================
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
   return(list(c(du, dv, dw, dz, dy)))
  })
}

#--------------------------------
# Boundaries
#--------------------------------
bound <- function(i,y,pars) {
  with (as.list(y), {
    if (i ==1) return (u-1)
    if (i ==2) return (v-1)
    if (i ==3) return (w-1)
    if (i ==4) return (z+10)
    if (i ==5) return (w-y)
 })
}

xguess <- seq(0, 1, len = 5)
yguess <- matrix(ncol = 5, (rep(c(1, 1, 1, -10, 0.91), times = 5)) )
rownames(yguess) <- c("u","v","w","z","y")

sol  <- bvpcol (bound = bound, x = seq(0, 1, by = 0.05), 
          leftbc = 4, func = musn, xguess = xguess, yguess = yguess)

mf <- par("mfrow")
plot(sol)
par(mfrow = mf)
}
\seealso{
  \code{\link{diagnostics.bvpSolve}}, for a description of diagnostic messages.

}

\keyword{ hplot }


\name{plot.bvpSolve}
\alias{plot.bvpSolve}

\title{
 Plot Method for Output of bvpSolve
}
\description{
Plot the output of boundary value solver routines.
}
\usage{
\method{plot}{bvpSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
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
  Based on the S3-method of package \code{bvpSolve}.

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

mf <- par("mfrow")
plot(sol)
par(mfrow=mf)
}

\keyword{ hplot }


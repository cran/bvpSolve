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
\method{plot}{bvpSolve}(x, ..., which = NULL, ask = NULL, 
                        obs = NULL, obspar= list())
\method{print}{bvpSolve}(x, ...)
}

\arguments{
  \item{x }{the output of \code{bvpSolve}, as returned by the
    boundary value solvers, and to be plotted.
    
    It is allowed to pass several objects of class
    \code{bvpSolve} after \code{x} (unnamed) - see second example.
  }
  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables, except the first column.
  }
  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.
  }
  \item{obs }{a \code{data.frame} or \code{matrix} with "observed data" that 
    will be added as \code{points} to the plots. \code{obs} can also be a 
    \code{list} with multiple data.frames and/or matrices containing 
    observed data.
    
    By default the first column of an observed data set should contain 
    the \code{time}-variable.  The other columns contain the observed values
    and they should have names that are known in \code{x}.
    
    If the first column of \code{obs} consists of factors or characters (strings),
    then it is assumed that the data are presented in long (database) format,
    where the first three columns contain (name, time, value).  
    
    If \code{obs} is not \code{NULL} and \code{which} is \code{NULL}, 
    then the variables, common to both \code{obs} and \code{x} will be plotted.
  }
  \item{obspar }{additional graphics arguments passed to \code{points}, for 
    plotting the observed data
  }
  \item{\dots}{additional arguments.
  
    The graphical arguments are passed to \code{\link{plot.default}}.

    The dots may also contain other objects of class 
    \code{bvpSolve}, as returned by the boundary value solvers, and 
    to be plotted on the same graphs as \code{x}  - see second example.
    In this case, \code{x} and and these other objects should be compatible, 
    i.e. the names should be the same and they should have same number of rows.
    
    The arguments after \ldots must be matched exactly.     
  } }

\details{
  \code{print.bvpSolve} prints the matrix without the attributes.

  \code{plot.bvpSolve} plots multiple figures on a page.

  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow = c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.
  Set \code{mfrow} equal to \code{NULL} to avoid the plotting function to 
  change user-defined \code{mfrow} or \code{mfcol} settings.

  Other graphical parameters can be passed as well. Parameters
  are vectorized, either according to the number of plots 
  (\code{xlab, ylab}, \code{main, sub}, \code{xlim, ylim}, \code{log},
   \code{asp, ann, axes, frame.plot},\code{panel.first,panel.last},
   \code{cex.lab,cex.axis,cex.main}) or 
  according to the number of lines within one plot (other parameters 
  e.g. \code{col}, \code{lty}, \code{lwd} etc.) so it is possible to
  assign specific axis labels to individual plots, resp. different plotting 
  style. Plotting parameter \code{ylim}, or \code{xlim} can also be a list 
  to assign different axis limits to individual plots.
  
  Similarly, the graphical parameters for observed data, as passed by 
  \code{obspar} can be vectorized, according to the number of observed 
  data sets.
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

musn <- function(t, Y, pars)  {
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
rownames(yguess) <- c("u", "v", "w", "z", "y")

sol  <- bvpcol (bound = bound, x = seq(0, 1, by = 0.05), 
          leftbc = 4, func = musn, xguess = xguess, yguess = yguess)

mf <- par("mfrow")
plot(sol)
par(mfrow = mf)

## =============================================================================
## Example 2. Example Problem 31 from Jeff Cash's website
## =============================================================================

Prob31 <- function(t, Y, pars)  {
  with (as.list(Y), {
    dy    <- sin(Tet)
    dTet  <- M
    dM    <- -Q/xi
    T     <- 1/cos (Tet) +xi*Q*tan(Tet)
    dQ    <- 1/xi*((y-1)*cos(Tet)-M*T)
    list(c( dy, dTet, dM, dQ))
  })
}

ini <- c(y = 0, Tet = NA, M = 0, Q = NA)
end <- c(y = 0, Tet = NA, M = 0, Q = NA)

# run 1
xi <-0.1
twp  <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob31, atol = 1e-10)

# run 2
xi <- 0.05
twp2 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob31, atol = 1e-10)

# run 3
xi <- 0.01
twp3 <- bvptwp(yini = ini, yend = end, x = seq(0, 1, by = 0.01),
               func = Prob31, atol = 1e-10)

# print all outputs at once
plot(twp, twp2, twp3, xlab = "x", ylab = names(ini))


# change type, colors, ...
plot(twp, twp2, twp3, type = c("l", "b", "p"), 
     main = paste ("State Variable", names(ini)), 
     col = c("red", "blue", "orange"), cex = 2)

## =============================================================================
## Assume we have two 'data sets':
## =============================================================================
# data set in 'wide' format
obs1 <- cbind(time = c(0, 0.5, 1), Tet = c(0.4, 0.0, -0.4))

# data set in 'long' format
obs2 <- data.frame(name = "Tet", time = c(0, 0.5, 1), value = c(0.35, 0.0, -0.35))

plot(twp, twp2, obs = obs1, obspar = list(pch = 16, cex = 1.5))

plot(twp, twp2, obs = list(obs1, obs2), 
    obspar = list(pch = 16, cex = 1.5))

plot(twp, twp2, obs = list(obs1, obs2), which = c("Tet", "Q"),
    obspar = list(pch = 16:17, cex = 1.5, col = c("red", "black"))
    )

}
\seealso{
  \code{\link{diagnostics.bvpSolve}}, for a description of diagnostic messages.

}

\keyword{ hplot }


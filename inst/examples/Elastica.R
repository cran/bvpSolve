#===============================================================================
# Elastica equation
#===============================================================================

Elastica <- function (t, y, pars) {
  return(list(c(cos(y[3]),
                sin(y[3]),
                y[4],
                y[5]*cos(y[3]),
                0)
         ))
}

require(bvpSolve)

x <- seq(0,0.5,len=21)
Sol <- bvptwp(yini = c(x=0,y=0,p=NA,k=0,F=NA), yend = c(NA,0,-pi/2,NA,NA),
              x = x, fun=Elastica, guess=c(0,0))
head(Sol)
pairs (Sol)

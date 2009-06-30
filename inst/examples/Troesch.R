#===============================================================================
# Troesch's equation
#===============================================================================

Troesch <- function (t, y, pars) {
  return(list(c(y[2],
                mu*sinh(mu*y[1]))
                ))
}

mu <- 5

x <- seq(0,1,len=11)
Sol <- bvptwp(yini = c(0,NA), yend = c(1,NA), x = x, fun=Troesch, guess=0)

plot(Sol)

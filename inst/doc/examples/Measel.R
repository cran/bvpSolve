## =============================================================================
##  PROBLEM measels
## Models the spread of measels in three equations
## U. M. Ascher, R. M. R. Mattheij, and R. D. Russell. Numerical Solution of
## Boundary Value Problems for Ordinary Differential Equations. Prentice{Hall,
## Englewood Cliffs, NJ, USA, 1988.
## =============================================================================

require(bvpSolve)

## =============================================================================
## Simple implementation, solved with bvpshoot
## =============================================================================

measel2<-function(t,y,pars,vv)   {
  bet<- 1575*(1+cos(2*pi*t))
  dy1<- mu-bet*y[1]*y[3]
  dy2<- bet*y[1]*y[3]-y[2]/lam
  dy3<-y[2]/lam-y[3]/vv
  
  list(c(dy1,dy2,dy3))
}
mu  <- 0.02
lam <- 0.0279
v   <- 0.1

guess <- c(0.01,0.01,0.01)

res <- function(y,yini,parms,vv) y-yini

Sol <- bvpshoot(fun=measel2, yini=c(y1=NA,y2=NA,y3=NA), yend=res, 
    x=(0:100)/100,vv=v, guess=rep(0.01,3))

plot(Sol)

## =============================================================================
## Complex implementation
## =============================================================================

measel<-function(t,y,pars,vv)  {
  bet<- 1575*(1+cos(2*pi*t))
  dy1<- mu-bet*y[1]*y[3]
  dy2<- bet*y[1]*y[3]-y[2]/lam
  dy3<-y[2]/lam-y[3]/vv
  dy4 <- 0
  dy5<-0
  dy6<-0
  
  list(c(dy1,dy2,dy3,dy4,dy5,dy6))
}
mu  <- 0.02
lam <- 0.0279
v   <- 0.1

guess <- c(0.01,0.01,0.01,0.01, 0.01, 0.01)

dmeasel<-function(t,y,pars,vv)
{
  df <- matrix (data=0,nrow=6,ncol=6)
  bet<- 1575*(1+cos(2*pi*t))
  df[1,1] <-  -bet*y[3]
  df[1,3] <-  -bet*y[1]

  df[2,1] <-  bet*y[3]
  df[2,2] <-  -1/lam
  df[2,3] <-  bet*y[1]

  df[3,2] <- 1/lam 
  df[3,3] <- -1/vv
  
  return(df)
}

bound <- function(i, y, pars,vv) {
  if ( i == 1 | i == 4) return(y[1]-y[4])
  if ( i == 2 | i == 5) return(y[2]-y[5])
  if ( i == 3 | i == 6) return(y[3]-y[6])  
  }

dbound <- function(i, y, pars,vv) {
  if ( i == 1 | i == 4) return(c(1,0,0,-1,0,0))
  if ( i == 2 | i == 5) return(c(0,1,0,0,-1,0))
  if ( i == 3 | i == 6) return(c(0,0,1,0,0,-1))
  }

print(system.time(
sola <- bvpshoot(fun=measel, bound=bound, 
    x=(0:100)/100,leftbc=3,vv=v, ncomp=6)
))

sol <- bvptwp(fun=measel, bound=bound, xguess=sola[,1], yguess=t(sola[,-1]),
    x=(0:100)/100,leftbc=3,vv=v, ncomp=6, nmax=100000, atol=1e-4)

plot(sol)

print(system.time(
Sol <- bvpshoot(fun=measel, bound=bound, jacbound = dbound, jacfunc=dmeasel, 
    x=(0:100)/100,leftbc=3,vv=v, ncomp=6)
))

print(system.time(
Sol2 <- bvpshoot(fun=measel, bound=bound, jacbound = dbound, jacfunc=dmeasel, 
    x=(0:100)/100,leftbc=3,vv=v, guess=rep(1,6))
))

plot(Sol2)


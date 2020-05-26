## DO NOT USE

## =============================================================================
## Example 1.
## Find the 4th eigenvalue of Mathieu's equation:
## y''+(lam-10cos2t)y=0   on the interval [0,pi]
## y(0)=1, y'(0)=0  and y'(pi)=0
##
## 2nd order problem is rewritten as:
## dy=y2
## dy2= -(y(3)-10cos(2t))*y
## dy3 = 0 
## =============================================================================

## ------------------------
## Problem specifications
## ------------------------

x      <- seq(0, pi, by = 0.01)
init   <- c(y = 1,dy = 0, lambda = NA)
xguess <- c(0, 1, 2*pi)
yguess <- matrix(nr = 3, rep(1, 9))
rownames(yguess) <- c("y", "dy", "lambda")

## ------------------------
## Inline fortran derivative
## ------------------------

Mathieu <- "
       F(1) = y(2)
       F(2) = -(y(3)-10*cos(2*x))*y(1)
       F(3) = 0
"
mathieu2 <- compile.bvp(Mathieu)   # compile it

## ------------------------
## Solution
## ------------------------

print(system.time(                  # run it
sol <- bvptwp(yini = init, yend = c(NA, 0, NA), x = x,
        func = mathieu2, xguess = xguess, yguess = yguess)
))
plot(sol)

## =============================================================================
## PROBLEM measels
## Models the spread of measels in three equations
## U. M. Ascher, R. M. R. Mattheij, and R. D. Russell. Numerical Solution of
## Boundary Value Problems for Ordinary Differential Equations. Prentice{Hall,
## Englewood Cliffs, NJ, USA, 1988.
## =============================================================================

## ------------------------
## Problem specifications
## ------------------------

mu  <- 0.02
lam <- 0.0279
v   <- 0.1
x <- seq (0, 1, by = 0.01)
yguess <- matrix(ncol = length(x), nrow = 6, data = 1)
rownames(yguess) <- paste("y", 1:6, sep = "")

## ------------------------
## Inline Fortran
## ------------------------

Measel <- "
   vv = rpar(1)
   mu = rpar(2)
   lam = rpar(3)
   bet = 1575d0*(1.+cos(2*pi*x))
   F(1) = mu - bet*y(1)*y(3)
   F(2) = bet*y(1)*y(3) - y(2)/lam
   F(3) = y(2)/lam-y(3)/vv
   F(4) = 0.d0
   F(5) = 0.d0
   F(6) = 0.d0
"


bound <- "
  if ( i == 1 .OR. i == 4) g = (y(1) - y(4))
  if ( i == 2 .OR. i == 5) g = (y(2) - y(5))
  if ( i == 3 .OR. i == 6) g = (y(3) - y(6))  
"
cMeasel <- compile.bvp(func = Measel,  bound = bound,
  header = "  DOUBLE PRECISION :: vv, mu, lam, bet \n  DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793116D0")

## ------------------------
## Solution
## ------------------------

print(system.time(
  sol <- bvptwp(func = cMeasel, 
    xguess = x, yguess = yguess,
    x=x, leftbc = 3, rpar = c(v, mu, lam), ncomp = 6, 
    nmax = 100000, atol = 1e-8)
))

plot(sol)
## =============================================================================
##   This is the example for MUSN in U. Ascher, R. Mattheij, and R. Russell,
##   Numerical Solution of Boundary Value Problems for Ordinary Differential
##   Equations, SIAM, Philadelphia, PA, 1995.  MUSN is a multiple shooting
##   code for nonlinear BVPs.  The problem is
##
##      u' =  0.5*u*(w - u)/v
##      v' = -0.5*(w - u)
##      w' = (0.9 - 1000*(w - y) - 0.5*w*(w - u))/z
##      z' =  0.5*(w - u)
##      y' = -100*(y - w)
##
##   The interval is [0 1] and the boundary conditions are
##
##      u(0) = v(0) = w(0) = 1,  z(0) = -10,  w(1) = y(1)
##
## note: there are two solutions...
## =============================================================================

## ------------------------
## Problem specifications
## ------------------------

x <- seq(0, 1, by = 0.05)
xguess <- seq(0, 1, len = 5)
yguess <- matrix(nc = 5,data = rep(c(1,1,1,-10,0.91),5))
rownames(yguess) <- c("u", "v", "w", "z", "y")

## ------------------------
## Inline Fortran
## ------------------------

fmusn = "
   F(1) = 0.5*y(1)*(y(3)-y(1))/y(2)
   F(2) = -0.5*(y(3)-y(1))
   F(3) = (0.9-1000*(y(3)-y(5))-0.5*y(3)*(y(3)-y(1)))/y(4)
   F(4) = 0.5*(y(3)-y(1))
   F(5) = -100*(y(5)-y(3))
"


# boundaries 
fbound <- "
    if (i == 1) g = (y(1)-1)
    if (i == 2) g = (y(2)-1)
    if (i == 3) g = (y(3)-1)
    if (i == 4) g = (y(4)+10)
    if (i == 5) g = (y(3)-y(5))
"
cmusn <- compile.bvp(func = fmusn, bound = fbound)

## ------------------------
## Solution
## ------------------------

print(system.time(
  sol <- bvptwp(yini = NULL, x = x, func = cmusn, 
                xguess = xguess, yguess = yguess, leftbc = 4, atol = 1e-10)
))
plot(sol)


## =============================================================================
##   This is a nerve impulse model considered in Example 7.1 of
##   R. Seydel, From Equilibrium to Chaos, Elsevier, 1988.  The
##   differential equations
##
##      y1' = 3*t*(y1 + y2 - 1/3*y1^3 + lambda)
##      y2' = -t*(y1 - 0.7 + 0.8*y2)/3
##
##   are to be solved subject to periodic boundary conditions
##
##      y1(0) = y1(1)
##      y2(0) = y2(1)
##
##
##   The parameter lambda has the value -1.3.  The
##   period T is unknown, i.e the length of the interval is not known.
##
##   An extra equation to estimate T is: -T(y1(0)-0.7+0.8*y2(0))/3=1
##
## =============================================================================

## ------------------------
## Problem specifications
## ------------------------

xguess <- seq(0, 1, by = 0.1)
yguess <- matrix(nr = 5, nc = length(xguess), 5.)
yguess[1,] <- sin(2*pi*xguess)
yguess[2,] <- cos(2*pi*xguess)

## ------------------------
## Inline Fortran
## ------------------------

fnerve = "
 F(1) = 3.d0*y(3)*(y(1) + y(2) - 1.d0/3*(y(1)**3) - 1.3)
 F(2) = (-1.d0/3)*y(3)*(y(1) - 0.7d0 + 0.8*y(2)) 
 F(3) = 0.d0
 F(4) = 0.d0
 F(5) = 0.d0"

# jacobian
dfnerve  = "
  df(1,1) = 3.0d0*y(3) -3.d0*y(3)*y(1)**2.
  df(1,2) = 3.0d0*y(3)
	df(1,3) = 3.0d0*(y(1) + y(2) - 1.d0/3*(y(1)**3) - 1.3d0)
  df(2,1) = (-1.0d0/3)*y(3)
	df(2,2) = 0.80d0* (-1.0d0/3)*y(3)
  df(2,3) = (-1.0d0/3)*(y(1) - 0.7d0 + 0.8d0*y(2)) 
"

# boundary function
boundnerve <- "
 if (i == 1) g = (y(3)*(-1.d0/3)*(y(1) - 0.7d0 + 0.8d0*y(2)) - 1.d0)
 if (i == 2) g = (y(1)-y(4))
 if (i == 3) g = (y(2)-y(5))
 if (i == 4) g = (y(1)-y(4))
 if (i == 5) g = (y(2)-y(5))
"
cfnerve <- compile.bvp(func = fnerve, jacfunc = dfnerve, bound = boundnerve)

## ------------------------
## Solution
## ------------------------

print (system.time(
Sol  <- bvpcol(func = cfnerve, 
        x = seq(0, 1, by = 0.01), 
        ynames = c("y", "dy", "T", "yi", "yj"),
        leftbc = 3, xguess = xguess, yguess = yguess)
))
plot(Sol)


## ==================  tubular reactor with axial dispersion ===================
## y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
## y'(0) = Pe (y(0)-1),
## y'(1)=0
##
## dy=y2
## dy2=Pe(dy-Ry^n)
##
## The initial condition y'(0) is a function of y(0)
## =============================================================================

Pe <- 1
R  <- 2
n  <- 2

## ------------------------
## Inline Fortran
## ------------------------

freac <- "
  Pe = rpar(1)
  R = rpar(2)

  F(1) = y(2)
  F(2) = Pe * (y(2) + R*(Y(1)**2.))
"

# 2nd order implementation
freac2 <- "
  Pe = rpar(1)
  R = rpar(2)

  F(1) = Pe * (y(2) + R*(Y(1)**2.))
"

jreac2 <- "
  Pe = rpar(1)
  R = rpar(2)

  dF(1, 1) = Pe * R*2 *Y(1)
  dF(1, 2) = Pe * y(2)
"

greac <- "
  Pe = rpar(1)
  if (i == 1) g = (y(2)-Pe*(y(1)-1.d0))
  if (i == 2) g = (y(2))
"
cfreac <- compile.bvp(func = freac, jacfunc = jreac2, bound = greac, 
    header = " DOUBLE PRECISION :: Pe, R")

cfreac2 <- compile.bvp(func = freac2, jacfunc = jreac2, bound = greac, 
    header = " DOUBLE PRECISION :: Pe, R")

## ------------------------
## Solution
## ------------------------

sol <- bvptwp(func = cfreac, x = seq(0, 1, by = 0.01),  
              leftbc = 1, ynames = c("y", "dy"), bound = cgreac, rpar = c(Pe, R))

sol2 <- bvpcol(func = cfreac2, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), bound = cgreac, rpar = c(Pe, R))

sol3 <- bvpcol(func = cfreac2, jacfunc = cjreac2, x = seq(0, 1, by = 0.01), order = 2,
              leftbc = 1, ynames = c("y", "dy"), bound = cgreac, rpar = c(Pe, R))

plot(sol, sol2, sol3)


## =============================================================================
## "Swirling Flow III", a test problem in Ascher, Mattheij,
## and Russell, Numerical Solution of Boundary Value Problems
## for Ordinary Differential Equations", Classics in Applied
## Mathematics Series, SIAM, Philadelphia, 1995].
## g'' = (g f' - f g')/eps
## f'''' = (-ff'''-gg')/eps
## g(0)=-1,f(0)=0,f'(0)=0, g(1)=1,f(1)=0,f'(1)=0
##
## 1st order system (y1=g, y3=f)
## y1' = y2
## y2' = (y1*y4 -y3*y2)/eps
## y3'=y4
## y4'=y5
## y5'=y6
## y6'=(-y3y6-y1y2)/eps
## y1(0)=-1,y3(0)=0,y4(0)=0, y1(1)=1,y3(1)=0,y4(1)=0
## =============================================================================

x       <- seq(0, 1, 0.01)
yini <- c(-1, NA, 0, 0, NA, NA)
yend <- c(1 , NA, 0, 0, NA, NA)

## ------------------------
## Inline Fortran
## ------------------------

fswirl <- "
  eps = rpar(1)
  f(1) = Y(2)
  f(2) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(3) = Y(4)
 	f(4) = Y(5)
 	f(5) = Y(6)
  f(6) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl <- compile.func(fswirl, header = "double precision:: eps")

fswirl2 <- "
  eps = rpar(1)
  f(1) = (Y(1)*Y(4) - Y(3)*Y(2))/eps
  f(2) = (-Y(3)*Y(6) - Y(1)*Y(2))/eps
"
cswirl2 <- compile.func(fswirl2, header = "double precision:: eps")

## ------------------------
## Solution
## ------------------------

print(system.time(sol <- bvptwp(x = x, func = cswirl, 
                  yini = yini, yend = yend, eps = 0.01, parms = 0.01)))

# For successively smaller values of eps:
print(system.time(sol2 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol[,1], yguess = t(sol[,-1]),
                  yini = yini, yend = yend, eps=0.001, epsini = 0.01,
                  parms=0.001)))

print(system.time(sol3 <- bvptwp(x = x, func = cswirl, 
                  xguess = sol2[,1], yguess = t(sol2[,-1]),
                  yini = yini, yend = yend, eps = 0.0001, 
                  epsini = 0.001, parms = 0.0001)))

## ================================================================     comments
## fluid injection problem
## Original definition:
## f'''- R[(f')^2 -f*f''] + A = 0
## h'' + R*f*h' + 1 = 0
## O'' + P*f*O' = 0
## A is unknown, P = 0.7*R
##
## rewritten as:
## df=f1                   #f
## df1=f2                  #f1=f'
## df2=R(f1^2-f*f2)-A      #f2=f''
## dh=h1
## dh1= -Rfh1-1
## dO=O1
## dO1 = O2
## dO2 = -P*f*O1
## dA = 0                  # the constant to be estimated
## ================================================================ end comments

times  <- seq(0, 1, by = 0.1)      
yini   <- c(f=0, f1=0, f2=NA, h=0, h1=NA, O=0, O1=NA, A=NA)
yend   <- c(f=1, f1=0, f2=NA, h=0, h1=NA, O=1, O1=NA, A=NA)

## ------------------------
## Inline Fortran
## ------------------------

ffluid <- "
 R    = rpar(1)
 P    = 0.7d0*R
 f(1) = y(2)
 f(2) = y(3)
 f(3) = R*(y(2)**2-y(1)*y(3))-y(8)
 f(4) = y(5)
 f(5) = -R*y(1)*y(5)-1.d0
 f(6) = y(7)
 f(7) = -P*y(1)*y(7)
 f(8) = 0.d0
"
cfluid <- compile.func(ffluid, header = "Double precision :: R, P")

soltwp <- bvptwp(func = cfluid, x = times, parms = NULL, rpar = 200,
                nmax = 195, cond = TRUE, yini = yini, yend = yend)
diagnostics(soltwp)

## ------------------------
## Solution
## ------------------------

soltwp2 <- bvptwp(func = cfluid, x = times, parms = NULL, rpar = 5000,
                nmax = 1100, cond = TRUE, verbose = TRUE,    #
                yini = yini, yend = yend, xguess = soltwp[,1],
                yguess = t(soltwp[,-1]))
diagnostics(soltwp2)

# plot the results
plot(soltwp, soltwp2, type = "l", lwd = 2)

## =============================================================================
## Standard linear problem with boundary layer at the origin
##
##
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
##
## Solved using shooting and bvptwp
## =============================================================================

## ------------------------
## Inline Fortran
## ------------------------
x <- seq(-1, 1, length.out = 100)
fbnd <- "
  F(1) = y(2)
  F(2) = - 3.d0 * rpar(1) * y(1)/(rpar(1) + x*x)**2.
"
fjac <- "
  df(1, 2) = 1.d0
  df(2, 1) = -3.d0*rpar(1)/(rpar(1) +x*x)**2.
"
fbound <- "
  if (i == 1) then
     g = (y(1) + 0.1/sqrt(rpar(1) + 0.01d0))
  else  
     g = (y(1) - 0.1/sqrt(rpar(1) + 0.01d0))
  end if
"
fjacbound <- "
  dg(1) = 1.d0
"

cbnd      <- compile.func(fbnd)
cjac      <- compile.jacfunc(fjac)
cbound    <- compile.bound(fbound)
cjacbound <- compile.jacbound(fjacbound)

## ------------------------
## Solution
## ------------------------

p    <-1e-5

print(system.time(Sol <- bvptwp(x = x, leftbc = 1, func = cbnd, 
        bound = cbound, jacbound = cjacbound, jacfunc = cjac, 
        ncomp = 2, verbose = TRUE, rpar = p)))
plot(Sol, which = 1)

for (pp in 0:9){
  Soln <- bvptwp(x = x, leftbc = 1, func = cbnd, 
        bound = cbound, jacbound = cjacbound, jacfunc = cjac, 
        ncomp = 2, rpar = 10^(-pp))
  lines(Soln[,1], Soln[,2])      
}  
        

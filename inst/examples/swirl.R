#==============================================================================
# "Swirling Flow III", a test problem in Ascher, Mattheij,
# and Russell, Numerical Solution of Boundary Value Problems
# for Ordinary Differential Equations", Classics in Applied
# Mathematics Series, SIAM, Philadelphia, 1995].
# g'' = (g f' - f g')/eps
# f'''' = (-ff'''-gg')/eps
# g(0)=-1,f(0)=0,f'(0)=0, g(1)=1,f(1)=0,f'(1)=0
#
# 1st order system (y1=g, y3=f)
# y1' = y2
# y2' = (y1*y4 -y3*y2)/eps
# y3'=y4
# y4'=y5
# y5'=y6
# y6'=(-y3y6-y1y2)/eps
# y1(0)=-1,y3(0)=0,y4(0)=0, y1(1)=1,y3(1)=0,y4(1)=0
#==============================================================================
require(bvpSolve)

# Declare global problem dependent parameters.
eps     <- 0.01

neqns   <- 6   #number of differential equations
leftbc  <- 3   #number of boundary conditions at left end of the interval
rightbc <-neqns-leftbc   #number of boundary conditions at the right end
  
fsub <- function (t,Y,pars,eps)
{ return(list(c(f1 = Y[2],
                f2 = (Y[1]*Y[4] - Y[3]*Y[2])/eps,
	              f3 = Y[4],
              	f4 = Y[5],
              	f5 = Y[6],
	              f6 = (-Y[3]*Y[6] - Y[1]*Y[2])/eps)))
}

# Solve the model. Simple call
# initial guesses...
nsub   <- 10
A      <- 0.
B      <- 1.
xguess <- seq(A, B, len=nsub)
yguess <- matrix(nr=neqns, nc=nsub, data=0.)
slope  <- 0.375-0.9129
yguess[1,] <- 2*xguess -1
yguess[2,] <- 2

x     <- seq(A,B,len=100)
#shooting does not work

#collocation does...
print(system.time(Sol <- bvptwp(atol=1e-5,x=x,func=fsub,guess= c(2,0,0),
                  yini=c(-1,NA,0,0,NA,NA), yend=c(1,NA,0,0,NA,NA),eps=eps)))

# xguess have equal number of nodes as x!
print(system.time(Sol <- bvptwp(atol=1e-5,
                  xguess=xguess,yguess=yguess,x=x,func=fsub,
                  yini=c(-1,NA,0,0,NA,NA), yend=c(1,NA,0,0,NA,NA),eps=eps)))
pairs(Sol,col="blue",lty=2)


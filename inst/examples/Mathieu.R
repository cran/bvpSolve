
################################################################################
# Find the 4th eigenvalue of Mathieu's equation:
# y''+(lam-10cos2t)y=0   on the interval [0,pi]
# y(0)=1, y'(0)=0  and y'(pi)=0
# 2nd order problem is rewritten as:
# dy=y2
# dy2= -(lam-10cos(2t))*y
################################################################################

require(bvpSolve)

mathieu<- function(t,y,lambda=15)
{
 list(c(y[2],-(lambda-10*cos(2*t))*y[1]))
}

#----------------------
# Solution method 1
#  **  shooting  **
#----------------------
x = seq(0,pi,by=0.01)

init <- c(1,0)
sol  <- bvpshoot(yini=init,yend=c(NA,0),x=x,
        func=mathieu, guess=NULL,  extra=15)
plot(sol)

#----------------------
# Solution method 2
# multiroot + collocation
#----------------------

cost <- function(X)
{  sol<- bvptwp(yini=c(1,NA), yend=c(NA,0),x=c(0,pi),parms=X,
        func=mathieu,guess=0)
  return(sol[2,3])  # y2[0]=0
}

# find the root
lam <- multiroot(f=cost,start=15)

# solve the mode with this root...
Sol<- bvptwp(yini=c(1,0), yend=c(NA,NA),x=x,parms=lam$root,
        func=mathieu,atol=1e-10)
lines(Sol,col="red")


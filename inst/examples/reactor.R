####################### PROBLEM 1 #########################
# tubular reactor with axial dispersion
# y''=Pe(y'+Ry^n) Pe=1,R=2,n=2
# y'0=Pe(y(0)-1), y'(1)=0

# dy=y2
# dy2=Pe(dy-Ry^n)
#
# Can be solved with bvpshoot and bvpcolmod, NOT with bvptwp
####################### PROBLEM 1 #########################

Reactor<-function(x,y,parms)
{
  list(c(y[2],Pe*(y[2]+R*(y[1]^n))))
}

res <- function(y,yini,pars)
   y[2]  # y'=0

Pe <- 1
R  <- 2
n  <- 2

yini <- function (x,parms) c(x,Pe*(x-1))

sol<-bvpshoot(func=Reactor,yend=c(NA,0),yini=yini,x=seq(0,1,by=0.01),extra=1)
plot(sol)



############################# TEST PROBLEM 3 ###################################
# Standard linear problem with boundary layer at the origin
#
# Example 1: simple standard problem
# solve the BVP ODE:
# d2y/dt^2=-3py/(p+t^2)^2
# y(t= -0.1)=-0.1/sqrt(p+0.01)
# y(t=  0.1)= 0.1/sqrt(p+0.01)
# where p = 1e-5
#
# analytical solution y(t) = t/sqrt(p + t^2).
#
# The problem is rewritten as a system of 2 ODEs:
# dy=y2
# dy2=-3p*y/(p+t^2)^2
#
# Solved using shooting and collocation method
################################################################################
 require(bvpSolve)

#--------------------------------
# Derivative function
#--------------------------------
fun <- function(t,y,pars)
{ dy1 <-y[2]
  dy2 <- - 3*p*y[1]/(p+t*t)^2
  return(list(c(dy1,
             dy2))) }

# parameter value
p    <-1e-5

# initial and final condition; second conditions unknown
init <- c(-0.1/sqrt(p+0.01), NA)
end  <- c(0.1/sqrt(p+0.01), NA)
x    <- seq(-0.1,0.1,by=0.001)

#---------------------
# Solution method 1
#  **  shooting  **
#---------------------

print(system.time(sol  <- as.data.frame(bvpshoot(yini=init,x=x,
       func=fun, yend=end, guess=1, atol=1e-10))))
plot(sol$time,sol[,2],type="l")

# add analytical solution
curve(x/sqrt(p+x*x),add=TRUE,type="p")

#---------------------
# Solution method 2
# MIRK method -simple input
#---------------------

print(system.time(Sol2<- as.data.frame(bvptwp(yini=init,x=x,
      func=fun, yend=end, guess=1,atol=1e-10))))
lines(Sol2[,1],Sol2[,2],type="l",col="red")

#-------------------------------
# Solution method 3
#  ** MIRK -full input **
#--------------------------------

# Now using bvptwp with complex input functions
# the jacobian
jacfun <- function(t,y,pars)
{ return(matrix(nr=2,nc=2,byrow=TRUE,data=c(0            ,1,
                                           -3*p/(p+t*t)^2,0)) +1e-10*(runif(4)-0.5))
}

# the boundaries
boundfun <- function (i,y,pars)
  if (i ==1) return(y[1] +0.1/sqrt(p+0.01)) else return(y[1]-0.1/sqrt(p+0.01))

# the jacobian of the boundaries
boundjac <- function (i,y,pars) return(c(1,0))
#

print(system.time(Sol <- as.data.frame(bvptwp(x=x,leftbc=1,func=fun, guess=1,
        bound=boundfun,jacbound=boundjac,jacfunc=jacfun), verbose=TRUE)))
lines(Sol[,1],Sol[,2],type="l",col="blue")



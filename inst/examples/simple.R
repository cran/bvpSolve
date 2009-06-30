
############################# TEST PROBLEM  ###################################
# Standard problem from the twpbvp source code
################################################################################
 require(bvpSolve)

#--------------------------------
# Derivative function
#--------------------------------
fun <- function(t,y,pars)
{
  return(list(c(
     f1 = y[2],
     f2 = (-exp(y[1])*y[2] +
            0.5*pi*sin(0.5*pi*t)*exp(2.*y[1]))/eps)))
}

# parameter value
eps    <-0.01

# initial and final condition; second conditions unknown
x <- seq(from=0,to=1,len=101)

yini <- c(0,NA)
yend <- c(0,NA)

#---------------------
# Solution method 1
#  **  shooting  **
#---------------------

print(system.time(sol  <- as.data.frame(bvpshoot(yini=yini,yend=yend,x=x,
       guess=1,func=fun, atol=1e-6))))
plot(sol$time,sol[,2],type="l")

#---------------------
# Solution method 2
# Collocation methods
#---------------------

print(system.time(sol2  <- as.data.frame(bvptwp(yini=yini,yend=yend,x=x,
       guess=1,func=fun, atol=1e-6))))
points(sol2[,1:2])


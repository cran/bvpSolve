################################################################################
# solve d2y/dx2 + 1/x*dy/dx + (1-1/(4x^2)y = sqrt(x)*cos(x), y(1)=1, y(6)=-0.5
# Analytic solution = y(x)=0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+1/4*x^(3/2)*sin(x)

# Write as set of 2 odes
# dy/dx = y2
# dy2/dx  = - 1/x*dy/dx - (1-1/(4x^2)y + sqrt(x)*cos(x)
################################################################################

f2 <- function(x,y,parms)
{
 dy  <- y[2]
 dy2 <- -1/x*y[2]-(1-1/(4*x^2))*y[1] + sqrt(x)*cos(x)
 list(c(dy,dy2))
}

x    <- seq(1,6,0.1)
print(system.time(sol  <- bvpshoot(yini=c(1,NA),yend=c(-0.5,NA),x=x,func=f2 ,guess=1)))
plot(sol)

# the analytic solution
curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+0.740071*sin(x)/sqrt(x)+
      1/4*x^(3/2)*sin(x),add=TRUE,type="l")

print(system.time(sol2  <- bvptwp(yini=c(1,NA),yend=c(-0.5,NA),x=x,func=f2 ,guess=1)))
points(sol2,col="red",pch="+")


## =============================================================================
##
## flat moon problem
## finding the optimal ascent trajectory for launch from a flat moon to 
## a 100 nautical mile circular orbat
##
## assumptions: constant thrust, constant mass, no drag
##  uniform flat-moon gravity
##
## implementation: Karline Soetaert
## =============================================================================

require(bvpSolve)

x    <- seq(0, 1, length.out = 100)   # non-dimensional time vector

h    <- 185.2e3 # metres, final altitude at 100nmi
VC   <- 1627.0  # m/s circular spped at 100nmi
g    <- 1.62    # m/s2 gravitational acceleration
t2w  <- 3       # thrust to weight ratio for ascent vehicle, in lunar g's
A    <- g* 3    # accelleration

# boundary conditions
yini <- c(x  = 0, y = 0, Vx = 0,  Vy = 0, l2 = NA, l4 = NA, tau =NA)
yend <- c(x = NA, y = h, Vx = VC, Vy = 0, l2 = NA, l4 = NA, tau = NA)


flatmoon <- function(t, y, p) {
   res <- c(y[3],                                        #dx /dtau
            y[4],                                        #dy /dtau
            A * (1/sqrt(1+y[6]^2)),                      #dVx/dtau
            A * (y[6]/sqrt(1+y[6]^2)) - g,               #dVy/dtau
            0.,                                          #dl2/dtau
            -y[5],                                       #dl4/dtau
            0.                                           #dt /dtau
           ) 
    list (res * y[7])
}

print (system.time(
sol <-  bvptwp(func = flatmoon, x = x, yini = yini, yend = yend, xguess = x,
  yguess = matrix (nrow = 7, ncol = length(x), data = c(0, 0, 0, 0, 0, 0, 700)))
))

sol.unscaled <- sol 
sol.unscaled[,1] <- sol[,1]*sol[,ncol(sol)]

plot(sol.unscaled)

fflat <- "
      g = rpar(1)
      A  = 3.d0 * g
      tf = y(7)

      F(1) = Y(3) * tf                             !x
      F(2) = Y(4) * tf                             !y
      F(3) = A*(1.d0/sqrt(1+y(6)**2))* tf          !vx
      F(4) = (A*(y(6)/sqrt(1+y(6)**2)) - g)* tf    !vy
      F(5) = 0.d0                                  !lambda2_bar_dot
      F(6) = -y(5) * tf                            !lambdar4bardot
      F(7) = 0.d0                         !tf
"
cflat <- compile.func(fflat, header = " double precision :: g, A, tf")

print (system.time(
sol2 <- bvptwp(func = cflat, x = x, yini = yini, yend = yend, xguess = x, 
  yguess = matrix (nrow = 7, ncol = length(x), data = c(0, 0, 0, 0, 0, 0, 700)),
  rpar = g,  nmax = 1e5)
))

plot(sol, sol2)
  
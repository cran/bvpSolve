
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using the shooting method
##==============================================================================

bvpshoot<- function(yini, x, func, yend, parms=NULL, guess=NULL, extra=NULL,
    atol=1e-8, rtol=1e-8, maxiter=100, positive =FALSE, method="lsoda", ...)  {


  init <- function(X) {  # initialises yini and parms..
    if (is.function(yini))
      y <<-yini(X,parms,...) else y <- yini
    if (lini>0)
      y[inix] <<- X[1:lini]
    if(lex>0)
      parms[1:lex] <<- X[(lini+1):(lini+lex)]
  }
  
  cost <- function(X,...)  {  # objective function to minimise
    times <- c(x[1], x[length(x)])
    init(X)
    out   <- ode(y=y, times=times, fun=func, parms=parms, method=method,
                 atol=atol, rtol=rtol, ...)
    if (is.function(yend) )
      Res   <- yend(out[nrow(out),2:(ly+1)], y, parms,...)
    else {
      Res <-yend - out[nrow(out),2:(ly+1)]
      Res <- Res[! is.na(Res)]
    }
    return(Res)
  }
  if (is.function(yini))
    y <- yini(extra,parms,...)
  else y <- yini

  ly <- length(y)
  
  inix       <- which (is.na(y))
  lini       <- length(inix)
  lex        <- length(extra)
  #Karline:check this...
  if (lini > 0 & is.null(guess))  {
    warning("estimates for unknown initial conditions not given ('guess'); assuming 0's")
    guess <- rep(0,lini)
  }

  if (lini != length(guess))  {
    if (is.null(extra))
      stop("length of guess should be equal to number of NAs in y") else
    if (lex > length(parms))
      stop("length of extra should be smaller than number of parameters")
  }
  if (lini > 0)
    y[inix] <- guess

  if (lini+lex==0)
    stop ("this is not a boundary value problem - use initial value problem solver instead")
    sol <- multiroot(start = c(guess,extra), cost, atol=atol, rtol=rtol,
                   maxiter=maxiter, positive =positive, ...)


  init(sol$root)
  out <- ode (t=x, fun=func, y=y, parms=parms, method=method,
              atol=atol, rtol=rtol, ...)
  attr(out,"istate") <- NULL
  attr(out,"rstate") <- NULL
  attr(out,"roots")  <- data.frame(root=sol$root,
                                   f.root=sol$f.root, iter=sol$iter)
  out
}

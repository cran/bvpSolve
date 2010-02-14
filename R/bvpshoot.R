
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using the shooting method
##==============================================================================

bvpshoot<- function(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL, 
    extra=NULL, jacfunc=NULL, bound=NULL, jacbound=NULL, 
    leftbc=NULL, ncomp = NULL, atol=1e-8, rtol=1e-8, maxiter=100, 
    positive =FALSE, method="lsoda", ...)  {

## ---------------------
## check input
## ---------------------
  if (is.null(yini)   && is.null(bound))
    stop("either 'yini' and 'yend' or 'bound' should be inputted")
  if (!is.null(yini)  && is.null(yend))
    stop("if 'yini' is given, 'yend' should also be given")
  if (!is.null(yini)  && !is.null(bound))
    stop("either 'yini' or bound should be given, not both")
  if (!is.null(bound)  && !is.null(extra))
    stop("cannot combine 'bound' with 'extra'; use 'yini' and 'yend' ")

  lex      <- length(extra)

## ---------------------
## yini or bound
## ---------------------
  if (! is.null(yini))  {    
    # initial value yini; a function or vector
    inity <- function(X,parms) {  # initialises yini and parms..
      if (is.function(yini))
        Y <-yini(X,Parms,...)
      else Y <- yini
      if (lini>0)
        Y[inix] <- X[1:lini]
      names(Y)<-Ynames
      Y
    }
    # parameters; some may be estimated (lex)
    initparms <- function (X) {
      Initparms <- parms
      if(lex>0)
        Initparms[1:lex] <- X[(lini+1):(lini+lex)]
      return(Initparms)  
    }  

    if (is.function(yini))
      y <- yini(extra,parms,...)
    else
      y <- yini
      
    # root function to solve for  
    rootfun <- function(X,...)  {  
      times <- c(x[1], x[length(x)])
      Parms <- initparms(X)
      Y     <- inity(X,Parms)
      out   <- ode(y=Y, times=times, fun=func, jacfunc=jacfunc, 
                 parms=Parms, method=method,
                 atol=atol, rtol=rtol, ...)
    # deviation with yend should be =0             
      if (is.function(yend) )
        Res   <- yend(out[nrow(out),2:(ly+1)], Y, Parms,...)
      else {
        Res <-yend - out[nrow(out),2:(ly+1)]
        Res <- Res[! is.na(Res)]
      }
      return(Res)
    }
    # the jacobian of root function: not specified 
    JacBound <- NULL   

  } else {                  # bound is specified   

    if (is.null(leftbc))
      stop("leftbc should be inputted if bound is given")
    y <- guess  
  
    rootfun <- function(X,...)  {  
      times <- c(x[1], x[length(x)])
      out   <- ode(y=X, times=times, fun=func, jacfunc=jacfunc, 
                   parms=parms, method=method,
                   atol=atol, rtol=rtol, ...)
      Yend <- out[nrow(out),2:(ly+1)]             
      Res <- vector(len=ly)
      for ( i in 1:leftbc) Res[i] <- bound(i,X,parms,...)
      if (leftbc<ly) for (i in  (leftbc+1):ly) Res[i]<- bound(i,Yend,parms,...)
      return(Res)
    }
    JacBound <- NULL   
    if (! is.null(jacbound)) {  
      JAC <- matrix(nr=length(y),nc=length(y),0)  
      JacBound <- function(x,...)  {
        for (i in 1:ly) 
          JAC[i,]<- jacbound(i,x,parms,...)
        return(JAC)  
      }    
    }    
  }
## ---------------------
## names, guess of y
## ---------------------
  
  Ynames <- attr(y,"names")
  if (is.null(Ynames) & ! is.null(yend)) 
    Ynames <- names(yend)
  if (is.null(y)) {
    y <- rep(0,ncomp)
    warning("estimates for initial conditions not given ('guess'); assuming 0's")
  }  
  if (is.null(y)) 
    stop ("either provide 'guess' for initial conditions or 'ncomp', number of compartments")

  ly <-  length(y)
  
  inix     <- which (is.na(y))
  lini     <- length(inix)
  if (lini > 0 & is.null(guess))  {
    warning("estimates for unknown initial conditions not given ('guess'); assuming 0's")
    guess <- rep(0,lini)
  }

  if (! is.null(yini) & lini != length(guess))  {
    if (is.null(extra))
      stop("length of guess should be equal to number of NAs in y") else
    if (lex > length(parms))
      stop("length of extra should be smaller than number of parameters")
  }
  if (lini > 0)
    y[inix] <- guess

  if (is.null(yini) & is.null(guess)) 
    guess <- y
  if (! is.null(yini) & lini+lex==0)
    stop ("this is not a boundary value problem - use initial value problem solver instead")

## ---------------------
## root solver: 
## ---------------------
  # find unknown initial condition + extra parameter  
  sol <- multiroot(start = c(guess,extra), f = rootfun, atol=atol, rtol=rtol,
                   maxiter=maxiter, jacfunc = JacBound, positive =positive, ...)
  
## ---------------------
## Output 
## ---------------------

  if (! is.null(yini) ) {
    Parms <- initparms(sol$root)
    Y     <- inity(sol$root,Parms)
  }  else {
    Parms <- parms 
    Y     <-  sol$root
  }
    
  out <- ode (t=x, fun=func, y=Y, parms=Parms, method=method, jacfunc=jacfunc, 
              atol=atol, rtol=rtol, ...)
              
  attr(out,"istate") <- NULL  # similar attributes of deSolve solvers
  attr(out,"rstate") <- NULL
  attr(out,"roots")  <- data.frame(root=sol$root,
                                   f.root=sol$f.root, iter=sol$iter)
  class(out) <- c("bvpSolve","matrix")  # a boundary value problem
  colnames(out)[1] <- "x"
  attr(out,"name") <- "bvpshoot"

  out
}

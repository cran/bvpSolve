
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using using collocation method "twpbvp"
##==============================================================================

bvptwp<- function(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL,
     xguess=NULL, yguess=NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
     leftbc=NULL, islin=FALSE, nmax=1000, colp=NULL, atol=1e-8, ...)   {

  rho <- environment(func)

  aleft  <- x[1]
  aright <- x[length(x)]
  fixpt  <- x[-c(1,length(x))]
##---------------------
## error checking...
##---------------------
  if (is.null(yini)   && is.null(bound))
    stop("either 'yini' and 'yend' or 'bound' should be inputted")
  if (!is.null(yini)  && is.null(yend))
    stop("if 'yini' is given, 'yend' should also be given")
  if (!is.null(yini)  && !is.null(bound))
    stop("either 'yini' or bound should be given, not both")

##---------------------
## Boundary conditions  specified by yini and yend
##---------------------

  if (! is.null(yini))  {    # yini is either a vector or a function
    y       <- yini
    Y       <- y
    inix    <- which (is.na(y))
    nas     <- length(inix)
    leftbc  <- length(which (!is.na(y)))   # NA when initial condition not known

    if (is.null(guess)&& ! is.null(yguess))
      guess <- yguess[1,inix]

    if (nas > 0 & is.null(guess))  {
      warning("estimates for unknown initial conditions not given ('guess'); assuming 0's")
      guess <- rep(0,nas)
    }

    if (nas != length(guess))
      stop("length of 'guess' should be equal to number of NAs in y")
    if (nas > 0)
      y[inix] <- guess
    inix   <- which (!is.na(yini))
    finalx <- which (!is.na(yend))
  } else   {               # boundaries specified by boundary function 'bound'
    if (is.null(leftbc))
      stop("leftbc should be inputted if bound is given")
      
    if (! is.null(yguess))
      guess <- yguess[1,]
    Y <- y <- guess  # trick
  }
  if (leftbc==length(y))
    stop ("this is not a boundary value problem - use initial value problem solver instead")
  Ynames <- attr(y,"names")
  if (is.null(Ynames) & is.matrix(yguess)) Ynames <- colnames(yguess)
  if (is.null(Ynames) & is.vector(yguess)) Ynames <- names(yguess)

  Func    <- function(x,state)  {
    attr(state,"names") <- Ynames
    func   (x,state,parms,...)[1]
  }
  Func2   <- function(x,state)  {
    attr(state,"names") <- Ynames
    func   (x,state,parms,...)
  }
  if (! is.null(jacfunc))
    JacFunc <- function(x,state) {
      attr(state,"names") <- Ynames
      jacfunc(x,state,parms,...)
    }
  if (! is.null(bound))
    Bound  <- function(i,state)  {
      attr(state,"names") <- Ynames
      bound   (i,state,parms,...)
    }
  if (! is.null(jacbound))
    JacBound   <- function(ii,state)  {
      attr(state,"names") <- Ynames
      jacbound(ii,state,parms,...)
    }

## function evaluation
  tmp <- eval(Func(x[1], y), rho)
  if (!is.list(tmp))
    stop("Model function must return a list\n")
  ncomp  <- length(tmp[[1]])    # number of differential equations

## in case jacobian function is not defined...
  if ( is.null(jacfunc)) {
    JAC      <- matrix(nr=ncomp,nc=ncomp)
    perturbfac  <- 1e-8
    JacFunc <- function (x, state) {
      state2 <- state
      tmp2   <- unlist(eval(Func(x, state), rho))
      for (i in 1:ncomp) {
        dstate   <-max(perturbfac,state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- unlist(eval(Func(x, state), rho))
        JAC[,i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(JAC)
    }
  }
## in case boundary function is not defined...
  if ( is.null(bound)) {
    iini <- sum(!is.na(Y))  
    iend <- sum(!is.na(yend))
    iibb <- c(which( !is.na(Y)),which(!is.na(yend)))
    bb    <- c(Y[!is.na(Y)],yend[!is.na(yend)])
    Bound  <- function(i,state)
          {
    if (is.function(yini))
       {y   <- yini(state,parms,...)
        bb  <- c(y[!is.na(y)],yend[!is.na(yend)])
       }
       return(state[iibb[i]]-bb[i]) }        # too simple for now..
   }       
  if ( is.null(jacbound)) {
    JacBound <- function (ii, state) {
      BJAC        <- numeric(ncomp)
      perturbfac  <- 1e-8
      state2 <- state
      tmp2     <- Bound(ii, state)
      for (i in 1:ncomp) {
        dstate   <-max(perturbfac,state[i]*perturbfac)
        state[i] <- state[i]+ dstate
        tmp      <- Bound(ii, state)
        BJAC[i]  <- (tmp-tmp2)/dstate
        state[i] <- state2[i]
      }
      return(BJAC)
    }
  }

  Nglobal <- if (length(tmp) > 1)
    length(unlist(tmp[-1]))
  else 0
  Nmtot <- attr(unlist(tmp[-1]), "names")

  storage.mode(y) <- storage.mode(x) <- "double"
  linear <- givmesh <- givu <-FALSE
  nmesh  <- 0
  Xguess <- xguess
  Yguess <- yguess
  if (! is.null (xguess))  {
    givmesh <- TRUE
    nmesh   <- length(xguess)
    if (nmesh<length(x)) {
      warning ("expanding xguess; must be provided in as many mesh points as x")
      Xguess <- x
    } else Xguess <- xguess
  }
  if (! is.null(yguess))  {
    if (is.null(xguess))
      stop ("xguess must be provided if yguess is")
    givu <- TRUE
    if (nmesh<length(x))  { # expand yguess
      Yguess <- NULL
      nmesh <- length(Xguess)
      for (i in 1:ncomp)
        Yguess <- rbind(Yguess,approx(xguess,yguess[i,],Xguess)$y)
    } else  Yguess <- yguess # ncomp,nmesh
    if (length(Yguess) != nmesh*ncomp) stop ("xguess and yguess not compatible")
  }

  lwrkfl <- nmax *(4*ncomp*ncomp+12*ncomp+3) + 5 *(ncomp*ncomp)-2*ncomp-9*ncomp
  lwrkin <- nmax*(ncomp+2)+ncomp
  if (length(atol) ==1)
    atol <- rep(atol,len=ncomp)
  else  if (length(atol) != ncomp)
    stop("tol must either be one number or a vector with length=number of state variables")

  out <- .Call("call_colmod",as.integer(ncomp),as.integer(leftbc),
            as.double(fixpt),as.double(aleft),as.double(aright),
            as.double(atol),as.integer(linear),
            as.integer(givmesh),as.integer(givu),as.integer(nmesh),
            as.integer(nmax),as.integer(lwrkfl),as.integer(lwrkin),
            as.double(Xguess), as.double(Yguess),
            Func, JacFunc, Bound, JacBound,
            rho, PACKAGE="bvpSolve")
  nn <- attr(out,"istate")
  mesh <- nn[3]
  attr(out,"istate") <- NULL

  nm <- c("time",
          if (!is.null(attr(y,"names"))) names(y) else as.character(1:ncomp))
  out <- cbind(out[1:mesh],matrix(nr=mesh,out[-(1:mesh)],byrow=TRUE))
  dimnames(out) <- list(NULL,nm)
  out
}

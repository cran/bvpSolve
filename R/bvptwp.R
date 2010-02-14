
##==============================================================================
## Solving boundary value problems of ordinary differential equations
## using MIRK method "twpbvp"
##==============================================================================

bvptwp<- function(yini=NULL, x, func, yend=NULL, parms=NULL, guess=NULL,
     xguess=NULL, yguess=NULL, jacfunc=NULL, bound=NULL, jacbound=NULL,
     leftbc=NULL, islin=FALSE, nmax=1000, atol=1e-8, cond = FALSE,
     allpoints=TRUE, dllname=NULL, initfunc=dllname, ncomp=NULL,
     forcings=NULL, initforc = NULL, fcontrol=NULL, verbose = FALSE, ...)   {

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

  if (!is.function(func) && !is.character(func))
    stop("`func' must be a function or character vector")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
  if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
    stop("`jacfunc' must be a function or character vector")
  if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
    stop("`jacfunc' must be a function or character vector")
  if (!is.null(bound) && !(is.function(bound) || is.character(bound)))
    stop("`bound' must be a function or character vector")
  if (!is.null(jacbound) && !(is.function(jacbound) || is.character(jacbound)))
    stop("`jacbound' must be a function or character vector")

  if (is.character(func)) {
      if (!is.character(jacfunc))
         stop("If 'func' is dynloaded, so must 'jacfunc' be")
      if (!is.character(bound))
         stop("If 'func' is dynloaded, so must 'bound' be")
      if (!is.character(jacbound))
         stop("If 'func' is dynloaded, so must 'jacbound' be")
  }

##---------------------
## Boundary conditions  specified by yini and yend
##---------------------

  if (! is.null(yini))  {    # yini is specified
    ncomp  <- length(yini) 
    y       <- yini
    Y       <- y
    inix    <- which (is.na(y))
    nas     <- length(inix)
    leftbc  <- length(which (!is.na(y)))   # NA when initial condition not known

    if (is.null(guess)&& ! is.null(yguess))
      guess <- yguess[inix,1]

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
    if (leftbc==length(y))
      stop ("this is not a boundary value problem - use initial value problem solver instead")
  } else   {               # boundaries specified by boundary function 'bound'
    if (is.null(leftbc))
      stop("leftbc should be inputted if bound is given")
      
    if (! is.null(yguess))
      guess <- yguess[,1]  #KS: or [,1]
    Y <- y <- guess  # trick
  }
  if (length(y) == 0 && is.null(ncomp))
    stop ("don't know the number of state variables - provide 'ncomp' or 'guess', a initial guess")
  Ynames <- attr(y,"names")
  if (is.null(Ynames) & ! is.null(yend))   Ynames <- names(yend)
  if (is.null(Ynames) & is.matrix(yguess)) Ynames <- colnames(yguess)
  if (is.null(Ynames) & is.vector(yguess)) Ynames <- names(yguess)
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL
  # The functions are in a DLL
  if (is.character(func)) {
    if (sum(duplicated (c(func,initfunc,jacfunc))) >0)
      stop("func, initfunc, or jacfunc cannot be the same")

    if (! is.null(initfunc))  # KS: ADDED THAT to allow absence of initfunc
      if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else if (initfunc != dllname && ! is.null(initfunc))
        stop(paste("'initfunc' not loaded ",initfunc))

    # Easier to deal with NA in C-code
    if (is.null(initfunc)) ModelInit <- NA

    funcname <- func
    ## get the pointer and put it in func
    if (is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
    } else stop(paste("bvp function not loaded",funcname))

    ## the Jacobian
    if (is.loaded(jacfunc, PACKAGE = dllname))  {
      JacFunc <- getNativeSymbolInfo(jacfunc, PACKAGE = dllname)$address
    } else
      stop(paste("jacobian function not loaded ",jacfunc))

    ## the boundary
    if (is.loaded(bound, PACKAGE = dllname))  {
      Bound <- getNativeSymbolInfo(bound, PACKAGE = dllname)$address
    } else
      stop(paste("boundary function not loaded ",bound))

    ## the boundary Jacobian
    if (is.loaded(jacbound, PACKAGE = dllname))  {
      JacBound <- getNativeSymbolInfo(jacbound, PACKAGE = dllname)$address
    } else
      stop(paste("boundary jac function not loaded ",jacbound))
    Nglobal <-  0
    Nmtot <- NULL

    if (! is.null(forcings))
      flist <- checkforcings(forcings,x,dllname,initforc,FALSE,fcontrol)

  } else {      # The functions are R-code
  
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
    if (is.null(y) & ! is.null(ncomp)) y<-runif(ncomp) 
    tmp <- eval(Func(x[1], y), rho)
    if (!is.list(tmp))
      stop("Model function must return a list\n")
   if (! is.null(ncomp)  )  {
    if (length(tmp[[1]])!= length(y))
      stop (" 'func' should return ",paste(ncomp)," elements")
  } else ncomp  <- length(tmp[[1]])    # number of differential equations

## in case jacobian function is not defined...
    if ( is.null(jacfunc)) {
      JAC      <- matrix(nrow=ncomp, ncol=ncomp)
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
      iibb <- c(which( !is.na(Y)),which(!is.na(yend)))
      bb    <- c(Y[!is.na(Y)],yend[!is.na(yend)])
      Bound  <- function(i,state)        {
         return(state[iibb[i]]-bb[i])
      }         
    }  else {
     tmp <- eval(bound(1, y), rho)
     if (length(tmp) > 1)
       stop ("function 'bound' should return only ONE value")
    }
# Jacobian of the boundary function
    if ( is.null(jacbound)) {
      if (is.null(bound)) {
        JMat <- matrix(data=0, nrow=ncomp, ncol = ncomp)
        iibb <- c(which( !is.na(Y)),which(!is.na(yend)))
        JMat[cbind(1:ncomp,iibb)] <- 1
        JacBound <- function (ii, state) 
           return(JMat[ii,])
        }
      else 
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
  }  # ! is.character func

  storage.mode(y) <- storage.mode(x) <- "double"
  linear <- givmesh <- givu <-FALSE
  nmesh  <- 0
  Xguess <- xguess
  Yguess <- yguess
# check dimensions
  if (! is.null(yguess) && ! is.null(yguess))
     if (length(xguess) != ncol(yguess))
       stop("yguess should have as many columns as length of xguess")
       
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
  ntol <- ncomp
  lwrkfl <- nmax *(6*ncomp*ncomp+22*ncomp+3) + 6 *(ncomp*ncomp)+22*ncomp+2*ntol
  lwrkin <- nmax*(2*ncomp+3)+2*ncomp
  if (length(atol) ==1)
    atol <- rep(atol,len=ncomp)
  else  if (length(atol) != ncomp)
    stop("tol must either be one number or a vector with length=number of state variables")
  Ipar <- 1
  Rpar <- 1.0
    if(is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
    else
      initpar <- as.double(parms)
  if(verbose) ow <- options("warn"=1)   # print as warnings are generated
      
  out <- .Call("call_bvptwp",as.integer(ncomp),as.integer(leftbc),
            as.double(fixpt),as.double(aleft),as.double(aright),
            as.double(atol),as.integer(linear), as.integer(verbose),
            as.integer(givmesh),as.integer(givu),as.integer(nmesh),
            as.integer(nmax),as.integer(lwrkfl),as.integer(lwrkin),
            as.double(Xguess), as.double(Yguess),
            as.double(Rpar), as.integer(Ipar), as.integer(cond),
            Func, JacFunc, Bound, JacBound, ModelInit, initpar,
            flist, rho, PACKAGE="bvpSolve")
  nn <- attr(out,"istate")
  rn <- attr(out,"rstate")
  mesh <- nn[3]
  attr(out,"istate") <- NULL

#(  if(verbose) print(names(baseenv()$last.warning))  # dangerous...
  if(verbose) options(ow)     # reset printing options

  nm <- c("x",
          if (!is.null(attr(y,"names"))) names(y) else as.character(1:ncomp))
  out <- cbind(out[1:mesh],matrix(data=out[-(1:mesh)],nrow=mesh,byrow=TRUE))
  # select only the rows corresponding to x-values
  if (! allpoints)
    if (nrow(out) > length(x))
      out <- out [which(out[,1]%in% x),]
  dimnames(out) <- list(NULL,nm)
  class(out) <- c("bvpSolve","matrix")  # a boundary value problem
  names(nn) <- c("flag","nmax","nmesh","nrwork","niwork")
  attr(out,"istate") <- nn 
  names(rn) <- c("ckappa1","gamma1","sigma","ckappa","ckappa2")
  attr(out,"rstate") <- rn 
  attr(out,"name") <- "bvptwp"
  out
}

 
 

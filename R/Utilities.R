### ============================================================================
### S3 methods
### ============================================================================

## An approximation function for bvpSolve objects that were solved with \
## bvpcol
approx <- function (x, ...) UseMethod("approx")

approx.default <- function (x, ...) {
if ("bvpSolve" %in% class (x))
  approx.bvpSolve(x,...)
else  
  stats::approx(x,...)
#nextmethod()  
}

### ============================================================================

approx.bvpSolve <- function(x, xout=NULL, ...){
 
  Attr <- attributes(x)
  if (Attr$name != "bvpcol")
    stop("can only use 'approx.bvpSolve' if problem was solved with 'bvpcol'")
  istate <- Attr$istate[-(1:6)]   # first 6 elements have nothing to do with continuation
  rstate <- Attr$rstate
  il <- length(xout)
  if (il <= 0)
    stop ("'approx' requires at least one value to approximate")   

  z <- rep(1, istate[4])
  appone <- function(x)
    .Fortran("appsln", as.double(x), 
            result = as.double(z), as.double(rstate), as.integer(istate))$result
  Out <- NULL
  for (i in 1:il)
    Out <- rbind(Out,c(xout[i],appone(xout[i])))
  colnames(Out) <- colnames(x)
  class (Out) <- c( "bvpSolve", "matrix" ) 
  attr(Out, "name") <-  "approx"
  dimnames(Out) <- dimnames(x)
  Out
}

### ============================================================================

print.bvpSolve <- function(x, ...)
   print(as.data.frame(x), ... )

### ============================================================================

plot.bvpSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
    t <- 1     # column with "times"
    var <- colnames(x)
    if (!is.numeric(which)) {
        ln <- length(which)
        Select <- which(var %in% which)
        if (length(Select) != ln)
            stop("not all variables in 'which' are in 'x'")
        which <- Select
    }
    else {
        which <- which + 1  # "which now refers to the column number
        if (max(which) > ncol(x))
            stop("index in 'which' too large")
        if (min(which) < 1)
            stop("index in 'which' should be > 0")
    }
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
        nc <- min(ceiling(sqrt(np)), 3)  # assume max 3 x 3 panels
        nr <- min(ceiling(np/nc), 3)
        mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else
        mfrow <- dots$mfrow

    if (!is.null(mfrow)) {
        mf <- par(mfrow = mfrow)
    }
    if (is.null(ask))
        ask <- prod(par("mfcol")) < length(which) && dev.interactive()
    ## interactively wait if there are remaining figures
    if (ask) {
        oask <- devAskNewPage(TRUE)
	      on.exit(devAskNewPage(oask))
    }

    Main <- is.null(dots$main)

    xxlab <- if (is.null(dots$xlab))  colnames(x)[t]  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""              else dots$ylab
    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        if (Main)
            dots$main <- colnames(x)[i]
            dots$xlab <- xxlab[i-1]
            dots$ylab <- yylab[i-1]
        do.call("plot", c(alist(x[, t], x[, i]), dots))
    }
}

### ============================================================================

diagnostics.bvpSolve<- function(obj, ...) {
  if (!"bvpSolve" %in% class(obj)) return(NULL)
  Attr <- attributes(obj)
  istate <- Attr$istate
  rstate <- Attr$rstate
  idid <- istate[1]

  if (is.null(istate) || is.null (rstate)) return(NULL)
  cat("\n--------------------\n")
  cat(paste( "solved with ",Attr$name))
  cat("\n--------------------\n")

  if (Attr$name == "bvpshoot") {
    cat("\n---------------------------------------------\n")
    cat("diagnostics of BVP solver ")
    cat("\n---------------------------------------------\n")
    df <- c( "The number of function evaluations              :", 
             "The number of jacobian evaluations +LU decomp   :",	
             "The number of steps                             :",
             "The number of calls to the ivp solver           :")
    printmessage(df, Attr$istate2)

    cat("\n---------------------------------------------\n")
    cat("diagnostics of the last run of the IVP solver ")
    cat("\n---------------------------------------------\n")
      diagnostics.deSolve(obj)
    
  } else if (Attr$name == "bvptwp") {
  
    if (idid ==0)  cat("  Integration was successful.\n") else
       cat("  Integration was NOT successful\n")
    df <- c( "The return code                             :",   #1
             "The number of function evaluations          :", 
             "The number of jacobian evaluations          :",	
             "The number of steps                         :", 
             "The number of boundary evaluations          :", 
             "The number of boundary jacobian evaluations :", 
             "The number of mesh resets                   :",
             "The maximal number of mesh points           :",   #2
             "The actual number of mesh points            :",
             "The size of the real work array             :",
             "The size of the integer work array          :"
             )


    printmessage(df, istate)

    cat("\n--------------------\n")
    cat(paste( "conditioning pars"))
    cat("\n--------------------\n")

    df <- c( "ckappa1  :",
             "gamma1   :",
             "sigma    :",
             "ckappa   :",
             "ckappa2  :")
    printmessage(df, rstate)
    
    
  } else if (Attr$name == "bvpcol") {
    if (idid ==1)  cat("  Integration was successful.\n") else
       cat("  Integration was NOT successful\n")
    df <- c( "The return code                                   :",   #1
             "The number of function evaluations                :", 
             "The number of jacobian evaluations                :",	
             "The number of steps                               :", 
             "The number of boundary evaluations                :", 
             "The number of boundary jacobian evaluations       :", 
             "The actual number of mesh points                  :",
             "The number of collocation points per subinterval  :",
             "The number of equations                           :",
             "The number of components (variables)              :",
             rep(
             "The order of each equation                        :",istate[4+5]))


    printmessage(df, istate[-c(11:13)])
    }    

}

### ============================================================================
## internal helper functions for printing solver return code messages
## these functions are not exported

printmessage <-function(message1, state, message2 = NULL, Nr = 1:length(message1)) {
  if (is.null(message2)) {
    cat("\n", paste(formatC(Nr, "##", width = 2), message1,
              signif(state, digits = getOption("digits")), "\n"), "\n")
  } else {
    cat("\n", paste(formatC(Nr, "##", width = 2), message1,
              signif(state, digits = getOption("digits")), message2, "\n"), "\n")

  }
}



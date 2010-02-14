### ============================================================================
### S3 methods
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

diagnostics.bvpSolve<- function(obj, ...) {
    if (!"bvpSolve" %in% class(obj)) return(NULL)
    Attr <- attributes(obj)
    istate <- Attr$istate
    rstate <- Attr$rstate

    if (is.null(istate) || is.null (rstate)) return(NULL)
    if (! Attr$name == "bvptwp") return(NULL)

    cat("\n--------------------\n")
    cat(paste( "return code"))
    cat("\n--------------------\n")

    idid <- istate[1]
    if (idid ==0)  cat("  Integration was successful.\n") else
       cat("  Integration was NOT successful\n")

    df <- c( "The return code                    :",   #1
             "The maximal number of mesh points  :",   #2
             "The actual number of mesh points   :",
             "The size of the real work array    :",
             "The size of the integer work array :")


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

}

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



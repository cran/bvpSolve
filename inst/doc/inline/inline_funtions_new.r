
## this is copied from R-package 'inline':
## as it is not exported, it is copied

compilecode <- function (f, code, args = c("ncomp", "x", "y", "F", "rpar", "ipar"),
    type = c("integer", rep("double", 4), "integer"), language, verbose    ) 
{
    wd <- getwd()
    on.exit(setwd(wd))
    if (.Platform$OS.type == "windows") {
        dir <- gsub("\\\\", "/", tempdir())
        libCFile <- paste(dir, "/", f, ".EXT", sep = "")
        libLFile <- paste(dir, "/", f, ".dll", sep = "")
        libLFile2 <- paste(dir, "/", f, ".dll", sep = "")
    }
    else {
        libCFile <- paste(tempdir(), "/", f, ".EXT", sep = "")
        libLFile <- paste(tempdir(), "/", f, .Platform$dynlib.ext, 
            sep = "")
        libLFile2 <- paste(tempdir(), "/", f, ".sl", sep = "")
    }
    extension <- switch(language, `C++` = ".cpp", C = ".c", Fortran = ".f", 
        F95 = ".f95", ObjectiveC = ".m", `ObjectiveC++` = ".mm")
    libCFile <- sub(".EXT$", extension, libCFile)
    write(code, libCFile)
    if (file.exists(libLFile)) 
        file.remove(libLFile)
    if (file.exists(libLFile2)) 
        file.remove(libLFile2)
    setwd(dirname(libCFile))
    errfile <- paste(basename(libCFile), ".err.txt", sep = "")
    cmd <- paste(R.home(component = "bin"), "/R CMD SHLIB ", 
        basename(libCFile), " 2> ", errfile, sep = "")
    if (verbose) 
        cat("Compilation argument:\n", cmd, "\n")
    compiled <- system(cmd, intern = !verbose)
    errmsg <- readLines(errfile)
    unlink(errfile)
    writeLines(errmsg)
    setwd(wd)
    if (!file.exists(libLFile) && file.exists(libLFile2)) 
        libLFile <- libLFile2
    if (!file.exists(libLFile)) {
        cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
        cat("\nProgram source:\n")
        code <- strsplit(code, "\n")
        for (i in 1:length(code[[1]])) cat(format(i, width = 3), 
            ": ", code[[1]][i], "\n", sep = "")
        stop(paste("Compilation ERROR, function(s)/method(s) not created!", 
            paste(errmsg, collapse = "\n")))
    }

    cleanup <- function(env) {
        if (f %in% names(getLoadedDLLs())) 
            dyn.unload(libLFile)
        unlink(libLFile)
    }
    reg.finalizer(environment(), cleanup, onexit = TRUE)
    DLL <- dyn.load(libLFile)
    fn <- function(arg) NULL

    F <- formals(fn)[rep(1, length(args))] 
    names(F) <- args
    formals(fn) <- F
    body <- quote(CONVENTION("EXTERNALNAME", ARG))[c(1:2, 
                rep(3, length(args)))]
    for (j in seq(along = args)) body[[j+2]] <-  as.name(type)
    body[[1]] <- get(convention)
    body[[2]] <- getNativeSymbolInfo(f, DLL)$address
    body(fn) <- body
    res <- new("CFunc", code = code)
    res@.Data <- fn

    if (verbose) {
        cat("Program source:\n")
        lines <- strsplit(code, "\n")
        for (i in 1:length(lines[[1]])) cat(format(i, width = 3), 
            ": ", lines[[1]][i], "\n", sep = "")
    }
    remove(list = c("F", "body", "fn"))
    return(res)
}


cfsub <- function (fsub, header = NULL, tail = NULL, language = "F95") {
  if (! is.character(fsub))
    stop ("'fsub' should be a character string ")
    
  if (! is.null(header))
    if (! is.character(header))
    stop ("'header' should be a character string or NULL")

  if (! is.null(tail))
    if (! is.character(tail))
    stop ("'tail' should be a character string or NULL")

  convention <- ".Fortran"  
  sig <- signature(XX = "integer", xvec = "numeric", 
    y = "numeric", F = "numeric", rpar = "numeric", ipar = "integer")
    
  signature <-"ncomp, x, y, F, rpar, ipar"
  f <- basename(tempfile())
  if (language == "Fortran") {
    head <- paste("       SUBROUTINE", f, "(", signature, ")\n        IMPLICIT NONE")
    head <- paste(head, "       INTEGER ncomp\n       DOUBLE PRECISION x, y(*), f(*), rpar(*)\n       integer ipar(*)\n")  
    fsub <- paste(head, fsub, tail, "\n       END SUBROUTINE")
    
  } else if (language == "F95") {
    head <- paste("subroutine", f, "(", signature, ")\n implicit none")
    head <- paste(head, " integer ncomp\n double precision x, y(*), f(*), rpar(*)\n integer ipar(*)\n")  
    fsub <- paste(head, fsub, tail, "\n end subroutine")
    
  } else  {
#    head <- "#include <R.h>\n#include <Rdefines.h>\n#include <R_ext/Error.h>\n"
    head <- paste(header, "\nvoid", f, "(", signature, ")\n{")
    fsub <- paste(head, fsub, tail, "\n }\n")
  }
  
  libLFile <- compilecode(f, fsub, language, TRUE)
  
  liblFile
}


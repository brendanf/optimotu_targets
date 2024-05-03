filterAndTrim <- function(fwd, filt, rev, filt.rev, ...) {
  if (length(fwd) > 0L) {
    # ensure files are created, even if they end up being empty
    file.create(c(filt, filt.rev))
    mycall <- match.call()
    mycall[[1]] <- dada2::filterAndTrim
    eval.parent(mycall)
    # return file names for samples where at least some reads passed
    purrr::keep(c(filt, filt.rev), file.exists)
  } else {
    character()
  }
}

derepFastq <- function(fls, n = 1e+06, verbose = FALSE, qualityType = "Auto",
                       names = fls) {
  if (length(fls) == 0) {
    out <- list()
  } else {
    out <- dada2::derepFastq(fls, n = n, verbose = verbose, qualityType = qualityType)
    if (methods::is(out, "derep")) {
      out <- list(out)
    }
    names(out) <- names
  }
  out
}

learnErrors <- function(fls, nbases = 1e+08, nreads = NULL,
                        errorEstimationFunction = dada2::loessErrfun,
                        multithread = FALSE, randomize = FALSE, MAX_CONSIST = 10,
                        OMEGA_C = 0, qualityType = "Auto", verbose = FALSE, ...) {
  if (length(fls) == 0) {
    NULL
  } else {
    mycall = match.call()
    mycall[[1]] <- dada2::learnErrors
    eval.parent(mycall)
  }
}

dada <- function(
    derep,
    err,
    errorEstimationFunction = dada2::loessErrfun,
    selfConsist = FALSE,
    pool = FALSE,
    priors = character(0),
    multithread = FALSE,
    verbose = TRUE,
    ...
) {
  if (is.null(err) && isFALSE(selfConsist)) {
    NULL
  } else if (length(derep) == 0) {
    list()
  } else {
    mycall <- match.call()
    mycall[[1]] <- dada2::dada
    out <- eval.parent(mycall)
    if (methods::is(out, "dada") && !methods::is(derep, "derep")) {
      out <- list(out)
      names(out) <- names(derep)
    }
    out
  }
}

mergePairs <- function(dadaF, derepF, dadaR, derepR,
                       minOverlap = 12,
                       maxMismatch = 0,
                       returnrejects = FALSE,
                       propagateCol = character(0),
                       justConcatenate = FALSE,
                       trimOverhang = FALSE,
                       verbose = FALSE,
                       ...) {
  if (length(dadaF) == 0 || length(dadaR) == 0) {
    list()
  } else {
    mycall <- match.call()
    mycall[[1]] <- dada2::mergePairs
    out <- eval.parent(mycall)
    if (is.data.frame(out) && !methods::is(dadaF, "dada")) {
      out <- list(out)
      names(out) <- names(dadaF)
    }
    out
  }
}

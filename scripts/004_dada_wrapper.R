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


# From dada2 commit https://github.com/benjjneb/dada2/commit/7714487b153ca133cb6f03cb01d09fc05be60159
makeBinnedQualErrfun <- function(binnedQ=c(2, 11, 25, 37)) {
  function(trans, binnedQuals=binnedQ) {
    qq <- as.numeric(colnames(trans))
    # Get min and max observed quality scores
    qmax <- max(qq[colSums(trans)>0])
    qmin <- min(qq[colSums(trans)>0])
    # Check for data consistency with provided binned qualities
    if(qmax > max(binnedQuals)) stop("Input data contains a higher quality score than the provided binned values.")
    if(qmin < min(binnedQuals)) stop("Input data contains a lower quality score than the provided binned values.")
    if(!qmax %in% binnedQuals) warning("Maximum observed quality score is not in the provided binned values.")
    if(!qmin %in% binnedQuals) warning("Minimum observed quality score is not in the provided binned values.")

    est <- matrix(0, nrow=0, ncol=length(qq))
    for(nti in c("A","C","G","T")) {
      for(ntj in c("A","C","G","T")) {
        if(nti != ntj) {
          errs <- trans[paste0(nti,"2",ntj),]
          tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
          p <- errs/tot
          df <- data.frame(q=qq, errs=errs, tot=tot, p=p)
          # Check and enforce that this q scores start at zero
          if(!all(df$q == seq(nrow(df))-1)) stop("Unexpected Q score series.") ###!
          pred <- rep(NA, nrow(df))
          for(i in seq(length(binnedQuals)-1)) {
            loQ <- binnedQuals[i]
            hiQ <- binnedQuals[i+1]
            loP <- df$p[loQ+1]
            hiP <- df$p[hiQ+1]
            # Linear interpolation between the binned Q scores observed in the data
            if(!is.na(loP) && !is.na(hiP)) {
              pred[(loQ+1):(hiQ+1)] <- exp(seq(log(loP), log(hiP), length.out=(hiQ-loQ+1)))
            }
          }

          maxrli <- max(which(!is.na(pred)))
          minrli <- min(which(!is.na(pred)))
          pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
          pred[seq_along(pred)<minrli] <- pred[[minrli]]
          est <- rbind(est, pred)
        } # if(nti != ntj)
      } # for(ntj in c("A","C","G","T"))
    } # for(nti in c("A","C","G","T"))

    # HACKY
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE

    # Expand the err matrix with the self-transition probs
    err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                 est[4,], 1-colSums(est[4:6,]), est[5:6,],
                 est[7:8,], 1-colSums(est[7:9,]), est[9,],
                 est[10:12,], 1-colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
    colnames(err) <- colnames(trans)
    # Return
    return(err)
  }
}

choose_dada_error_function <- function(fls, ...) {
  bins <- optimotu.pipeline:::fastq_qual_bins(fls, ...)
  if (length(bins) < 10) {
    makeBinnedQualErrfun(bins)
  } else {
    dada2::loessErrFun
  }
}

learnErrors <- function(fls, nbases = 1e+09, nreads = NULL,
                        errorEstimationFunction = choose_dada_error_function,
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
    errorEstimationFunction = choose_dada_error_function,
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

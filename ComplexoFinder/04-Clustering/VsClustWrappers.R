library(data.table)
library(vsclust)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(parallel)

# PrepareForVsClust
PrepareForVsClust_custom <- function(
  dat,
  NumReps,
  NumCond,
  isPaired = FALSE,
  isStat
) {
  qvals <- statFileOut <- Sds <- NULL
  tdat <- NULL

  # convert to matrix
  dat <- as.matrix(dat)

  # Run statistical testing
  if (isStat) {
    if (ncol(dat) != NumReps * NumCond) {
      stop(
        "Number of data columns must correspond to product of conditions
             and replicates!"
      )
    }
    if (isPaired) {
      ttt <- SignAnalysisPaired(dat, NumCond, NumReps)
    } else {
      ttt <- SignAnalysis(dat, NumCond, NumReps)
    }

    Sds <- ttt$Sds
    qvals <- ttt$qvalues

    # Generate extended letter labels for more than 26 conditions
    extended_letters <- c(LETTERS, paste0("A", LETTERS), paste0("B", LETTERS))
    colnames(qvals) <-
      paste("qvalue ", extended_letters[2:(NumCond)], "vsA", sep = "")

    tdat <- averageCond(dat, NumReps, NumCond)
  } else {
    Sds <- dat[, ncol(dat)]
    tdat <- dat[, seq_len(ncol(dat) - 1)]
    NumReps <- 1
    NumCond <- ncol(dat) - 1
    dat <- tdat
  }

  if (isStat) {
    statFileOut <- cbind(tdat, Sds, qvals)
  } else {
    statFileOut <- cbind(tdat, Sds)
  }

  ## Preparing output
  Out <-
    list(dat = cbind(tdat, Sds), qvals = qvals, statFileOut = statFileOut)
  Out
}

# RunClustWrapper
runClustWrapper_custom <- function(
  dat,
  NClust,
  proteins = NULL,
  VSClust = TRUE,
  scaling = "standardize",
  constraints = NULL,
  cores,
  verbose = FALSE
) {
  tData <- dat[, seq_len(ncol(dat) - 1)]
  sds <- dat[, ncol(dat)]

  #Standardize
  # scale standard deviations by the ones in the actual data to cope for the
  # following standardization
  if (!any(scaling == c("standardize", "center", "none"))) {
    stop("parameter scaling needs to be standardize, center or none!")
  }
  if ((scaling == "standardize")) {
    sds <- sds / rowSds(as.matrix(tData), na.rm = TRUE)
  }
  tData <- t(scale(
    t(tData),
    center = (scaling != "none"),
    scale = (scaling == "standardize")
  ))
  if (is.null(rownames(tData))) {
    rownames(tData) <- seq_len(nrow(tData))
  }
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(
    cl = cl,
    varlist = c("ClustComp", "vsclust_algorithm"),
    envir = environment()
  )

  clustout <- ClustComp(
    tData,
    NClust = NClust,
    Sds = sds,
    constraints = constraints,
    NSs = 16,
    cl = cl,
    verbose = verbose
  )

  if (VSClust) {
    Bestcl <- clustout$Bestcl
  } else {
    Bestcl <- clustout$Bestcl2
  }

  # Handle single cluster case - all peptides in one cluster
  if (NClust == 1) {
    # Set all cluster assignments to 1
    Bestcl$cluster <- rep(1, nrow(tData))
    names(Bestcl$cluster) <- rownames(tData)

    # Set all memberships to 1 (full membership in single cluster)
    Bestcl$membership <- matrix(1, nrow = nrow(tData), ncol = 1)
    rownames(Bestcl$membership) <- rownames(tData)
    colnames(Bestcl$membership) <- "membership of cluster 1"

    # Calculate center as mean of all data
    Bestcl$centers <- matrix(colMeans(tData, na.rm = TRUE), nrow = 1)
    rownames(Bestcl$centers) <- "Cluster 1"
    colnames(Bestcl$centers) <- colnames(tData)

    # Create output
    outFileClust <- tData
    if (!is.null(proteins)) {
      outFileClust <- cbind(
        outFileClust,
        names = as.character(proteins[rownames(outFileClust)])
      )
    }

    ClustInd <- data.frame(
      Cluster = "1",
      Members = nrow(tData)
    )

    Out <- list(
      dat = tData,
      Bestcl = Bestcl,
      outFileClust = outFileClust,
      ClustInd = ClustInd
    )
    return(Out)
  }

  # For NClust > 1, do normal processing
  Bestcl <- SwitchOrder(Bestcl, NClust)

  # sorting for membership values (globally)
  Bestcl$cluster <-
    Bestcl$cluster[order(rowMaxs(Bestcl$membership, na.rm = TRUE))]
  Bestcl$membership <-
    Bestcl$membership[order(rowMaxs(Bestcl$membership, na.rm = TRUE)), ]
  tData <- tData[names(Bestcl$cluster), ]

  colnames(Bestcl$membership) <-
    paste("membership of cluster", colnames(Bestcl$membership))
  outFileClust <- tData
  if (!is.null(proteins)) {
    outFileClust <-
      cbind(
        outFileClust,
        names = as.character(proteins[rownames(outFileClust)])
      )
  }

  rownames(Bestcl$centers) <-
    paste("Cluster", rownames(Bestcl$centers))
  ClustInd <-
    as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership) > 0.5]))
  if (ncol(ClustInd) == 2) {
    colnames(ClustInd) <- c("Cluster", "Members")
  } else {
    ClustInd <-
      cbind(seq_len(max(Bestcl$cluster)), rep(0, max(Bestcl$cluster)))
  }

  ## Output
  Out <-
    list(
      dat = tData,
      Bestcl = Bestcl,
      outFileClust = outFileClust,
      ClustInd = ClustInd
    )
  return(Out)
}

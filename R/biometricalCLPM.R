# ==============================================================================
# biometricalCLPM.R
#
# Reusable OpenMx builders for the Biometrical Cross-Lagged Panel Model (CLPM).
#
# Singh M., Hunter M.D., Assary E., Verhulst B., Peterson R.E., Maes H.H.M.,
# Dolan C.V., Eley T.C., & Neale M.C. (2025). Causation Between Smoking Quantity
# and Depressive Symptoms in Young Adults: Evidence from Novel Cross-Lagged Twin
# Models. medRxiv. https://doi.org/10.1101/2025.11.18.25340516
#
# These functions build exactly the models specified in
#   Model/Biometrical_CLPM_Model_MultipleTimePoints.R  and
#   TEDS/4_*.R, TEDS/5_*.R
# and use identical parameter labels, so results are directly comparable with
# the published analysis scripts.
#
# Requires: OpenMx (>= 2.21.0). No other dependencies.
#
# Usage:  source("R/biometricalCLPM.R")
# ==============================================================================

if (!requireNamespace("OpenMx", quietly = TRUE)) {
  stop("The Biometrical CLPM functions require the OpenMx package: ",
       "install.packages('OpenMx')")
}
library(OpenMx)


# ------------------------------------------------------------------------------
# Naming conventions
# ------------------------------------------------------------------------------
#
# Two phenotypes, x and y, measured at Ts occasions indexed by `waveLabels`
# (integers, not necessarily consecutive: the two-wave TEDS analysis uses
# waveLabels = c(1, 6) so that labels match the six-wave analysis).
#
# Manifest variables : paste0(xVar, w) and paste0(yVar, w) for w in waveLabels,
#                      suffixed by twinSuffix (default "_Tw1" / "_Tw2").
#
# Parameter labels (s and t are consecutive wave labels, s before t):
#
#   Phenotypic regressions (matrix BE)
#     AR_bx{t}x{s}   phenotypic autoregression of x
#     AR_by{t}y{s}   phenotypic autoregression of y
#     CL_by{t}x{s}   distal (lagged) effect  x_s -> y_t
#     CL_bx{t}y{s}   distal (lagged) effect  y_s -> x_t
#     Inst_by{w}x{w} proximal (cross-sectional) effect  x_w -> y_w
#     Inst_bx{w}y{w} proximal (cross-sectional) effect  y_w -> x_w
#
#   Innovation (co)variances of the latent A / C / E factors (matrices PSA/PSC/PSE)
#     VAx{w}, VAy{w}, covAxy{w}      (likewise VC*, covCxy*, VE*, covExy*)
#
#   Autoregressions and cross-lags of the latent factors (matrices BEA/BEC/BEE)
#     AR_Ax{t}x{s}, AR_Ay{t}y{s}     (likewise AR_C*, AR_E*)
#     CL_Ay{t}x{s}, CL_Ax{t}y{s}     (likewise CL_C*, CL_E*)
#
#   Means / intercepts (matrix b0)
#     b0_{manifest name}
#
#   Thresholds for ordinal phenotypes (matrix threG)
#     th{k}_{manifest name}
#
# `clpmLabels()` returns all of these programmatically.
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Default start values for the Biometrical CLPM
#'
#' @param vax,vay Innovation variance of the A factors of x and y.
#' @param vcx,vcy Innovation variance of the C factors of x and y.
#' @param vex,vey Innovation variance of the E factors of x and y.
#' @param raxy,rcxy,rexy Innovation correlations between x and y within A, C, E.
#' @param arPheno Phenotypic autoregression.
#' @param clPheno Phenotypic cross-lag.
#' @param arLatent Start value used when freeing latent (A/C/E) autoregressions.
#' @param clLatent Start value used when freeing latent (A/C/E) cross-lags.
#' @param instPheno Start value used when freeing proximal (DoC) effects.
#' @param thresholds Start values for ordinal thresholds (length `nThresholds`).
#' @return A named list of start values, passed to `biometricalCLPM(start = )`.
# ------------------------------------------------------------------------------
clpmStart <- function(vax = 0.4, vay = 0.4, raxy = 0.2,
                      vcx = 0.2, vcy = 0.2, rcxy = 0.2,
                      vex = 0.4, vey = 0.4, rexy = 0.2,
                      arPheno = 0.6, clPheno = 0.1,
                      arLatent = 0.5, clLatent = 0.1,
                      instPheno = 0.2,
                      mean = 0,
                      thresholds = c(0.5, 1)) {
  list(vax = vax, vay = vay, raxy = raxy,
       vcx = vcx, vcy = vcy, rcxy = rcxy,
       vex = vex, vey = vey, rexy = rexy,
       covaxy = sqrt(vax * vay) * raxy,
       covcxy = sqrt(vcx * vcy) * rcxy,
       covexy = sqrt(vex * vey) * rexy,
       arPheno = arPheno, clPheno = clPheno,
       arLatent = arLatent, clLatent = clLatent,
       instPheno = instPheno,
       mean = mean, thresholds = thresholds)
}


# ------------------------------------------------------------------------------
#' Manifest variable names for a Biometrical CLPM
#'
#' @param xVar,yVar Base names of the two phenotypes (the wave label is pasted
#'   directly onto these, e.g. xVar = "DepSxResRN_T" gives "DepSxResRN_T1").
#' @param waveLabels Integer vector of wave labels, length Ts.
#' @param twinSuffix Length-2 character vector of twin-1 / twin-2 suffixes.
#' @return A list with `labs` (one row of the twin pair, x and y interleaved by
#'   wave), `selVars` (both twins), `x`, `y`, `waveLabels`, `Ts`.
# ------------------------------------------------------------------------------
clpmVarNames <- function(xVar = "X", yVar = "Y", waveLabels = 1:2,
                         twinSuffix = c("_Tw1", "_Tw2")) {
  stopifnot(length(twinSuffix) == 2L, length(waveLabels) >= 2L)
  xn <- paste0(xVar, waveLabels)
  yn <- paste0(yVar, waveLabels)
  labs <- as.vector(rbind(xn, yn))            # x1, y1, x2, y2, ...
  list(labs       = labs,
       selVars    = c(paste0(labs, twinSuffix[1]), paste0(labs, twinSuffix[2])),
       x          = xn,
       y          = yn,
       waveLabels = waveLabels,
       Ts         = length(waveLabels),
       twinSuffix = twinSuffix)
}


# ------------------------------------------------------------------------------
#' All parameter labels used by the Biometrical CLPM
#'
#' Useful for `omxSetParameters()`, `mxCI()`, and for reading the output tables.
#'
#' @param waveLabels Integer vector of wave labels.
#' @return A nested list: `$pheno$ar`, `$pheno$cl`, `$pheno$inst`, and for each
#'   component in A/C/E: `$A$var`, `$A$cov`, `$A$ar`, `$A$cl`.
# ------------------------------------------------------------------------------
clpmLabels <- function(waveLabels = 1:2) {
  w  <- waveLabels
  ws <- utils::head(w, -1L)   # "from" waves
  wt <- utils::tail(w, -1L)   # "to" waves

  comp <- function(k) list(
    var = c(paste0("V", k, "x", w), paste0("V", k, "y", w)),
    cov = paste0("cov", k, "xy", w),
    ar  = c(paste0("AR_", k, "x", wt, "x", ws), paste0("AR_", k, "y", wt, "y", ws)),
    cl  = c(paste0("CL_", k, "y", wt, "x", ws), paste0("CL_", k, "x", wt, "y", ws))
  )

  list(
    pheno = list(
      ar   = c(paste0("AR_bx", wt, "x", ws), paste0("AR_by", wt, "y", ws)),
      cl   = c(paste0("CL_by", wt, "x", ws), paste0("CL_bx", wt, "y", ws)),
      clyx = paste0("CL_by", wt, "x", ws),   # x -> y  (distal)
      clxy = paste0("CL_bx", wt, "y", ws),   # y -> x  (distal)
      inst = c(paste0("Inst_by", w, "x", w), paste0("Inst_bx", w, "y", w)),
      instyx = paste0("Inst_by", w, "x", w), # x -> y  (proximal)
      instxy = paste0("Inst_bx", w, "y", w)  # y -> x  (proximal)
    ),
    A = comp("A"), C = comp("C"), E = comp("E"),
    waveLabels = w
  )
}


# ---- internal helpers --------------------------------------------------------

.clpmPsi <- function(vn, comp, start) {
  # Innovation (co)variance matrix of the latent A, C or E factors.
  nv   <- 2L * vn$Ts
  dn   <- paste0(comp, vn$labs)
  vals <- matrix(0, nv, nv, dimnames = list(dn, dn))
  labs <- free <- vals
  labs[] <- NA_character_
  free[] <- FALSE
  storage.mode(free) <- "logical"

  vx  <- start[[paste0("v", tolower(comp), "x")]]
  vy  <- start[[paste0("v", tolower(comp), "y")]]
  cxy <- start[[paste0("cov", tolower(comp), "xy")]]

  for (i in seq_along(vn$waveLabels)) {
    w  <- vn$waveLabels[i]
    ax <- paste0(comp, vn$x[i]); ay <- paste0(comp, vn$y[i])
    vals[ax, ax] <- vx;  vals[ay, ay] <- vy
    vals[ax, ay] <- vals[ay, ax] <- cxy
    labs[ax, ax] <- paste0("V", comp, "x", w)
    labs[ay, ay] <- paste0("V", comp, "y", w)
    labs[ax, ay] <- labs[ay, ax] <- paste0("cov", comp, "xy", w)
    free[ax, ax] <- free[ay, ay] <- free[ax, ay] <- free[ay, ax] <- TRUE
  }
  mxMatrix(type = "Full", nrow = nv, ncol = nv, free = free, values = vals,
           labels = labs, dimnames = list(dn, dn), name = paste0("PS", comp))
}

.clpmLatentBeta <- function(vn, comp) {
  # Autoregression / cross-lag matrix of the latent A, C or E factors.
  # Built with labels but all paths fixed at 0; freed by useBiometricalAR() and
  # useLiabilityCrossLags().
  nv <- 2L * vn$Ts
  dn <- paste0(comp, vn$labs)
  labs <- matrix(NA_character_, nv, nv, dimnames = list(dn, dn))

  for (i in seq_len(vn$Ts - 1L)) {
    s <- vn$waveLabels[i]; t <- vn$waveLabels[i + 1L]
    xs <- paste0(comp, vn$x[i]);      ys <- paste0(comp, vn$y[i])
    xt <- paste0(comp, vn$x[i + 1L]); yt <- paste0(comp, vn$y[i + 1L])
    labs[xt, xs] <- paste0("AR_", comp, "x", t, "x", s)
    labs[yt, ys] <- paste0("AR_", comp, "y", t, "y", s)
    labs[yt, xs] <- paste0("CL_", comp, "y", t, "x", s)
    labs[xt, ys] <- paste0("CL_", comp, "x", t, "y", s)
  }
  mxMatrix(type = "Full", nrow = nv, ncol = nv, free = FALSE, values = 0,
           labels = labs, dimnames = list(dn, dn), name = paste0("BE", comp))
}

.clpmPhenoBeta <- function(vn, start) {
  # Phenotypic regression matrix: autoregression + distal cross-lags (free),
  # proximal cross-sectional DoC paths (labelled but fixed at 0).
  nv   <- 2L * vn$Ts
  vals <- matrix(0, nv, nv, dimnames = list(vn$labs, vn$labs))
  free <- matrix(FALSE, nv, nv, dimnames = list(vn$labs, vn$labs))
  labs <- matrix(NA_character_, nv, nv, dimnames = list(vn$labs, vn$labs))

  for (i in seq_len(vn$Ts - 1L)) {
    s <- vn$waveLabels[i]; t <- vn$waveLabels[i + 1L]
    xs <- vn$x[i];      ys <- vn$y[i]
    xt <- vn$x[i + 1L]; yt <- vn$y[i + 1L]
    vals[xt, xs] <- vals[yt, ys] <- start$arPheno
    vals[yt, xs] <- vals[xt, ys] <- start$clPheno
    free[xt, xs] <- free[yt, ys] <- free[yt, xs] <- free[xt, ys] <- TRUE
    labs[xt, xs] <- paste0("AR_bx", t, "x", s)
    labs[yt, ys] <- paste0("AR_by", t, "y", s)
    labs[yt, xs] <- paste0("CL_by", t, "x", s)
    labs[xt, ys] <- paste0("CL_bx", t, "y", s)
  }
  for (i in seq_along(vn$waveLabels)) {
    w <- vn$waveLabels[i]
    labs[vn$y[i], vn$x[i]] <- paste0("Inst_by", w, "x", w)
    labs[vn$x[i], vn$y[i]] <- paste0("Inst_bx", w, "y", w)
  }
  mxMatrix(type = "Full", nrow = nv, ncol = nv, free = free, values = vals,
           labels = labs, dimnames = list(vn$labs, vn$labs), name = "BE")
}

.clpmData <- function(d, what) {
  if (is(d, "MxData")) {
    if (!is.null(d$observedStats)) {
      return(mxData(type = "none", observedStats = d$observedStats,
                    numObs = d$numObs))
    }
    return(d)
  }
  if (is.data.frame(d) || is.matrix(d)) return(mxData(as.data.frame(d), type = "raw"))
  stop(what, " must be a data.frame or an MxData object (e.g. from ",
       "omxAugmentDataWithWLSSummary()).", call. = FALSE)
}


# ------------------------------------------------------------------------------
#' Build a twin cross-lagged panel model
#'
#' Returns the *standard twin CLPM* (Figure 1A of the paper): phenotypic
#' autoregression and phenotypic cross-lags, with an ACE decomposition of the
#' innovation (co)variances at each wave. This is the baseline from which the
#' Biometrical CLPM is reached by freeing latent autoregressions with
#' `useBiometricalAR()`; alternative sources of the cross-lagged association are
#' tested with `useLiabilityCrossLags()` and `addProximalDoC()`.
#'
#' The model is a multi-group (MZ, DZ) model. All matrices and the standardising
#' algebras are also placed at the top level, so `fit$V`, `fit$rA`, `fit$stdBeta`
#' etc. are available on the fitted object.
#'
#' @param xVar,yVar Base names of the two phenotypes; the wave label is pasted
#'   directly onto these (see `clpmVarNames`).
#' @param waveLabels Integer vector of wave labels, length Ts. Need not be
#'   consecutive: use `c(1, 6)` to label a two-wave subset of a six-wave design
#'   so that parameter labels match across analyses.
#' @param mzData,dzData Data for the two zygosity groups. Either raw
#'   `data.frame`s whose columns include `clpmVarNames(...)$selVars`, or `MxData`
#'   objects carrying `observedStats` (from `omxAugmentDataWithWLSSummary()`)
#'   for weighted least squares.
#' @param estimator `"ML"` (raw full-information maximum likelihood) or `"WLS"`
#'   (weighted least squares on summary statistics, as used in the paper).
#' @param components Variance components to estimate; a subset of `c("A","C","E")`.
#'   Components omitted here are fixed at zero.
#' @param ordinal Which phenotypes are ordinal: any of `c("x","y")`. Ordinal
#'   variables must be `mxFactor()`s in raw data.
#' @param nThresholds Number of thresholds per ordinal variable (levels - 1).
#' @param freeThresholds Estimate thresholds? Defaults to `TRUE` for ML and
#'   `FALSE` for WLS (where they are carried in `observedStats`).
#' @param equalMeansByTwin Constrain twin 1 and twin 2 to the same means/intercepts.
#' @param allContinuousMethod Passed to `mxFitFunctionWLS()`; the paper used
#'   `"marginals"`.
#' @param start Start values, see `clpmStart()`.
#' @param name Model name.
#'
#' @return An `MxModel`, not yet run. Fit with `mxRun()` or `mxTryHard()`.
#'
#' @examples
#' \dontrun{
#' vn  <- clpmVarNames("DepSxResRN_T", "CigDay3L_T", waveLabels = c(1, 6))
#' m0  <- biometricalCLPM("DepSxResRN_T", "CigDay3L_T", waveLabels = c(1, 6),
#'                        mzData = MZsumStat, dzData = DZsumStat,
#'                        estimator = "WLS", components = c("A", "E"),
#'                        ordinal = "y", nThresholds = 2)
#' fit0 <- mxRun(m0)                       # standard twin CLPM
#' fit1 <- mxRun(useBiometricalAR(fit0))   # Biometrical CLPM
#' mxCompare(fit1, fit0)
#' }
# ------------------------------------------------------------------------------
biometricalCLPM <- function(xVar = "X", yVar = "Y",
                            waveLabels = 1:2,
                            mzData, dzData,
                            estimator = c("ML", "WLS"),
                            components = c("A", "C", "E"),
                            ordinal = character(0),
                            nThresholds = 2L,
                            freeThresholds = NULL,
                            twinSuffix = c("_Tw1", "_Tw2"),
                            equalMeansByTwin = TRUE,
                            allContinuousMethod = "marginals",
                            start = clpmStart(),
                            name = "CLPM_phenoAR") {

  estimator  <- match.arg(estimator)
  components <- match.arg(components, c("A", "C", "E"), several.ok = TRUE)
  ordinal    <- if (length(ordinal)) {
    match.arg(ordinal, c("x", "y"), several.ok = TRUE)
  } else character(0)
  if (!"A" %in% components || !"E" %in% components) {
    stop("`components` must include at least \"A\" and \"E\".", call. = FALSE)
  }
  if (is.null(freeThresholds)) freeThresholds <- (estimator == "ML")

  vn  <- clpmVarNames(xVar, yVar, waveLabels, twinSuffix)
  nv  <- 2L * vn$Ts

  ## ---- means / intercepts ----
  B0 <- mxMatrix(type = "Full", nrow = 1, ncol = nv, free = TRUE,
                 values = start$mean, labels = paste0("b0_", vn$labs),
                 name = "b0", dimnames = list("b0", vn$labs))
  ExpMean <- if (equalMeansByTwin) {
    mxAlgebra(cbind(b0, b0), name = "expMean")
  } else {
    stop("equalMeansByTwin = FALSE is not implemented; supply your own ",
         "means matrix if twin-order differences are of interest.", call. = FALSE)
  }

  ## ---- thresholds for ordinal phenotypes ----
  ordNames <- character(0)
  if (length(ordinal)) {
    ordBase  <- unlist(vn[ordinal], use.names = FALSE)
    ordBase  <- vn$labs[vn$labs %in% ordBase]            # keep model ordering
    ordNames <- c(paste0(ordBase, twinSuffix[1]), paste0(ordBase, twinSuffix[2]))
    thLabs   <- paste0("th", seq_len(nThresholds), "_", rep(ordBase, each = nThresholds))
    threG <- mxMatrix(type = "Full", nrow = nThresholds, ncol = length(ordNames),
                      free = freeThresholds,
                      values = rep(start$thresholds[seq_len(nThresholds)],
                                   length(ordNames)),
                      labels = rep(thLabs, 2L),
                      dimnames = list(paste0("th", seq_len(nThresholds)), ordNames),
                      name = "threG")
  }

  ## ---- structural matrices ----
  BE  <- .clpmPhenoBeta(vn, start)
  I4  <- mxMatrix(type = "Iden", nrow = nv, ncol = nv, name = "I4",
                  dimnames = list(vn$labs, vn$labs))

  psi   <- lapply(c("A", "C", "E"), .clpmPsi, vn = vn, start = start)
  lbeta <- lapply(c("A", "C", "E"), .clpmLatentBeta, vn = vn)
  names(psi) <- names(lbeta) <- c("A", "C", "E")

  ## Fix omitted components at zero
  for (k in setdiff(c("A", "C", "E"), components)) {
    psi[[k]]@free[]   <- FALSE
    psi[[k]]@values[] <- 0
  }

  inv <- list(
    mxAlgebra(solve(I4 - BEA), name = "BAi"),
    mxAlgebra(solve(I4 - BEC), name = "BCi"),
    mxAlgebra(solve(I4 - BEE), name = "BEi")
  )
  vcomp <- list(
    mxAlgebra(BAi %*% PSA %*% t(BAi), name = "VA",
              dimnames = list(paste0("A", vn$labs), paste0("A", vn$labs))),
    mxAlgebra(BCi %*% PSC %*% t(BCi), name = "VC",
              dimnames = list(paste0("C", vn$labs), paste0("C", vn$labs))),
    mxAlgebra(BEi %*% PSE %*% t(BEi), name = "VE",
              dimnames = list(paste0("E", vn$labs), paste0("E", vn$labs)))
  )

  covP1 <- mxAlgebra(VA + VC + VE, name = "V1", dimnames = list(vn$labs, vn$labs))
  iBE   <- mxAlgebra(solve(I4 - BE), name = "iBE")
  covP  <- mxAlgebra(iBE %*% V1 %*% t(iBE), name = "V", dimnames = list(vn$labs, vn$labs))
  covMZ <- mxAlgebra(iBE %*% (VA + VC) %*% t(iBE), name = "cMZ")
  covDZ <- mxAlgebra(iBE %*% (0.5 %x% VA + VC) %*% t(iBE), name = "cDZ")
  expCovMZ <- mxAlgebra(rbind(cbind(V, cMZ), cbind(t(cMZ), V)), name = "expCovMZ")
  expCovDZ <- mxAlgebra(rbind(cbind(V, cDZ), cbind(t(cDZ), V)), name = "expCovDZ")

  ## ---- data, expectation, fit function ----
  dataMZ <- .clpmData(mzData, "mzData")
  dataDZ <- .clpmData(dzData, "dzData")

  expArgs <- list(means = "expMean", dimnames = vn$selVars)
  if (length(ordNames)) {
    expArgs$thresholds  <- "threG"
    expArgs$threshnames <- ordNames
  }
  expMZ <- do.call(mxExpectationNormal, c(list(covariance = "expCovMZ"), expArgs))
  expDZ <- do.call(mxExpectationNormal, c(list(covariance = "expCovDZ"), expArgs))

  fitFun <- if (estimator == "WLS") {
    mxFitFunctionWLS(type = "WLS", allContinuousMethod = allContinuousMethod)
  } else {
    mxFitFunctionML()
  }

  ## ---- standardised algebras ----
  calc <- .clpmStdAlgebras(vn, components)

  pars <- c(list(B0, BE, I4, iBE, covP1, covP),
            if (length(ordNames)) list(threG) else NULL,
            psi, lbeta, inv, vcomp)

  modelMZ <- mxModel(pars, ExpMean, covMZ, expCovMZ, dataMZ, expMZ, fitFun, name = "MZ")
  modelDZ <- mxModel(pars, ExpMean, covDZ, expCovDZ, dataDZ, expDZ, fitFun, name = "DZ")
  multi   <- mxFitFunctionMultigroup(c("MZ", "DZ"))

  mxModel(name = name, pars, calc, modelMZ, modelDZ, multi)
}


.clpmStdAlgebras <- function(vn, components) {
  aL <- paste0("A", vn$labs); cL <- paste0("C", vn$labs); eL <- paste0("E", vn$labs)
  out <- list(
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(VA))) %&% VA, name = "rA", dimnames = list(aL, aL)),
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(VE))) %&% VE, name = "rE", dimnames = list(eL, eL)),
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(PSA))) %&% PSA, name = "psiA", dimnames = list(aL, aL)),
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(PSE))) %&% PSE, name = "psiE", dimnames = list(eL, eL)),
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(V))) %*% BE %*% vec2diag(sqrt(diag2vec(V))),
              name = "stdBeta", dimnames = list(vn$labs, vn$labs)),
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(VA))) %*% BEA %*% vec2diag(sqrt(diag2vec(VA))),
              name = "stdBetaA", dimnames = list(aL, aL)),
    mxAlgebra(vec2diag(1 / sqrt(diag2vec(VE))) %*% BEE %*% vec2diag(sqrt(diag2vec(VE))),
              name = "stdBetaE", dimnames = list(eL, eL))
  )
  if ("C" %in% components) {
    out <- c(out, list(
      mxAlgebra(vec2diag(1 / sqrt(diag2vec(VC))) %&% VC, name = "rC", dimnames = list(cL, cL)),
      mxAlgebra(vec2diag(1 / sqrt(diag2vec(PSC))) %&% PSC, name = "psiC", dimnames = list(cL, cL)),
      mxAlgebra(vec2diag(1 / sqrt(diag2vec(VC))) %*% BEC %*% vec2diag(sqrt(diag2vec(VC))),
                name = "stdBetaC", dimnames = list(cL, cL))
    ))
  }
  nv    <- 2L * vn$Ts
  rowUS <- paste0("US_", vn$labs)
  colUS <- paste0(rep(c("VA_", "VC_", "VE_", "SA_", "SC_", "SE_"), each = nv), vn$labs)
  c(out, list(mxAlgebra(cbind(VA, VC, VE, VA / V, VC / V, VE / V),
                        name = "US", dimnames = list(rowUS, colUS))))
}


# ---- label extraction from a built model -------------------------------------

.clpmMatLabels <- function(model, matName, pattern) {
  m <- model[[matName]]
  if (is.null(m)) stop("Matrix '", matName, "' not found; is this a model built ",
                       "by biometricalCLPM()?", call. = FALSE)
  l <- as.vector(m$labels)
  unique(l[!is.na(l) & grepl(pattern, l)])
}


# ------------------------------------------------------------------------------
#' Replace phenotypic autoregression with latent (biometrical) autoregression
#'
#' Turns the standard twin CLPM into the **Biometrical CLPM** (Figure 1B): the
#' phenotypic autoregressive paths are fixed at zero and the autoregressions of
#' the latent A, C and E factors are freed, so that stability is modelled as
#' separate genetic and environmental developmental processes (a biometrical
#' simplex). The standard CLPM is nested within the result.
#'
#' @param model A model from `biometricalCLPM()` (fitted or not).
#' @param components Which latent autoregressions to free; defaults to whichever
#'   of A, C, E are estimated in `model`.
#' @param value Start value for the freed latent autoregressions.
#' @param keepPhenotypic Keep the phenotypic autoregressions free as well?
#'   Not identified in most designs; `FALSE` by default.
#' @param name Name for the new model.
#' @return A new `MxModel`.
# ------------------------------------------------------------------------------
useBiometricalAR <- function(model, components = NULL, value = 0.5,
                             keepPhenotypic = FALSE, name = "CLPM_aceAR") {
  components <- .clpmActiveComponents(model, components)
  out <- mxModel(model, name = name)
  if (!keepPhenotypic) {
    out <- omxSetParameters(out, labels = .clpmMatLabels(model, "BE", "^AR_b"),
                            free = FALSE, values = 0)
  }
  for (k in components) {
    out <- omxSetParameters(out,
                            labels = .clpmMatLabels(model, paste0("BE", k), paste0("^AR_", k)),
                            free = TRUE, values = value)
  }
  out
}


# ------------------------------------------------------------------------------
#' Free cross-lagged paths between latent A, C and/or E factors
#'
#' Tests cross-lagged *confounding* by correlated genetic and environmental
#' liabilities (Figure 1C) as an alternative to phenotypic causation. Set
#' `dropPhenotypic = TRUE` to fit the pure confounding model, in which the
#' phenotypic cross-lags are fixed at zero; leave it `FALSE` to fit the hybrid
#' model with both.
#'
#' @param model A model from `biometricalCLPM()` / `useBiometricalAR()`.
#' @param components Which latent cross-lags to free; a subset of `c("A","C","E")`.
#' @param value Start value.
#' @param dropPhenotypic Fix the phenotypic cross-lags at zero.
#' @param name Name for the new model.
#' @return A new `MxModel`.
# ------------------------------------------------------------------------------
useLiabilityCrossLags <- function(model, components = "A", value = 0.1,
                                  dropPhenotypic = FALSE, name = NULL) {
  components <- match.arg(components, c("A", "C", "E"), several.ok = TRUE)
  if (is.null(name)) {
    name <- paste0("CLPM_CL", paste(components, collapse = ""),
                   if (dropPhenotypic) "_noPheno" else "")
  }
  out <- mxModel(model, name = name)
  for (k in components) {
    out <- omxSetParameters(out,
                            labels = .clpmMatLabels(model, paste0("BE", k), paste0("^CL_", k)),
                            free = TRUE, values = value)
  }
  if (dropPhenotypic) out <- setPhenotypicCrossLags(out, free = FALSE, name = name)
  out
}


# ------------------------------------------------------------------------------
#' Free or fix the phenotypic (distal) cross-lagged paths
#'
#' @param model A model from `biometricalCLPM()`.
#' @param free Free (`TRUE`) or fix at zero (`FALSE`).
#' @param direction `"both"`, `"xy"` (y -> x only) or `"yx"` (x -> y only).
#' @param waves Restrict to cross-lags arriving at these wave labels; default all.
#' @param value Start value when freeing.
#' @param name Name for the new model.
#' @return A new `MxModel`.
# ------------------------------------------------------------------------------
setPhenotypicCrossLags <- function(model, free = FALSE,
                                   direction = c("both", "xy", "yx"),
                                   waves = NULL, value = 0.1, name = NULL) {
  direction <- match.arg(direction)
  labs <- switch(direction,
                 both = .clpmMatLabels(model, "BE", "^CL_b"),
                 yx   = .clpmMatLabels(model, "BE", "^CL_by"),
                 xy   = .clpmMatLabels(model, "BE", "^CL_bx"))
  if (!is.null(waves)) {
    keep <- paste0("^CL_b[xy](", paste(waves, collapse = "|"), ")[xy]")
    labs <- labs[grepl(keep, labs)]
  }
  if (is.null(name)) {
    name <- paste0(model$name, if (free) "_CL" else "_noCL",
                   if (direction == "both") "" else direction)
  }
  out <- mxModel(model, name = name)
  omxSetParameters(out, labels = labs, free = free,
                   values = if (free) value else 0)
}


# ------------------------------------------------------------------------------
#' Add proximal (cross-sectional) direction-of-causation effects
#'
#' Applies the Direction-of-Causation model within the CLPM (Figure 1D): a
#' cross-sectional causal path is freed at a given wave, and the corresponding
#' within-wave A and/or E covariance is fixed at zero so that the model remains
#' identified. Compare the resulting models by `mxCompare()` / pseudo-AIC to
#' choose among the causal directions.
#'
#' In a hybrid DoC model at a wave with A, C and E variance, at most three of the
#' five sources of within-wave covariance (two causal paths, covA, covC, covE)
#' can be estimated; with AE variance, at most two of four.
#'
#' @param model A model from `biometricalCLPM()` / `useBiometricalAR()`.
#' @param wave The wave *label* at which to estimate proximal effects.
#' @param direction `"x->y"`, `"y->x"` or `"both"`.
#' @param dropCov Which within-wave innovation covariances to fix at zero at
#'   `wave`; a subset of `c("A","C","E")`. Fix one per freed causal path.
#' @param value Start value for the freed proximal path(s).
#' @param name Name for the new model.
#' @return A new `MxModel`.
#'
#' @examples
#' \dontrun{
#' # Smoking -> depression at wave 6, identified by dropping covA at wave 6
#' fitDoC <- mxRun(addProximalDoC(fitNoLag, wave = 6,
#'                                direction = "y->x", dropCov = "A"))
#' }
# ------------------------------------------------------------------------------
addProximalDoC <- function(model, wave,
                           direction = c("both", "x->y", "y->x"),
                           dropCov = "A", value = 0.2, name = NULL) {
  direction <- match.arg(direction)
  dropCov   <- match.arg(dropCov, c("A", "C", "E"), several.ok = TRUE)
  stopifnot(length(wave) == 1L)

  instYX <- paste0("Inst_by", wave, "x", wave)   # x -> y
  instXY <- paste0("Inst_bx", wave, "y", wave)   # y -> x
  free   <- switch(direction, both = c(instYX, instXY),
                   `x->y` = instYX, `y->x` = instXY)

  known <- .clpmMatLabels(model, "BE", "^Inst_b")
  if (!all(free %in% known)) {
    stop("No proximal path labelled ", setdiff(free, known)[1],
         " in this model. Available waves: ",
         paste(unique(sub("^Inst_b[xy]([0-9]+)[xy].*$", "\\1", known)), collapse = ", "),
         call. = FALSE)
  }
  if (length(dropCov) < length(free)) {
    warning("Freeing ", length(free), " proximal path(s) but dropping only ",
            length(dropCov), " covariance(s); the model may not be identified. ",
            "Check with mxCheckIdentification().", call. = FALSE)
  }

  if (is.null(name)) {
    name <- paste0("CLPM_DoC", wave, "_",
                   switch(direction, both = "bidir", `x->y` = "yx", `y->x` = "xy"),
                   "_r", paste(dropCov, collapse = ""))
  }
  out <- mxModel(model, name = name)
  out <- omxSetParameters(out, labels = paste0("cov", dropCov, "xy", wave),
                          free = FALSE, values = 0)
  omxSetParameters(out, labels = free, free = TRUE, values = value)
}


.clpmActiveComponents <- function(model, components = NULL) {
  if (!is.null(components)) {
    return(match.arg(components, c("A", "C", "E"), several.ok = TRUE))
  }
  active <- character(0)
  for (k in c("A", "C", "E")) {
    m <- model[[paste0("PS", k)]]
    if (!is.null(m) && any(m$free)) active <- c(active, k)
  }
  active
}


# ------------------------------------------------------------------------------
#' Tidy the standardised estimates of a fitted Biometrical CLPM
#'
#' @param fit A fitted model from `mxRun()`.
#' @param what Which algebra to extract: `"stdBeta"` (phenotypic), `"stdBetaA"`,
#'   `"stdBetaE"`, `"stdBetaC"`, `"psiA"`, `"psiE"`, `"psiC"`, `"rA"`, `"rE"`, `"rC"`.
#' @param se Compute standard errors with `mxSE()` (slower).
#' @param digits Rounding.
#' @return A `data.frame` of non-zero elements with `from`, `to`, `estimate`
#'   and (optionally) `se`.
# ------------------------------------------------------------------------------
clpmStdEstimates <- function(fit, what = "stdBeta", se = FALSE, digits = 3) {
  alg <- fit[[what]]
  if (is.null(alg)) stop("Algebra '", what, "' not found in the fitted model.",
                         call. = FALSE)
  est <- alg$result
  idx <- which(est != 0, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(from = character(0), to = character(0),
                                    estimate = numeric(0)))
  out <- data.frame(
    from     = colnames(est)[idx[, "col"]],
    to       = rownames(est)[idx[, "row"]],
    estimate = round(est[idx], digits),
    stringsAsFactors = FALSE
  )
  if (se) {
    s <- mxSE(what, model = fit, silent = TRUE)
    out$se <- round(s[idx], digits)
  }
  out[order(out$to, out$from), ]
}


# ------------------------------------------------------------------------------
#' Build WLS summary statistics for one zygosity group
#'
#' Convenience wrapper around `omxAugmentDataWithWLSSummary()` reproducing the
#' preparation in `TEDS/3_Get_DepSx_CigDay_SumStats_for_WLS.R`. The raw data are
#' stripped from the returned object, so the summary statistics can be shared
#' without sharing individual-level data.
#'
#' @param data A `data.frame` in family-wise (wide) format for one zygosity group.
#' @param selVars Variables to retain, in model order (`clpmVarNames(...)$selVars`).
#' @param ordVars Which of `selVars` are ordinal; converted with `mxFactor()`.
#' @param nThresholds Number of thresholds for the ordinal variables.
#' @param allContinuousMethod Passed through; the paper used `"marginals"`.
#' @param dropRaw Remove the individual-level data from the returned object.
#' @return An `MxData` object with `observedStats`, ready for `biometricalCLPM()`.
# ------------------------------------------------------------------------------
clpmSumStats <- function(data, selVars, ordVars = character(0), nThresholds = 2L,
                         allContinuousMethod = "marginals", dropRaw = TRUE) {
  d <- as.data.frame(data)[, selVars, drop = FALSE]
  d <- d[rowSums(is.na(d)) < ncol(d), , drop = FALSE]     # drop empty rows
  if (length(ordVars)) {
    d[ordVars] <- mxFactor(d[ordVars], levels = 0:nThresholds)
  }
  out <- omxAugmentDataWithWLSSummary(mxData(d, type = "raw"),
                                      allContinuousMethod = allContinuousMethod,
                                      type = "WLS")
  if (dropRaw) out@observed <- NULL
  out
}

# END ==========================================================================

#' Generate Synthetic Multi-Modal Data for HMM-Based Supervised Biclustering
#'
#' Simulates multi-modality feature matrices, covariates, and an outcome under a
#' cluster/state structure used by \code{\link{iBiclust}}. For each modality,
#' features are assigned to one of three states within each cluster:
#' over-expressed (\code{1}), under-expressed (\code{2}), or baseline (\code{3}).
#' Outcomes can be continuous, binary, or censored.
#'
#' @details
#' \strong{Cluster assignment:} Samples are assigned to \code{K} clusters in
#' approximately equal-sized blocks.
#'
#' \strong{Feature membership (\code{kappajk}):} For each modality \eqn{m}, a
#' \eqn{p_m \times K} matrix is generated with entries in \code{\{1,2,3\}}.
#' If \code{OrderFeatured = TRUE}, features 1--100 are set to state \code{1}
#' (over-expressed) and features 101--200 to state \code{2} (under-expressed).
#' Otherwise, 100 over- and 100 under-expressed features are sampled at random
#' for each modality (requires \code{p[m] >= 200}).
#'
#' \strong{Outcome generation:} A latent mean \code{MeanY[Zk[i]]} is assigned to
#' each sample and combined with a linear covariate effect. Continuous outcomes
#' add Gaussian noise; binary outcomes threshold a noisy latent score; censored
#' outcomes return an \eqn{n \times 2} matrix \code{cbind(time, censor)},
#' where the second column is a censoring indicator (1 = censored, 0 = event observed).
#'
#' @param OrderFeatured Logical; if \code{TRUE}, uses a deterministic ordering for
#'   feature membership in each modality (features 1--100 over-expressed and
#'   101--200 under-expressed). If \code{FALSE}, feature membership is randomized.
#' @param K Integer; number of clusters.
#' @param NbrPlat Integer; number of modalities/platforms/views.
#' @param OutcType Outcome type: one of \code{"continuous"}, \code{"binary"}, or \code{"censored"}.
#' @param MeanY Numeric vector of length \code{K}; cluster-specific mean outcome levels.
#' @param pc Integer; number of covariates to generate. Must be a positive integer.
#' @param betacov Numeric; magnitude of the covariate coefficients. Coefficients
#'   are assigned as \code{+(betacov)} for the first half and \code{-(betacov)}
#'   for the remainder.
#' @param MeanOverExpre Numeric vector of length \code{K}; means for over-expressed
#'   features in each cluster.
#' @param MeanUnderExpre Numeric vector of length \code{K}; means for under-expressed
#'   features in each cluster.
#' @param p Integer vector of length \code{NbrPlat}; number of features per modality.
#' @param seed Integer; random seed used for simulating \code{X} and outcomes.
#' @param n Integer; number of samples.
#' @param SigmaOE Numeric vector of length \code{K}; standard deviations for over-expressed features.
#' @param SigmaUE Numeric vector of length \code{K}; standard deviations for under-expressed features.
#'
#' @returns A list with components:
#' \describe{
#'   \item{Y}{Outcome vector of length \code{n} for continuous/binary outcomes, or
#'     an \eqn{n \times 2} matrix for censored outcomes \code{cbind(time, censor)}, where the second column is the
#'     censoring indicator (1 = censored, 0 = event observed).}
#'   \item{X}{List of length \code{NbrPlat}; each element is an \eqn{n \times p_m}
#'     matrix of simulated features.}
#'   \item{kappajk}{List of length \code{NbrPlat}; each element is a \eqn{p_m \times K}
#'     matrix with entries in \code{\{1,2,3\}} indicating over/under/baseline feature
#'     membership for each cluster.}
#'   \item{Zk}{Integer vector of length \code{n} giving true cluster assignments in \code{1:K}.}
#'   \item{Cov}{Numeric matrix of dimension \eqn{n \times pc} containing generated covariates.}
#' }
#'
#' @examples
#' set.seed(1)
#' K <- 2
#' p <- c(300, 500)
#' n <- 100
#' dat <- iGenDataHMM(
#'   OrderFeatured = TRUE, K = K, NbrPlat = 2, OutcType = "continuous",
#'   MeanY = c(1, -1), pc = 4, betacov = 1.0,
#'   MeanOverExpre = (1:K) / 2, MeanUnderExpre = -(1:K) / 2,
#'   p = p, seed = 1, n = n,
#'   SigmaOE = rep(1, K), SigmaUE = rep(1, K)
#' )
#' str(dat)
#'
#' @export

iGenDataHMM <- function(
    OrderFeatured = FALSE, K = K, NbrPlat = 2, OutcType = "continuous",
    MeanY = MeanY, pc = NULL, betacov = 1, MeanOverExpre = 1:K,
    MeanUnderExpre = -c(1:K), p = p, seed = 1, n = n,
    SigmaOE = rep(1, K), SigmaUE = rep(1, K)) {
  Zk <- rep(1, n)
  for (i in 1:n) {
    for (k in 1:K) {
      # assign 1/K percent samples in cluster K
      if ((1 + n * (k - 1) / K <= i) && (i <= n * k / K)) {
        Zk[i] <- k
      }
    }
  }
  ### Create kappa, feature membership
  KappaKnown <- OrderFeatured
  set.seed(1)
  kappajk <- list()
  for (m in 1:NbrPlat) {
    # create a p[m] x K matrix filled with 3
    kappajk[[m]] <- matrix(3, p[m], K)
    if (KappaKnown == TRUE) {
      kappajk[[m]][1:100, ] <- 1
      kappajk[[m]][101:200, ] <- 2
    } else {
      # random sample under and over expressed genes in each platform
      x1 <- sample(p[m], 100)
      kappajk[[m]][x1, ] <- 1
      x2 <- sample(c(1:p[m])[-x1], 100)
      kappajk[[m]][x2, ] <- 2
    }
  }

  X <- list()
  set.seed(seed)
  for (m in 1:NbrPlat) {
    X[[m]] <- matrix(0, n, p[m])
    for (j in 1:p[m]) { # for each genes
      for (i in 1:n) {
        for (k in 1:K) {
          if (Zk[i] == k) {
            if ((kappajk[[m]][j, k] == 1)) {
              X[[m]][i, j] <- rnorm(1, mean = MeanOverExpre[k], sd = SigmaOE[k])
            } else if ((kappajk[[m]][j, k] == 2)) {
              X[[m]][i, j] <- rnorm(1, mean = MeanUnderExpre[k], sd = SigmaUE[k])
            } else {
              X[[m]][i, j] <- rnorm(1, mean = 0, sd = 1)
            }
          }
        }
      }
      # K is the number of clusters
    }
  }
  Xcov <- scale(matrix(rnorm(pc * n), nrow = n))
  Beta <- c(rep(betacov, floor(pc / 2)), rep(-betacov, pc - floor(pc / 2)))
  # print(Beta)

  ### generate the response
  Y <- rep(0, n)
  for (i in 1:n) {
    for (k in 1:K) {
      if (Zk[i] == k) {
        Y[i] <- MeanY[k]
      }
    }
  }
  Y <- Y + Xcov %*% Beta
  if (OutcType == "continuous") {
    Youtcome <- Y + rnorm(n)
  } else if (OutcType == "binary") {
    Youtcome <- 1 * (Y + rnorm(n) > 0)
  } else if (OutcType == "censored") {
    U <- Y + runif(n, min = 0, max = 1) ### 30% censored data
    Y <- Y + rnorm(n)
    Youtcome <- as.matrix(cbind(pmin(U, Y), 1 * (U < Y)))
  }
  return(list(Y = Youtcome, X = X, kappajk = kappajk, Zk = Zk, Cov = Xcov))
}

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
#' If \code{OrderFeatured = TRUE}, the first nActive features are set to state \code{1}
#' (over-expressed) and the next nActive features to state \code{2} (under-expressed).
#' Otherwise, nActive over- and nActive under-expressed features are sampled at random
#' for each modality (requires \code{p[m] >= 2 * nActive[m]}).
#'
#' \strong{Outcome generation:} A latent mean \code{MeanY[Zk[i]]} is assigned to
#' each sample and combined with a linear covariate effect. Continuous outcomes
#' add Gaussian noise; binary outcomes threshold a noisy latent score; censored
#' outcomes return an \eqn{n \times 2} matrix \code{cbind(time, censor)},
#' where the second column is a censoring indicator (1 = censored, 0 = event observed).
#'
#' @param OrderFeatured Logical; if \code{TRUE}, uses a deterministic ordering for
#'   feature membership in each modality (the first nActive features over-expressed and
#'   the next nActive under-expressed). If \code{FALSE}, feature membership is randomized.
#' @param nActive Positive integer or length-NbrPlat vector; number of features
#'   in each active state per modality (default 100).
#' @param K Integer; number of clusters.
#' @param NbrPlat Integer; number of modalities/platforms/views.
#' @param OutcType Outcome type: one of \code{"continuous"}, \code{"binary"}, or \code{"censored"}.
#' @param MeanY Numeric scalar, length-\code{K} vector, or NbrPlat by K matrix; cluster-specific mean outcome levels.
#' @param pc Integer; number of covariates to generate. Must be a positive integer.
#' @param betacov Numeric; magnitude of the covariate coefficients. Coefficients
#'   are assigned as \code{+(betacov)} for the first half and \code{-(betacov)}
#'   for the remainder.
#' @param MeanOverExpre Numeric scalar, length-\code{K} vector, or NbrPlat by K matrix; means for over-expressed
#'   features in each cluster.
#' @param MeanUnderExpre Numeric scalar, length-\code{K} vector, or NbrPlat by K matrix; means for under-expressed
#'   features in each cluster.
#' @param p Integer vector of length \code{NbrPlat}; number of features per modality.
#' @param seed Integer; random seed used for simulating \code{X} and outcomes.
#' @param n Integer; number of samples.
#' @param SigmaOE Numeric scalar, length-\code{K} vector, or NbrPlat by K matrix; standard deviations for over-expressed features.
#' @param SigmaUE Numeric scalar, length-\code{K} vector, or NbrPlat by K matrix; standard deviations for under-expressed features.
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
#'   OrderFeatured = TRUE, nActive = 100, K = K, NbrPlat = 2,
#'   OutcType = "continuous", MeanY = c(1, -1), pc = 4, betacov = 1.0,
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
    SigmaOE = rep(1, K), SigmaUE = rep(1, K), nActive = 100L) {
  # Expression means and SDs accept either a shared length-K vector or an
  # NbrPlat x K matrix: rows are modalities, columns are sample clusters.
  # Scalars are also accepted. Existing vector inputs keep their behavior.
  as_modality_cluster <- function(x, name, positive = FALSE) {
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop(name, " must contain finite numeric values")
    }
    if (is.matrix(x)) {
      if (!identical(dim(x), as.integer(c(NbrPlat, K)))) {
        stop(name, " must have NbrPlat rows and K columns")
      }
    } else {
      if (!is.null(dim(x)) || !length(x) %in% c(1L, K)) {
        stop(name, " must be a scalar, length-K vector, or NbrPlat x K matrix")
      }
      x <- matrix(rep(rep(x, length.out = K), each = NbrPlat), NbrPlat, K)
    }
    if (positive && any(x <= 0)) stop(name, " must be positive (standard deviations)")
    x
  }
  MeanOverExpre <- as_modality_cluster(MeanOverExpre, "MeanOverExpre")
  MeanUnderExpre <- as_modality_cluster(MeanUnderExpre, "MeanUnderExpre")
  SigmaOE <- as_modality_cluster(SigmaOE, "SigmaOE", positive = TRUE)
  SigmaUE <- as_modality_cluster(SigmaUE, "SigmaUE", positive = TRUE)
  # Number of features in EACH active state, shared across clusters within a
  # modality. The default retains the original 100 over / 100 under features.
  if (length(p) != NbrPlat || any(!is.finite(p)) || any(p < 1 | p != floor(p))) {
    stop("p must contain one positive integer per modality")
  }
  if (!length(nActive) %in% c(1L, NbrPlat) || any(!is.finite(nActive)) ||
      any(nActive < 1 | nActive != floor(nActive))) {
    stop("nActive must be a positive integer or a length-NbrPlat vector")
  }
  nActive <- rep(nActive, length.out = NbrPlat)
  if (any(2 * nActive > p)) stop("Each p must be at least twice nActive")
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
      kappajk[[m]][seq_len(nActive[m]), ] <- 1
      kappajk[[m]][nActive[m] + seq_len(nActive[m]), ] <- 2
    } else {
      # Sample disjoint over- and underexpressed features in each modality.
      x1 <- sample(p[m], nActive[m])
      kappajk[[m]][x1, ] <- 1
      remaining <- seq_len(p[m])[-x1]
      x2 <- remaining[sample.int(length(remaining), nActive[m])]
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
              X[[m]][i, j] <- rnorm(1, mean = MeanOverExpre[m, k], sd = SigmaOE[m, k])
            } else if ((kappajk[[m]][j, k] == 2)) {
              X[[m]][i, j] <- rnorm(1, mean = MeanUnderExpre[m, k], sd = SigmaUE[m, k])
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
  list(Y = Youtcome, X = X, kappajk = kappajk, Zk = Zk, Cov = Xcov)
}

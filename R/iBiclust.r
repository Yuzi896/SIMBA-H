#' Bayesian Supervised Biclustering with SIMBA-H or SIMBA for Multi-Modal Data
#'
#' Fits a Bayesian model-based supervised biclustering model for multi-modality
#' datasets using a Gibbs sampling algorithm. SIMBA-H incorporates a Hidden
#' Markov Model (HMM) structure across ordered features; SIMBA does not incorporate
#' the HMM structure. Each modality/platform is provided as a matrix
#' with samples in rows and features in columns; internally, each matrix is
#' optionally standardized (centered/scaled) and transposed before being passed
#' to the native routine.
#'
#' @details
#' \strong{Data format:}
#' \itemize{
#'   \item \code{dataList} must be a list of length \code{M} (number of modalities/platforms).
#'   \item Each element \code{dataList[[t]]} is an \eqn{n \times p_t} matrix:
#'         \eqn{n} samples (rows) and \eqn{p_t} features (columns).
#'   \item If \code{normalized = TRUE}, each modality in \code{dataList} is centered and
#'         scaled feature-wise using its own mean and standard deviation. If
#'         \code{dataListTest} is provided, the test data are standardized using the
#'         corresponding training means and standard deviations (not re-estimated
#'         from the test data).
#' }
#'
#' \strong{Ordering requirement (HMM):}
#' If \code{Method = "SIMBA-H"} (the default), the features within each
#' modality in \code{dataList} are assumed to be pre-ordered according to
#' their natural or genomic sequence. The function does not reorder features
#' internally. For \code{Method = "SIMBA"}, no feature ordering is required.
#'
#'
#' \strong{Outcome types:}
#' \itemize{
#'   \item \code{"continuous"}: \code{OutcomeType = 1}
#'   \item \code{"binary"}: \code{OutcomeType = 2}
#'   \item \code{"censored"}: \code{OutcomeType = 3}, where \code{Youtcome} must be an
#'         \eqn{n \times 2} matrix (e.g., value and censoring indicator).
#' }
#'
#' \strong{Test data:}
#' If \code{dataListTest} is provided, it must have the same list structure as
#' \code{dataList}. If \code{normalized = TRUE}, test data are standardized using
#' the training means/SDs.
#'
#'
#' @param Method Character string specifying the modeling approach: either
#'   \code{"SIMBA-H"} (default), which incorporates a Hidden Markov Model (HMM)
#'   structure across ordered features, or \code{"SIMBA"}, which does not
#'   incorporate the HMM structure.
#' @param dataList A list of modality matrices. Each \code{dataList[[t]]} is an
#'   \eqn{n \times p_t} matrix with samples in rows and features in columns.
#' @param OutcType Outcome type: one of \code{"continuous"}, \code{"binary"}, or \code{"censored"}.
#' @param Youtcome Outcome vector/matrix. For censored outcomes,
#'   this must be an \eqn{n \times 2} matrix \code{cbind(time, censor)},
#'   where the second column represents the \emph{censoring indicator}
#'   (1 = censored, 0 = event observed).
#' @param Cov Covariate matrix of dimension \eqn{n \times p_c}.
#' @param sample Number of MCMC samples to retain.
#' @param burnin Number of burn-in iterations.
#' @param K Number of sample clusters used in the model.
#' @param normalized Logical; if \code{TRUE}, each modality is centered and scaled feature-wise.
#' @param seed Integer random seed passed to the native routine.
#' @param hyperSigma Numeric vector of length 2. The prior distribution on the
#'  state-specific variances \eqn{\sigma^2_{d k 1}} and \eqn{\sigma^2_{d k 2}}
#'  for \eqn{\mu_{d,j,k,l}}.
#' @param alphaP Numeric vector of length \code{K} specifying the Dirichlet prior
#'  concentration parameters for the cluster probability vector. These parameters
#'  control the prior probability that a sample belongs to each cluster \eqn{k}.
#' @param delt Numeric vector of length 3 giving the Dirichlet prior concentration
#'  parameters \eqn{(\delta_1, \delta_2, \delta_3)}. For
#'  \code{Method = "SIMBA-H"}, these specify the prior for each row of the
#'  transition matrix \eqn{\xi}. For \code{Method = "SIMBA"}, they specify
#'  the prior for the feature-state probability vector.
#' @param TruncVal Positive numeric truncation threshold \eqn{t} used in the Gibbs sampler.
#' @param dataListTest Optional test set list with same structure as \code{dataList}.
#' @param CovTest Optional test covariate matrix of dimension \eqn{n_{test} \times p_c}.
#'
#' @param return_samples Logical; additionally return ZSample, probZSample,
#'   MuSample, and probKappaSample when TRUE (default FALSE).
#' @returns A list with posterior summaries, coefficient draws, test labels and
#'   predictions, and log-posterior values. Only ZSample, probZSample, MuSample,
#'   and probKappaSample require \code{return_samples = TRUE}. Label switching
#'   is checked automatically. If label switching occurs, cluster labels are
#'   aligned using iterative ECR:
#' \describe{
#'   \item{ZSample, probZSample, MuSample, probKappaSample}{Retained labels, cluster
#'     probabilities, latent means, and feature-state probabilities when requested.}
#'   \item{DIC}{Deviance Information Criterion returned by the sampler.}
#'   \item{probZMean}{Matrix of posterior mean cluster membership probabilities (\eqn{n \times K}).}
#'   \item{MuMean}{Array of posterior means of the latent mean parameters \eqn{\mu_{dkl}} for the over- and under-expressed
#'   states (\eqn{l = 1,2}), indexed by modality \eqn{d} and cluster \eqn{k}.}
#'   \item{CoefySample}{cluster-specific outcome means with dimension \eqn{sample \times K}.}
#'   \item{sig2yM}{Posterior mean of outcome variance parameter.}
#'   \item{KappaMean}{List of arrays of posterior means for biclustering indicators per modality.}
#'   \item{TransMatMean}{Array of posterior mean transition matrices per modality
#'     for \code{Method = "SIMBA-H"}. For \code{Method = "SIMBA"}, each
#'     matrix has identical rows containing posterior mean feature-state probabilities.}
#'   \item{Betasample}{Matrix of sampled covariate coefficients (\eqn{sample \times p_c}).}
#'   \item{ZtestSample}{Matrix of sampled test cluster labels (\eqn{sample \times n_{test}}) if test data provided.}
#'   \item{YpredTestSample}{Matrix of sampled test predictions (\eqn{sample \times n_{test}}) if test data provided.}
#'   \item{logpostsample}{Vector of log-posterior values of length \code{sample + burnin}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Generate synthetic multi-modality data
#' K <- 2
#' p <- c(300, 500)   # smaller for an example
#' n <- 100
#' nbrcov <- 4
#' OutcType <- "continuous"
#'
#' dat <- iGenDataHMM(
#'   OrderFeatured = TRUE, nActive = 100, K = K, NbrPlat = 2, OutcType = OutcType,
#'   MeanY = c(rep(1, floor(K/2)), rep(-1, K - floor(K/2))),
#'   pc = nbrcov, betacov = 1.0,
#'   MeanOverExpre = (1:K) / 2, MeanUnderExpre = -(1:K) / 2,
#'   p = p, seed = 1, n = n,
#'   SigmaOE = rep(1, K), SigmaUE = rep(1, K)
#' )
#'
#' ## Generate test data
#' ntest <- 100
#' datTest <- iGenDataHMM(
#'   OrderFeatured = TRUE, nActive = 100, K = K, NbrPlat = 2, OutcType = OutcType,
#'   MeanY = c(rep(1, floor(K/2)), rep(-1, K - floor(K/2))),
#'   pc = nbrcov, betacov = 1.0,
#'   MeanOverExpre = (1:K) / 2, MeanUnderExpre = -(1:K) / 2,
#'   p = p, seed = 2, n = ntest,
#'   SigmaOE = rep(1, K), SigmaUE = rep(1, K)
#' )
#'
#' ## Fit SIMBA-H with the HMM structure (keep MCMC small in examples)
#' ## Use Method = "SIMBA" to fit without the HMM structure.
#' result <- iBiclust(
#'   Method = "SIMBA-H",
#'   dataList = dat$X,
#'   Youtcome = dat$Y,
#'   OutcType = OutcType,
#'   Cov = dat$Cov,
#'   sample = 200,
#'   burnin = 100,
#'   K = K,
#'   normalized = FALSE,
#'   seed = 1,
#'   hyperSigma = c(1, 1),
#'   alphaP = rep(1, K),
#'   delt = rep(1, 3),
#'   TruncVal = 0.1,
#'   dataListTest = datTest$X,
#'   CovTest = datTest$Cov
#' )
#'
#' str(result)
#' }

#'
#' @export

iBiclust <- function(
    Method = "SIMBA-H", dataList = dataList, OutcType = "continuous",
    Youtcome = Youtcome, Cov = NULL, sample = 2000, burnin = 1000, K = 2,
    normalized = TRUE, seed = 1, hyperSigma = c(1, 1), alphaP = rep(1, K),
    delt = rep(1, 3), TruncVal = .1, dataListTest = NULL, CovTest = NULL,
    return_samples = FALSE) {
  if (!is.character(Method) || length(Method) != 1L || is.na(Method) ||
      !Method %in% c("SIMBA-H", "SIMBA")) {
    stop('Method must be "SIMBA-H" or "SIMBA"')
  }

  # The native sampler uses 1 for SIMBA-H (HMM) and 2 for SIMBA (no HMM).
  MethodCode <- if (Method == "SIMBA-H") 1L else 2L

  # Latent draws are used for label switching; return those only on request.
  if (!is.logical(return_samples) || length(return_samples) != 1L || is.na(return_samples)) {
    stop("return_samples must be TRUE or FALSE")
  }
  if (OutcType == "continuous") {
    OutcomeType <- 1
  } else if (OutcType == "binary") {
    OutcomeType <- 2
  } else if (OutcType == "censored") {
    nc <- ncol(Youtcome)
    if (nc != 2) {
      print("Enter the outcome as an N times 2 matrix")
      stop("Enter the right outcome for censored outcome data")
    }
    OutcomeType <- 3
  }
  NbrPlat <- length(dataList) # number of platform but what is platform mean???
  # seems like platform only differ by the number of features

  P <- NULL # number of pathway; empty will be filled out later
  n <- nrow(dataList[[1]]) # number of samples
  Mean <- list()
  Sd <- list()
  for (t in 1:NbrPlat) { # for each platform
    dataList[[t]] <- as.matrix(dataList[[t]])
    P[t] <- ncol(dataList[[t]])
    Mean[[t]] <- apply(dataList[[t]], 2, mean)
    Sd[[t]] <- apply(dataList[[t]], 2, sd)
    if (normalized == TRUE) {
      dataList[[t]] <- scale(as.matrix(dataList[[t]]), center = TRUE, scale = TRUE)
    }
    dataList[[t]] <- t(dataList[[t]])
  }
  ntest <- 0
  if (!is.null(dataListTest)) {
    ntest <- nrow(dataListTest[[1]])
    # Mean <- list()
    # Var <- list()
    for (t in 1:NbrPlat) {
      dataListTest[[t]] <- as.matrix(dataListTest[[t]])
      if (normalized == TRUE) {
        for (j in 1:P[t]) {
          dataListTest[[t]][, j] <- (dataListTest[[t]][, j] - Mean[[t]][j]) / Sd[[t]][j]
        }
      }
      # Mean[[t]] <- apply(dataList[[t]], 2, mean)
      # Var[[t]] <- apply(dataList[[t]], 2, mean)
      dataListTest[[t]] <- as.matrix(t(dataListTest[[t]]))
    }
  }
  pc <- 0
  Cov1 <- NULL
  CovTest1 <- NULL
  if (!is.null(Cov)) {
    Cov <- as.matrix(Cov, nrow = n)
    pc <- ncol(Cov)
    Cov1 <- as.vector(t(Cov))
    if (!is.null(CovTest)) {
      CovTest <- as.matrix(CovTest, nrow = ntest)
      CovTest1 <- as.vector(t(CovTest))
    }
  }
  # Mean1 <- unlist(Mean)
  # Var1 <- unlist(Var)


  datasets <- unlist(dataList)
  datasetsTest <- unlist(dataListTest)
  truncatedvalue <- TruncVal
  sigma20Mu1 <- 100
  # void MCMCfunction(int *Method, int *M1, int *p1, int *nn, double *Data1, double *Y, int *KK, double *truncatedvalue,
  # double *rhosigma, double *gamsigma, double *alphaP, double *delt, double *sigma20Mu1, int *Nsample1,
  # int *burnin1, int *seed1,double *probZMean1, double *probkappajk1, double *TransMatMean1, double *MuMean1, d
  # ouble *CoefYMean, double *sig2yMean, float *logpostsample, float *DIC)
  sampler <- getNativeSymbolInfo("MCMCfunctionWithProbAndZSamples", PACKAGE = "SIMBAH")
  res <- .C(sampler,
            Method = as.integer(MethodCode), OutcomeType = as.integer(OutcomeType),
            M1 = as.integer(NbrPlat),
            p1 = as.integer(P), nn = as.integer(n), Data1 = as.double(datasets),
            Cov1 = as.double(Cov1), pc = as.integer(pc),
            Youtcome = as.double(Youtcome), KK = as.integer(K),
            truncatedvalue = as.double(truncatedvalue),
            rhosigma = as.double(hyperSigma[1]), gamsigma = as.double(hyperSigma[2]),
            alphaP = as.double(alphaP), delt = as.double(delt),
            sigma20Mu1 = as.double(sigma20Mu1), Nsample1 = as.integer(sample),
            burnin1 = as.integer(burnin), seed1 = as.integer(seed),
            probZMean1 = as.double(rep(0, n * K)),
            probkappajk1 = as.double(rep(0, sum(P) * K * 3)),
            TransMatMean1 = as.double(rep(0, NbrPlat * 3 * 3)),
            MuMean1 = as.double(rep(0, K * NbrPlat * 2)),
            CoefYsample = as.double(rep(0, K * sample)), sig2yMean = as.double(0),
            logpostsample = as.double(rep(0, sample + burnin)),
            Betasample = as.double(rep(0, sample * pc)),
            DIC = as.double(0), nnTest = as.integer(ntest),
            DataTest1 = as.double(datasetsTest), CovTest1 = as.double(CovTest1),
            YpredTest = as.double(rep(0, sample * ntest)),
            Ztest = as.integer(rep(0, sample * ntest)),
            MuSample = double(sample * NbrPlat * K * 2),
            probZSample = double(sample * n * K),
            probKappaSample = double(sample * sum(P) * K * 3),
            ZSample = integer(sample * n)
  )



  Betasample1 <- matrix(res$Betasample, nrow = sample, ncol = pc, byrow = TRUE)

  # SampleCoefy=matrix(res$SampleCoefy,sample,K,byrow = T)

  KappaMean <- list()
  probKappaSample <- vector("list", NbrPlat)
  kappaOffset <- 0
  MuMean <- aperm(array(res$MuMean1, c(2, K, NbrPlat))) # is the dimension flipped?
  TransMatMean <- aperm(array(res$TransMatMean1, c(3, 3, NbrPlat)))

  m1 <- 1
  for (t in 1:NbrPlat)
  {
    probk <- res$probkappajk1[m1:(m1 - 1 + P[t] * K * 3)]
    KappaMean[[t]] <- aperm(array(probk, c(3, P[t], K)))
    m1 <- m1 + P[t] * K * 3
    # Conditional probabilities: [iteration, feature, cluster, state].
    count <- sample * P[t] * K * 3
    draws <- res$probKappaSample[kappaOffset + seq_len(count)]
    probKappaSample[[t]] <- aperm(array(draws, c(3, K, P[t], sample)), c(4, 3, 2, 1))
    kappaOffset <- kappaOffset + count
  }
  probZMean <- matrix(res$probZMean1, n, byrow = TRUE)

  # Every post-burn-in iteration is retained in these arrays.
  # MuSample: [iteration, platform, cluster, expression state].
  MuSample <- aperm(array(res$MuSample, c(2, K, NbrPlat, sample)), c(4, 3, 2, 1))
  # probZSample: [iteration, individual, cluster].
  probZSample <- aperm(array(res$probZSample, c(K, n, sample)), c(3, 2, 1))
  # Actual sampled labels: [iteration, individual], aligned with probZSample.
  ZSample <- matrix(res$ZSample, nrow = sample, ncol = n, byrow = TRUE)

  # Keep the iteration dimension even when there are no test individuals.
  ZtestSample <- matrix(res$Ztest, nrow = sample, ncol = ntest, byrow = TRUE)
  CoefySample <- matrix(res$CoefYsample, nrow = sample, ncol = K, byrow = TRUE)
  Betasample <- Betasample1
  YpredTestSample <- matrix(res$YpredTest, nrow = sample, ncol = ntest, byrow = TRUE)

  # Each row lists the old cluster labels in their new order.
  run <- label.switching::ecr.iterative.2(z = ZSample, K = K, p = probZSample)
  permutations <- run$permutations
  identity_permutations <- matrix(seq_len(K), nrow = sample, ncol = K, byrow = TRUE)
  if (all(permutations == identity_permutations)) {
    print("No Label switching is detected")
  } else {
    print("Label switching detected")
    for (s in seq_len(sample)) {
      perm <- permutations[s, ]
      ZSample[s, ] <- match(ZSample[s, ], perm)
      probZSample[s, , ] <- probZSample[s, , perm]
      MuSample[s, , , ] <- MuSample[s, , perm, ]
      CoefySample[s, ] <- CoefySample[s, perm]
      for (d in seq_len(NbrPlat)) {
        probKappaSample[[d]][s, , , ] <- probKappaSample[[d]][s, , perm, ]
      }
      if (ntest > 0) {
        ZtestSample[s, ] <- match(ZtestSample[s, ], perm)
      }
    }
    probZMean <- apply(probZSample, c(2, 3), mean)
    MuMean <- apply(MuSample, c(2, 3, 4), mean)
    KappaMean <- lapply(probKappaSample, function(x) {
      aperm(apply(x, c(2, 3, 4), mean), c(2, 1, 3))
    })
    # Covariate coefficients have no cluster axis. Predictions are unchanged
    # when cluster labels and their coefficients are permuted together.
    print("Label switching addressed")
  }

  result <- list(
    DIC = res$DIC,
    probZMean = probZMean, MuMean = MuMean, KappaMean = KappaMean,
    sig2yM = res$sig2yMean, TransMatMean = TransMatMean,ZtestSample = ZtestSample,
    YpredTestSample = YpredTestSample,CoefySample = CoefySample,
    Betasample = Betasample, logpostsample = res$logpostsample
  )
  if (return_samples) {
    result <- c(result, list(
      ZSample = ZSample, probZSample = probZSample, MuSample = MuSample, probKappaSample = probKappaSample
    ))
  }
  result
}

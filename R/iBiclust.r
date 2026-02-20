#' Bayesian Supervised Biclustering via Hidden Markov Models for Multi-Modal Data
#'
#' Fits a Bayesian model-based supervised biclustering model for multi-modality
#' datasets using an Gibb sampling algorithm. Each modality/platform is provided as a matrix
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
#' @param Method Integer specifying the modeling approach.
#'   If \code{Method = 1}, the model incorporates a Hidden Markov Model (HMM)
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
#'  parameters \eqn{(\delta_1, \delta_2, \delta_3)} for each row of the transition
#'  matrix \eqn{\xi}.
#' @param TruncVal Positive numeric truncation threshold \eqn{t} used in the Gibbs sampler.
#' @param dataListTest Optional test set list with same structure as \code{dataList}.
#' @param CovTest Optional test covariate matrix of dimension \eqn{n_{test} \times p_c}.
#'
#' @returns A list with posterior summaries and samples:
#' \describe{
#'   \item{DIC}{Deviance Information Criterion returned by the sampler.}
#'   \item{probZMean}{Matrix of posterior mean cluster membership probabilities (\eqn{n \times K}).}
#'   \item{MuMean}{Array of posterior means of the latent mean parameters \eqn{\mu_{dkl}} for the over- and under-expressed
#'   states (\eqn{l = 1,2}), indexed by modality \eqn{d} and cluster \eqn{k}.}
#'   \item{CoefySample}{cluster-specific outcome means with dimension \eqn{sample \times K}.}
#'   \item{sig2yM}{Posterior mean of outcome variance parameter.}
#'   \item{KappaMean}{List of arrays of posterior means for biclustering indicators per modality.}
#'   \item{TransMatMean}{Array of posterior mean transition matrices per modality.}
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
#' p <- c(30, 50)   # smaller for an example
#' n <- 50
#' nbrcov <- 4
#' OutcType <- "continuous"
#'
#' dat <- iGenDataHMM(
#'   OrderFeatured = TRUE, K = K, NbrPlat = 2, OutcType = OutcType,
#'   MeanY = c(rep(1, floor(K/2)), rep(-1, K - floor(K/2))),
#'   pc = nbrcov, betacov = 1.0,
#'   MeanOverExpre = (1:K) / 2, MeanUnderExpre = -(1:K) / 2,
#'   p = p, seed = 1, n = n,
#'   SigmaOE = rep(1, K), SigmaUE = rep(1, K)
#' )
#'
#' ## Generate test data
#' ntest <- 20
#' datTest <- iGenDataHMM(
#'   OrderFeatured = TRUE, K = K, NbrPlat = 2, OutcType = OutcType,
#'   MeanY = c(rep(1, floor(K/2)), rep(-1, K - floor(K/2))),
#'   pc = nbrcov, betacov = 1.0,
#'   MeanOverExpre = (1:K) / 2, MeanUnderExpre = -(1:K) / 2,
#'   p = p, seed = 2, n = ntest,
#'   SigmaOE = rep(1, K), SigmaUE = rep(1, K)
#' )
#'
#' ## Fit model (keep MCMC small in examples)
#' result <- iBiclust(
#'   Method = 1,
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
    Method = 1, dataList = dataList, OutcType = "continuous",
    Youtcome = Youtcome, Cov = NULL, sample = 2000, burnin = 1000, K = 2,
    normalized = TRUE, seed = 1, hyperSigma = c(1, 1), alphaP = rep(1, K),
    delt = rep(1, 3), TruncVal = .1, dataListTest = NULL, CovTest = NULL) {
   # dyn.load("~/Desktop/~/RA_Thierry/bicluster/CodeBiolo/IBiclust.so")
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
  res <- .C("MCMCfunction",
            Method = as.integer(Method), OutcomeType = as.integer(OutcomeType),
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
            Ztest = as.integer(rep(0, sample * ntest))
  )



  Betasample1 <- t(matrix(res$Betasample, nrow = pc))

  # SampleCoefy=matrix(res$SampleCoefy,sample,K,byrow = T)

  KappaMean <- list()
  MuMean <- aperm(array(res$MuMean1, c(2, K, NbrPlat))) # is the dimension flipped?
  TransMatMean <- aperm(array(res$TransMatMean1, c(3, 3, NbrPlat)))

  m1 <- 1
  for (t in 1:NbrPlat)
  {
    probk <- res$probkappajk1[m1:(m1 - 1 + P[t] * K * 3)]
    KappaMean[[t]] <- aperm(array(probk, c(3, P[t], K)))
    m1 <- m1 + P[t] * K * 3
  }
  probZMean <- matrix(res$probZMean1, n, byrow = TRUE)

  list(
    DIC = res$DIC, probZMean = probZMean, MuMean = MuMean,
    CoefySample = t(matrix(res$CoefYsample, nrow = K)), sig2yM = res$sig2yMean,
    KappaMean = KappaMean, TransMatMean = TransMatMean,
    Betasample = Betasample1,
    ZtestSample = t(matrix(res$Ztest, nrow = ntest)),
    YpredTestSample = t(matrix(res$YpredTest, nrow = ntest)),
    logpostsample = res$logpostsample
  )
}

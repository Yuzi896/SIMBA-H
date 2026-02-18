#' Evaluate Test-Set Predictions
#'
#' Computes an outcome-specific performance metric using posterior predictive
#' samples returned by \code{\link{iBiclust}}.
#'
#' @param result A fitted model object returned by \code{\link{iBiclust}} containing
#'   \code{YpredTestSample} (posterior predictive samples for test outcomes).
#' @param TestY Test-set response variable. For continuous or binary outcomes,
#'   this should be a numeric vector. For censored outcomes, this must be an
#'   \eqn{n \times 2} matrix of the form \code{cbind(time, censor)}, where the
#'   second column is a censoring indicator (1 = censored, 0 = event observed).
#' @param OutcType Outcome type: one of \code{"continuous"}, \code{"binary"}, or \code{"censored"}.
#' @param threshold Numeric; classification threshold for binary outcomes applied to the
#'   posterior predictive mean. Default is 0.
#' @param return_all Logical; if \code{TRUE}, returns a list including intermediate
#'   predictions. If \code{FALSE}, returns only the metric(s).
#'
#' @returns
#' For \code{OutcType = "continuous"}: a scalar mean square error \code{MSE}.
#' For \code{OutcType = "binary"}: a scalar misclassification error \code{MCE}.
#' For \code{OutcType = "censored"}: a scalar concordance index \code{Cindex}.
#'
#' @examples
#' \dontrun{
#' ## Assume 'result' is a fitted iBiclust model with test predictions
#' ## and 'datTest$Y' is the corresponding test response.
#'
#' metric <- iEvaluateTest(
#'   result = result,
#'   TestY = datTest$Y,
#'   OutcType = "continuous"
#' )
#'
#' metric
#' }
#'
#' @export
EvaluateTest <- function(result, TestY, OutcType,
                          threshold = 0, return_all = FALSE) {

  if (is.null(result$YpredTestSample)) {
    stop("`result$YpredTestSample` is missing. Fit the model with test data to enable prediction evaluation.")
  }
  if (is.null(TestY)) {
    stop("`TestY` is missing.")
  }

  # Posterior predictive summary per test sample
  pred_mean <- apply(result$YpredTestSample, 2, mean)

  if (OutcType == "continuous") {
    y_true <- as.numeric(TestY)

    # TRUE mean squared error
    mse <- mean((pred_mean - y_true)^2)

    if (return_all) return(list(MSE = mse, pred_mean = pred_mean))

    return(mse)
  }

  if (OutcType == "binary") {
    y_true <- as.integer(TestY)
    y_hat <- as.integer(pred_mean > threshold)

    # Misclassification error
    mce <- mean(y_hat != y_true)

    if (return_all) return(list(MCE = mce, pred_mean = pred_mean, pred_class = y_hat))
    return(mce)
  }

  if (OutcType == "censored") {
    # censored
    y_true <- TestY
    if (!is.matrix(y_true) || ncol(y_true) != 2) {
      stop("For censored outcomes, `TestY` must be an n x 2 matrix: cbind(time, status).")
    }
    t_test <- as.numeric(y_true[, 1])
    delta_test <- as.integer(y_true[, 2])

    cindex <- Cindex(Tpredict = pred_mean, Ttest = t_test, Status = 1-delta_test)

    if (return_all) return(list(Cindex = cindex, pred_mean = pred_mean))
    return(cindex)
  }
}

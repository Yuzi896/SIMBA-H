#' Concordance Index (C-index) for Censored Outcomes
#'
#' Computes the concordance index using predicted times, observed times, and
#' censoring indicators.
#'
#' @param Tpredict Numeric vector of predicted times.
#' @param Ttest Numeric vector of observed times.
#' @param Status Integer/numeric vector of event indicators
#'   (1 = event observed, 0 = censored).

#' @returns Numeric scalar concordance index.
#' @export
#'
Cindex <- function(Tpredict, Ttest, Status) {
  N <- length(Tpredict)
  ConcordDeno <- 0
  ConcordNum <- 0

  for (i in 1:N) {
    time1 <- Ttest[i]
    timePred1 <- Tpredict[i]
    for (j in 1:N) {
      if (i != j) {
        time2 <- Ttest[j]
        timePred2 <- Tpredict[j]

        ConcordNum <- ConcordNum +
          (timePred2 > timePred1) * (time2 > time1) * (Status[i] == 1) +
          (timePred2 < timePred1) * (time2 < time1) * (Status[j] == 1) +
          0.5 * ((timePred2 == timePred1) || (time2 == time1)) * (Status[i] == 1) * (Status[j] == 0) +
          0.5 * ((timePred2 == timePred1) || (time2 == time1)) * (Status[j] == 1) * (Status[i] == 0)

        ConcordDeno <- ConcordDeno +
          (time2 > time1) * (Status[i] == 1) +
          (time2 < time1) * (Status[j] == 1) +
          (time2 == time1) * (Status[i] == 1) * (Status[j] == 0) +
          (time2 == time1) * (Status[i] == 0) * (Status[j] == 1)
      }
    }
  }

  ConcordNum / ConcordDeno
}

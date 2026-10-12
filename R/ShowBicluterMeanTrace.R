#' Plot MCMC Traces of Bicluster Means
#'
#' Creates trace plots for the overexpressed and underexpressed bicluster means
#' retained by \code{\link{iBiclust}}. The plot shows retained iteration on the
#' horizontal axis and mean parameter value on the vertical axis, with colors
#' identifying clusters and panels identifying modalities and expression states.
#'
#' @param result Results returned by \code{\link{iBiclust}} with
#'   \code{return_samples = TRUE}. Its \code{MuSample} component must be an array
#'   with dimensions retained iterations by modalities by clusters by expression
#'   states, where state 1 is overexpressed and state 2 is underexpressed.
#'
#' @returns A list with two components:
#' \describe{
#'   \item{data}{A long-format data frame containing the retained iteration
#'     \code{id}, cluster \code{k}, modality label \code{d}, expression state
#'     \code{name}, sampled mean \code{value}, modality/cluster label
#'     \code{group}, and \code{chain} (set to 1).}
#'   \item{figure}{A \code{ggplot} object. Use \code{print()} to display it.}
#' }
#'
#' @importFrom dplyr %>% left_join mutate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line facet_grid labs theme_bw theme
#'   element_rect element_text element_line
#' @export
ShowBicluterMeanTrace <- function(result){
  if(is.null(result$MuSample)){
    stop("Bicluster means not found. Fit iBiclust with return_samples = TRUE.")
  }
  dims <- dim(result$MuSample)
  if (!is.numeric(result$MuSample) || length(dims) != 4L ||
      any(dims < 1L) || dims[4] != 2L) {
    stop("MuSample must be a numeric array with dimensions iterations by modalities by clusters by two expression states.")
  }
  overexpressed <- c()
  underexpressed <- c()
  N <- dim(result$MuSample)[1]
  D <- dim(result$MuSample)[2]
  K <- dim(result$MuSample)[3]
  for(d in 1:D){
    for(k in 1:K){
      overexpressed <- rbind.data.frame(overexpressed,cbind.data.frame(over=result$MuSample[,d,k,1],k=k,d=d,id=c(1:N)))
      underexpressed <- rbind.data.frame(underexpressed,cbind.data.frame(under=result$MuSample[,d,k,2],k=k,d=d,id=c(1:N)))
    }
  }

  mcmc1 <- underexpressed%>%left_join(overexpressed, by = c("k", "d", "id"))%>%
    pivot_longer(c("over", "under"))%>%
    mutate(group = paste0("D = ", .data$d, " K = ", .data$k),
           name = ifelse(.data$name == "over", "Overexpressed","Underexpressed"),
           d = paste0("D = ", .data$d),
           chain = 1)

  image <- mcmc1%>%
    ggplot(aes(.data$id, .data$value))+
    geom_line(aes(color=factor(.data$k))) + facet_grid(d~name) + labs(x="Iteration", color="Biclusters")+
    theme_bw(base_size = 20) +
    theme(
      strip.background = element_rect(
        fill = "grey80", color = "grey65", linewidth = 0.4
      ),
      strip.text = element_text(size = 20, color = "black"),
      panel.border = element_rect(color = "grey65", linewidth = 0.4),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      panel.spacing = grid::unit(0.15, "cm"),
      axis.text = element_text(size = 20, color = "grey20"),
      axis.title = element_text(size = 20),
      legend.position = "bottom"
    )

  return(list(data = mcmc1, figure = image))
}

#' Theil-Sen Slope Estimator
#'
#' Estimate the slope of a time-series by the
#' median of all pairwise 1-step slopes.
#'
#' @param x Vector of values
#' @param t Optional vector of times, numerical
#'
#' @export
#' @useDynLib ConMK

theil_sen <- function(x, t = 1:length(x), ...) {
  delta <- c_theil_sen_vector(as.numeric(x), t)
  median(delta, na.rm = TRUE)
}

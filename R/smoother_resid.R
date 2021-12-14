
#' Residuals of the Smoother
#'
#' Calculate the residuals of the smoother to the data points.
#'
#' @inheritParams edecob
#' @param smoother_pts A data frame containing the smoother. Preferably the
#'   output of one of the smoother functions included in this package.
#'
#' @return A vector of the same length as \code{data} containing the residuals.
#' @export
smoother_resid <- function(data, smoother_pts) {

  resid <- numeric(nrow(data))

  if (nrow(data) > 0 && nrow(smoother_pts) > 0) {
    # calculate residuals for those time_points where there is both a
    # smoother and a data point
    for (ii in 1:nrow(smoother_pts)) {
      dataset_ind <- as.logical(data$time_point == smoother_pts$time_point[ii])
      resid[dataset_ind] <- data$value[dataset_ind] - smoother_pts$value[ii]
    }

    return(resid)

  } else {
    warning("Cannot calculate residuals (no data or smoother points)")
  }
}

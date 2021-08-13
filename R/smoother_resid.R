
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
#'
#' @examples
smoother_resid <- function(data, study_day, smoother_pts) {

  resid <- numeric(length(data))

  if (length(data) > 0 && nrow(smoother_pts) > 0) {
    # calculate residuals for those study_days where there is both a
    # smoother and a data point
    for (ii in 1:nrow(smoother_pts)) {
      dataset_ind <- as.logical(study_day == smoother_pts$study_day[ii])
      resid[dataset_ind] <- data[dataset_ind] - smoother_pts$pts[ii]
    }

    return(resid)

  } else {
    warning("Cannot calculate residuals (no data or smoother points)")
  }
}

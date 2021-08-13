
# Moving Median; Helper function for the bootstrap step.
fun_running_median <- function(data, win_beg_day, bt_rep, bt_Y, smoother_pts, width){


  window_ind <- as.logical((data$study_day >= win_beg_day) *
                           (data$study_day < win_beg_day + width))

  if (sum(window_ind) > 0) {
    return(list("value" = stats::median(bt_Y[window_ind], na.rm = TRUE),
                "study_day" = win_beg_day + width / 2,
                "bt_rep" = bt_rep)
    )
  } else {
    return(list("value" = NA,
                "study_day" = NA,
                "bt_rep" = bt_rep))
  }
}


# Perform one Bootstrap Step
#
# Helper function to bootstrap the epsilon (error of AR) and then reconstruct
# smoother (currently the moving median) using AR model and bootstrapped
# epsilon
bt_eps <- function(data, bt_rep, smoother_pts, resid, width){

  # bootstrap error
  bt_Y <- numeric(nrow(data))
  bt_eta <- numeric(nrow(data))
  bt_epsilon <- numeric(nrow(data))

  # fit autoregression model on residuals
  data_ind <- !is.na(resid)
  if (sum(data_ind) > 1 && stats::var(resid[data_ind]) != 0) {
    ar_resid <- stats::ar(resid[data_ind])

    # remove NAs and bootstrap the residuals of the AR model
    ar_resid_non_na <- ar_resid$resid[which(!is.na(ar_resid$resid))]
    bt_epsilon <- sample(ar_resid_non_na, length(ar_resid$resid), replace = T)

    # calculate residuals of the smoother with bootstrapped errors
    bt_eta <- sapply(1:length(data$study_day), function(x, ar_resid, resid) {
      # if first k data points where k order of AR (not enough previous data pts for AR)
      if (x <= ar_resid$order) {
        return(resid[min(which(!is.na(resid))) - 1 + x] + bt_epsilon[x])
      } else {
        pred_helper <- rep(NA, ar_resid$order)
        for (kk in 1:ar_resid$order) {
          pred_helper[kk] <-
            resid[min(which(!is.na(resid))) - 1 + x - kk]
        }
        return(stats::predict(ar_resid, newdata = pred_helper)$pred[[1]] + bt_epsilon[x])
      }
    }, ar_resid, resid)

    # calculate Y* (the bootstrapped data points using the AR model)
    bt_Y <- vapply(1:max(which(data$study_day <= max(smoother_pts$study_day))), function(x, bt_eta) {
      return(smoother_pts$pts[smoother_pts$study_day == data$study_day[x]] + bt_eta[x])
    }, numeric(1), bt_eta)

    # calculate S_star (the bootstrapped median)
    win_beg_day <- min(data$study_day)
    last_data_study_day <- max(data$study_day)

    if (last_data_study_day > win_beg_day + width) {

      S_star_one_bt <- as.data.frame(do.call(rbind, lapply(
        seq(win_beg_day, last_data_study_day - width - 1, 1),
        fun_running_median, data, bt_rep, bt_Y, smoother_pts, width
      )))
      return(S_star_one_bt)
    } else {
      return(data.frame("value" = NA,
                        "study_day" = NA,
                        "bt_rep" = bt_rep))
    }
  } else {
    warning("AR model cannot be fitted (variance of residuals is 0 or residuals are NA)")
  }
}

#' Bootstrap the Smoother
#'
#' Fit an autoregressive model on the residuals of the smoother and bootstrap
#' the error of the autoregressive model. Reconstruct the smoother by adding the
#' bootstrapped error, the fitted autoregressive model, and the smoother to form
#' the bootstrapped smoother.
#'
#' An autoregressive (AR) model is used for the residuals of the smoother:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sum^{p}_{j = 1} \phiⱼ \eta(t - j) + \epsilon}
#' where \eqn{t} is the point in time, \eqn{Y(t)} the data point,
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the residual of the smoother, \eqn{p}
#' the order of the AR model, \eqn{\phiⱼ} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model.
#'
#' The bootstrap procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(tᵢ) = Y(tᵢ) - S(tᵢ)}.
#'   \item Fit an AR(p) model to \eqn{\eta(tᵢ)} to obtain the coefficients
#'     \eqn{\phi₁, \dots, \phiₚ} and residuals \eqn{\epsilon(tᵢ) = \eta(tᵢ) -
#'     \sum^{p}_{j = 1} \phiⱼ \eta(tᵢ - tᵢ₋ⱼ)}.
#'   \item Resample \eqn{\epsilon(tᵢ)*} from \eqn{\epsilon(tₚ₊₁), \dots,
#'     \epsilon(tₙ)} to obtain \deqn{Y(tᵢ)* = S(tᵢ) + \eta(tᵢ)*,} where  \deqn{\eta(tᵢ)* = \sum^{p}_{j=1} \phiⱼ \eta(tᵢ₋ⱼ)*+ \epsilon(tᵢ₋ⱼ)*.}
#'   \item Compute \eqn{S(.)* = g(Y(t₁), \dots, Y(tₙ))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times.
#' }
#'
#'
#' @param smoother_pts A data frame containing the smoother. Preferably the
#'   output of one of the smoother functions included in this package.
#' @param resid A vector containing the residuals of the smoother to the data
#'   points.
#' @inheritParams edecob
#'
#' @return A data frame containing the bootstrap repetitions of the smoother.
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#' @examples
bt_smoother <- function(data, smoother_pts, resid, width, bt_tot_rep) {
  bt_Y <- do.call(rbind, lapply(1:bt_tot_rep, bt_eps, data, smoother_pts, resid, width))
  bt_Y$study_day <- unlist(bt_Y$study_day)
  bt_Y$value <- unlist(bt_Y$value)
  bt_Y$bt_rep <- unlist(bt_Y$bt_rep)
  bt_Y <- bt_Y[!is.na(bt_Y$value), ]
  return(bt_Y)
}

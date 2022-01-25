
# Perform one Bootstrap Step
#
# Helper function to bootstrap the epsilon (error of AR) and then reconstruct
# smoother (currently the moving median) using AR model and bootstrapped
# epsilon
bt_eps <- function(bt_rep, data, smoother, resample_method, smoother_pts, resid, ar_resid, ...){

  # bootstrap error
  bt_eta <- numeric(nrow(data))
  bt_epsilon <- numeric(nrow(data))
  if (smoother == "mov_med") {
    med_win <- list(...)$med_win
  }

  data_ind <- !is.na(resid)
  bt_Y <- numeric(length(data_ind))
  time_points_no_na <- data$time_point[which(!is.na(resid))]
  #print(time_points_no_na)
  # if ar_restr then bt only from window
  if ("ar_restr" %in% names(list(...))) {
    ar_restr = list(...)$ar_restr

    if (ar_restr == "win") {
      for (ii in 1:length(data_ind)) {
        ii_time_point <- time_points_no_na[ii]
        win_ind <- as.logical((time_points_no_na <= ii_time_point + med_win[2]) *
                              (time_points_no_na >= ii_time_point + med_win[1]))
        if (stats::var(resid[win_ind]) != 0) {

          # if order given in function call
          if ("order" %in% names(match.call())) {
            order <- match.call()$order
            ar_resid <- stats::ar(resid[win_ind], aic = FALSE, order.max = order)
          } else {
            ar_resid <- stats::ar(resid[win_ind])
          }

          # remove NAs and bootstrap the residuals of the AR model
          ar_resid_non_na <- ar_resid$resid[which(!is.na(ar_resid$resid))]

          # resample epsilon
          bt_epsilon <- sample(ar_resid_non_na, length(ar_resid$resid), replace = T)

          # calculate residuals of the smoother with bootstrapped errors
          bt_eta <- sapply(1:length(resid[win_ind]), function(x, ar_resid, resid) {
            # if first k data points where k order of AR (not enough previous data pts for AR)
            if (x <= ar_resid$order) {
              return(resid[min(which(!is.na(resid))) - 1 + x] + bt_epsilon[x])
            } else {
              pred_helper <- rep(NA, ar_resid$order)
              for (kk in 1:ar_resid$order) {
                if (kk > 0) {
                  pred_helper[kk] <-
                    resid[min(which(!is.na(resid))) + x - kk]
                }
              }
              return(stats::predict(ar_resid, newdata = pred_helper)$pred[[1]] + bt_epsilon[x])
            }
          }, ar_resid, resid)

          # calculate Y* (the bootstrapped data points using the AR model)
          bt_Y_one <- vapply(1:length(time_points_no_na[win_ind][which(time_points_no_na[win_ind] <= max(smoother_pts$time_point))]), function(x, bt_eta) { # data$time_point[which(data$time_point <= max(smoother_pts$time_point))])
            return(smoother_pts$value[smoother_pts$time_point == time_points_no_na[win_ind][x]] + bt_eta[x])
          }, numeric(1), bt_eta)
          bt_Y[ii] <- bt_Y_one[which(time_points_no_na[win_ind] == ii_time_point)]
        } else {
          bt_Y[ii] <- NA
        }
      }

      # calculate bootstrapped smoother S*
      win_beg_day <- min(data$time_point)
      last_data_time_point <- max(data$time_point)

      if (smoother == "mov_med") {
        med_win <- list(...)$med_win

        if (last_data_time_point > win_beg_day + med_win[2]) {
          S_star_one_bt <-
            mov_med(data.frame(source = rep(data$source[1], length(bt_Y)),
                               time_point = time_points_no_na,#[which(time_points_no_na <= max(smoother_pts$time_point))],
                               value = bt_Y),
                    med_win)
          S_star_one_bt <- S_star_one_bt[, c("source", "time_point", "value")]
          S_star_one_bt$bt_rep <- rep(bt_rep, nrow(S_star_one_bt))
          return(S_star_one_bt)
        }
      }
    }
  }

  # remove NAs and bootstrap the residuals of the AR model
  ar_resid_non_na <- ar_resid$resid[which(!is.na(ar_resid$resid))]

  # resample epsilon depending on method
  if (resample_method == "all") {
    bt_epsilon <- sample(ar_resid_non_na, length(ar_resid_non_na), replace = T)
  }

  if (resample_method == "past") {
    bt_epsilon <- numeric(length(ar_resid_non_na))
    for (ii in 1:length(ar_resid_non_na)) {
      bt_epsilon[ii] <- sample(ar_resid_non_na[1:ii], 1)
    }
  }

  if (resample_method == "window") {
    bt_epsilon <- numeric(length(ar_resid_non_na))
    time_points_non_na <- data$time_point[which(!is.na(ar_resid$resid))]
    for (ii in 1:length(ar_resid_non_na)) {
      ii_time_point <- time_points_non_na[ii]
      resample_win_ind <- as.logical((time_points_non_na <= ii_time_point + med_win[2]) *
                                     (time_points_non_na >= ii_time_point + med_win[1]))
      bt_epsilon[ii] <- sample(ar_resid_non_na[resample_win_ind], 1)
    }
  }

  # calculate residuals of the smoother with bootstrapped errors
  bt_eta <- sapply(1:length(data$time_point), function(x, ar_resid, resid) {
    # if first k data points where k order of AR (not enough previous data pts for AR)
    if (x <= ar_resid$order) {
      return(resid[min(which(!is.na(resid))) - 1 + x] + bt_epsilon[x])
    } else {
      pred_helper <- rep(NA, ar_resid$order)
      for (kk in 1:ar_resid$order) {
        if (kk > 0) {
          pred_helper[kk] <-
            resid[min(which(!is.na(resid))) + x - kk]

        }
      }
      return(stats::predict(ar_resid, newdata = pred_helper)$pred[[1]] + bt_epsilon[x])
    }
  }, ar_resid, resid)

  # calculate Y* (the bootstrapped data points using the AR model)
  bt_Y <- vapply(1:length(data$time_point[which(data$time_point <= max(smoother_pts$time_point))]), function(x, bt_eta) {
    return(smoother_pts$value[smoother_pts$time_point == data$time_point[x]] + bt_eta[x])
  }, numeric(1), bt_eta)

  # calculate bootstrapped smoother S*
  win_beg_day <- min(data$time_point)
  last_data_time_point <- max(data$time_point)

  if (smoother == "mov_med") {
    med_win <- list(...)$med_win

    if (last_data_time_point > win_beg_day + med_win[2]) {
      S_star_one_bt <-
        mov_med(data.frame(source = rep(data$source[1], length(bt_Y)),
                           time_point = data$time_point[which(data$time_point <= max(smoother_pts$time_point))],
                           value = bt_Y),
                med_win)
      S_star_one_bt <- S_star_one_bt[, c("source", "time_point", "value")]
      S_star_one_bt$bt_rep <- rep(bt_rep, nrow(S_star_one_bt))
      return(S_star_one_bt)

    } else {
      return(data.frame("source" = data$source[1],
                        "time_point" = NA,
                        "value" = NA,
                        "bt_rep" = bt_rep
                        ))
    }
  }

}

#' Bootstrap the Smoother
#'
#' First, fit an autoregressive model on the residuals of the smoother. Then bootstrap
#' the errors of the autoregressive model. Afterwards, reconstruct the measurements
#' by adding the bootstrapped error, the autoregressive model, and the smoother.
#' We can again calculate the smoother using these reconstructed measurements to
#' obtain the bootstrapped smoother (which can later be used to construct the
#' simultaneous confidence bounds). For details see below.
#'
#' An autoregressive (AR) model is used for the residuals of the smoother:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sum^{p}_{j = 1} \phi_j \eta(t - j) + \epsilon}
#' where \eqn{t} is the point in time, \eqn{Y(t)} the data point,
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the residual of the smoother, \eqn{p}
#' the order of the AR model, \eqn{\phi_j} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model.
#'
#' The bootstrap procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(t_i) = Y(t_i) - S(t_i)}.
#'   \item Fit an AR(p) model to \eqn{\eta(t_i)} to obtain the coefficients
#'     \eqn{\phi_1, \dots, \phi_p} and residuals \eqn{\epsilon(t_i) = \eta(t_i) -
#'     \sum^{p}_{j = 1} \phi_j \eta(t_i - t_{i-j})}.
#'   \item Resample \eqn{\epsilon(t_i)*} from \eqn{\epsilon(t_{p+1}), \dots,
#'     \epsilon(t_n)} to obtain \deqn{Y(t_i)* = S(t_i) + \eta(t_i)*,} where  \deqn{\eta(t_i)* = \sum^{p}_{j=1} \phi_j \eta(t_{i-j})*+ \epsilon(t_{i-j})*.}
#'   \item Compute \eqn{S(.)* = g(Y(t_1), \dots, Y(t_n))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times.
#' }
#'
#'
#' @param smoother_pts A data frame containing the smoother with columns time_point
#'   and value. Preferably the output of one of the smoother functions within
#'   this package.
#' @param resid A vector of the same length as the number of rows of data
#'   containing the difference between the smoother and the
#'   measurements. Preferably the output of \code{\link{smoother_resid}}.
#' @inheritParams edecob
#'
#' @return A data frame containing the bootstrap repetitions of the smoother.
#'   The column are subject identifier, time point, value, and the bootstrap
#'   repetition the value corrsponds to.
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#'
bt_smoother <- function(data, smoother, resample_method, smoother_pts, resid, bt_tot_rep, ...) {

  bt_eta <- numeric(nrow(data))
  bt_epsilon <- numeric(nrow(data))
  if (smoother == "mov_med") {
    med_win <- list(...)$med_win
  }

  data_ind <- !is.na(resid)
  bt_Y <- numeric(length(data_ind))
  time_points_no_na <- data$time_point[which(!is.na(resid))]

  # fit AR model if possible
  if (sum(data_ind) > 1 && stats::var(resid[data_ind]) != 0) {

    # if order given in function call
    if ("order" %in% names(match.call())) {
      order <- match.call()$order
      ar_resid <- stats::ar(resid[data_ind], aic = FALSE, order.max = order)
    } else {
      ar_resid <- stats::ar(resid[data_ind])
    }

    # bootstrap
    bt_Y <- do.call(rbind, lapply(1:bt_tot_rep, bt_eps, data, smoother, resample_method, smoother_pts, resid, ar_resid, ...))
    bt_Y <- bt_Y[!is.na(bt_Y$value), ]
    return(bt_Y)

  } else {
    warning(paste("AR model cannot be fitted for \"", data$source[1], "\" (variance of residuals is 0 or residuals are NA)", sep = ""))#, immediate. = TRUE)
    return(data.frame("source" = data$source[1],
                      "time_point" = NA,
                      "value" = NA,
                      "bt_rep" = 1:bt_tot_rep))
  }
}

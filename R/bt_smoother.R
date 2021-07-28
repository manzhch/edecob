
fun_running_median <- function(win_beg_day, bt_rep, bt_Y, date, smoother_pts, width){


  window_ind <- as.logical((date >= win_beg_day) *
                           (date < win_beg_day + width))

  if (sum(window_ind) > 0) {
    return(list("value" = stats::median(bt_Y[window_ind], na.rm = TRUE),
                "date" = win_beg_day + width / 2,
                "bt_rep" = bt_rep)
    )
  } else {
    return(list("value" = NA,
                "date" = NA,
                "bt_rep" = bt_rep))
  }
}

# one bootstrap step
# bootstrap the epsilon (error of AR) and then reconstruct smoother
# (currently the moving median) using AR model and epsilon*

bt_eps <- function(bt_rep, smoother_pts, resid, date, width){

  # bootstrap error
  bt_Y <- numeric(length(date))
  bt_eta <- numeric(length(date))
  bt_epsilon <- numeric(length(date))

  # fit autoregression model on residuals
  data_ind <- !is.na(resid)
  if (sum(data_ind) > 1 && stats::var(resid[data_ind]) != 0) {
    ar_resid <- stats::ar(resid[data_ind])

    # remove NAs and bootstrap the residuals of the AR model
    ar_resid_non_na <- ar_resid$resid[which(!is.na(ar_resid$resid))]
    bt_epsilon <- sample(ar_resid_non_na, length(ar_resid$resid), replace = T)

    # calculate residuals of the smoother with bootstrapped errors
    bt_eta <- sapply(1:length(date), function(x, ar_resid, resid) {
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
    bt_Y <- vapply(1:max(which(date <= max(smoother_pts$date))), function(x, bt_eta) {
      return(smoother_pts$pts[smoother_pts$date == date[x]] + bt_eta[x])
    }, numeric(1), bt_eta)

    # calculate S_star (the bootstrapped median)
    win_beg_day <- min(date)
    last_data_date <- max(date)

    if (last_data_date > win_beg_day + width) {

      S_star_one_bt <- as.data.frame(do.call(rbind, lapply(
        seq(
          win_beg_day,
          last_data_date - width - as.difftime(1, units = "days"),
          as.difftime(1, units = "days")
        ),
        fun_running_median, bt_rep, bt_Y, date, smoother_pts, width
      )))
      return(S_star_one_bt)
    } else {
      return(data.frame("value" = NA,
                        "date" = NA,
                        "bt_rep" = bt_rep))
    }
  } else {
    warning("AR model cannot be fitted (variance of residuals is 0 or residuals are NA)")
  }
}

bt_smoother <- function(smoother_pts, resid, date, width, bt_tot_rep) {
  bt_Y <- do.call(rbind, lapply(1:bt_tot_rep, bt_eps, smoother_pts, resid, date, width))
  bt_Y$date <- as.Date(unlist(bt_Y$date), origin = "1970-01-01")
  bt_Y$value <- unlist(bt_Y$value)
  bt_Y$bt_rep <- unlist(bt_Y$bt_rep)
  bt_Y <- bt_Y[!is.na(bt_Y$value), ]
  return(bt_Y)
}

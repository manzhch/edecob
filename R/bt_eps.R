
fun_running_median <- function(win_beg_day, bt_rep, bt_Y, date, smoother_pts, med_win_size){


  window_ind <- as.logical(
    (date >= win_beg_day) *
      (date < win_beg_day + med_win_size))

  if (sum(window_ind) > 0) {
    return(list("value" = stats::median(bt_Y[window_ind]),
                "date" = win_beg_day + med_win_size / 2,
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

bt_eps <- function(bt_rep, smoother_pts, resid, ar_resid, date, med_win_size){

  # bootstrap error
  bt_Y <- numeric(length(date))
  bt_eta <- numeric(length(date))
  bt_epsilon <- numeric(length(date))

  # if an AR model was fitted (i.e. data had > 1 rows and non-zero variance)
  if (!is.na(ar_resid)) {
    bt_epsilon <- sample(ar_resid$epsilon, length(ar_resid$epsilon), replace = T)

    # calculate residuals* with bootstrapped errors
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
    bt_Y <- sapply(1:length(date), function(x, bt_eta) {
      return(smoother_pts$pts[as.logical(smoother_pts$date == date[x])] +
               bt_eta[x])
    }, bt_eta)

    # calculate S_star (the bootstrapped median)
    win_beg_day <- min(date)
    last_data_date <- max(date)

    if (last_data_date > win_beg_day + med_win_size) {
      S_star_one_bt <- as.data.frame(do.call(rbind, lapply(
        seq(
          win_beg_day,
          last_data_date - med_win_size - as.difftime(1, units = "days"),
          as.difftime(1, units = "days")
        ),
        fun_running_median, bt_rep, bt_Y, date, smoother_pts, med_win_size
      )))
      return(S_star_one_bt)
    } else {
      return(data.frame("value" = NA,
                        "date" = NA,
                        "bt_rep" = bt_rep))
    }
  }



}

bt_S <- function(resid, ar_resid, date, med_win_size, bt_tot_rep) {
  bt_Y <- do.call(rbind, lapply(1:bt_tot_rep, bt_eps, resid, ar_resid, date, med_win_size))
  return(bt_Y)
}

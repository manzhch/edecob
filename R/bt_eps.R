
# function for calculating median
fun_running_median <- function(win_beg_day, bt_rep, bt_Y, date, med_pts, med_win_size){

  last_data_date <- max(date)
  window_ind <- as.logical(
    (date >= win_beg_day) *
      (date < win_beg_day + med_win_size))
  # print(median(bt_Y[window_ind]))
  if (sum(window_ind) > 0) {
    return(list("value" = stats::median(bt_Y[window_ind]),
                "date" = win_beg_day + med_win_size / 2,
                "USUBJID" = my_patient,
                "bt_rep" = bt_rep,
                "diff_from_first_approx" = stats::median(bt_Y[window_ind]) -
                  med_pts[as.logical(seq(min(date),
                                        max(date) - med_win_size/2 + as.difftime(1, units = "days"),
                                        as.difftime(1, units = "days"))
                                     == win_beg_day + med_win_size / 2)])
    )
  } else {
    return(list("value" = NA,
                "date" = NA,
                "USUBJID" = my_patient,
                "bt_rep" = bt_rep,
                "diff_from_first_approx" = NA))
  }
}


# bootstrap the epsilon (error of AR) and then reconstruct smoother
# using AR model and epsilon*

bt_eps <- function(med_pts, resid, ar_resid, date, med_win_size, bt_rep){

  # bootstrap error
  bt_Y <- numeric(length(patient_residual_ind[pt_ar_ind]))
  bt_eta <- numeric(length(patient_residual_ind[pt_ar_ind]))
  bt_epsilon <- numeric(length(patient_residual_ind[pt_ar_ind]))
  # print(pt_ar_ind)

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
      return(
        med_pts[as.logical(seq(min(date),
                               max(date) - med_win_size/2 + as.difftime(1, units = "days"),
                               as.difftime(1, units = "days"))
                           == date[x])]
        + bt_eta[x])
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
        fun_running_median, bt_rep, bt_Y, date, med_pts, med_win_size
      )))
      return(S_star_one_bt)
    } else {
      return(data.frame("value" = NA,
                        "date" = NA,
                        "USUBJID" = my_patient,
                        "bt_no" = bt_no,
                        "diff_from_first_approx" = NA))
    }
  }



}

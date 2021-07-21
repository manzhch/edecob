
# takes the bootstrapped smoother and checks what the ratio is of
# bootstrapped smoothers in the pointwise alpha_p quantiles
ratio_in_ci <- function(alpha_p,
                        bt_smoother,
                        smoother_pts, # bt_smoother/bt_tot_rep and smoother_pts should have the same size
                        bt_tot_rep,
                        alpha = 0.05) {
  # initialize variables
  uniq_date <- unique(date)
  my_quantile_lower <- numeric(length(uniq_date))
  my_quantile_upper <- numeric(length(uniq_date))
  bt_smoother$init_est <- bt_smoother - rep(smoother_pts, bt_tot_rep)

  # calculate quantiles
  for (ii in 1:length(uniq_date)) {
    my_quantile_lower[ii] <-
      stats::quantile(bt_smoother$init_est[bt_smoother$date == uniq_date[[ii]]], alpha_p)
    my_quantile_upper[ii] <-
      stats::quantile(bt_smoother$init_est[bt_smoother$date == uniq_date[[ii]]],
               1 - alpha_p)
  }
  names(my_quantile_lower) <- uniq_date
  names(my_quantile_upper) <- uniq_date

  # determine for each data point whether it is between the quantiles
  bt_smoother$in_interval <- logical(nrow(bt_smoother))
  bt_smoother$in_interval <-
    as.logical(pmin(bt_smoother$init_est <= rep(my_quantile_upper, bt_tot_rep),
                    bt_smoother$value >= rep(my_quantile_lower, bt_tot_rep)))

  # determine for each bootstrap whether all data points are between the quantiles
  bt_smoother$bt_in_interval <- logical(nrow(bt_smoother))
  for (ii in 1:bt_tot_rep) {
    bt_smoother$bt_in_interval[bt_smoother$bt_no == ii] <-
      as.logical(min(bt_smoother$in_interval[bt_smoother$bt_no == ii]))
  }

  bt_in_interval_ratio <-
    sum(bt_smoother$bt_in_interval[bt_smoother$date == min(unlist(bt_smoother$date))]) /
    bt_tot_rep

  return(bt_in_interval_ratio - 1 + alpha)
}

# looking for the right alpha_p such that 95% of the bootstrap curves
# are within CI
find_ptw_alpha <- function(bt_smoother,
                           smoother_pts,
                           bt_tot_rep,
                           alpha){

  if (nrow(bt_smoother) > bt_tot_rep &&
      nrow(bt_smoother[!is.na(bt_smoother$value), ]) > 0 &&
      sign(ratio_in_ci(0, bt_smoother, smoother_pts, bt_tot_rep, alpha)) !=
      sign(ratio_in_ci(0.05, bt_smoother, smoother_pts, bt_tot_rep, alpha))) {

    return(stats::uniroot(ratio_in_ci, lower = 0, upper = 0.05,
                          bt_smoother, smoother_pts, date, bt_tot_rep, alpha))

  } else if (nrow(bt_smoother) > bt_tot_rep &&
             nrow(bt_smoother[!is.na(bt_smoother$value), ]) > 0 &&
             sign(ratio_in_ci(0, bt_smoother, smoother_pts, bt_tot_rep, alpha)) ==
             sign(ratio_in_ci(0.05, bt_smoother, smoother_pts, bt_tot_rep, alpha))) {
    # print(2)
    return(NA)
    # print("sign same")
  } else if (nrow(bt_smoother) > bt_tot_rep &&
             nrow(bt_smoother[!is.na(bt_smoother$value), ]) == 0) {
    return(NA)
    # print("many NA")
  } else {
    return(NA)
  }
}



conf_band <- function(bt_smoother,
                      smoother_pts,
                      bt_tot_rep,
                      alpha){

  #bt_smoother$date <- as.Date(unlist(bt_smoother$date), format = "%Y-%m-%d", origin = "1970-01-01")

  # calculate pointwise quantile
  dates <- sort(unlist(unique(bt_smoother$date)))
  my_quantile_lower <- numeric(length(dates))
  my_quantile_upper <- numeric(length(dates))

  ptw_alpha <- find_ptw_alpha(bt_smoother, smoother_pts, bt_tot_rep, alpha)

  # calculate quantiles
  if (nrow(bt_smoother) > 0) {
    both_quantiles <-
      sapply(split(bt_smoother, factor(unlist(bt_smoother$date))),
             function(x){
               if (sum(!is.na(x)) > 0 &&
                   !is.na(ptw_alpha) &&
                   !is.null(ptw_alpha$root)) {

                 return(c(stats::quantile(unlist(x$init_est), ptw_alpha$root),
                          stats::quantile(unlist(x$init_est), 1 - ptw_alpha$root)))
               } else {
                  return(NA)
               }
             })
  } else {
    both_quantiles <- NA
  }

  if (!is.na(both_quantiles)) {
    my_quantile_lower <- both_quantiles[1, ]
    my_quantile_upper <- both_quantiles[2, ]
  }

  if (!is.null(date)) {
    names(my_quantile_lower) <- as.Date(dates, format = "%Y-%m-%d", origin = "1970-01-01")
    names(my_quantile_upper) <- as.Date(dates, format = "%Y-%m-%d", origin = "1970-01-01")

  }


  for (ii in unique(as.character(sort(smoother_pts$date[
    smoother_pts$date %in% unique(bt_smoother$date)])))) {

    ii <- as.Date(ii, "%Y-%m-%d")
    current_ind <- smoother_pts$date == ii
    smoother_pts$interval_upper[current_ind] <-
      smoother_pts$pts[current_ind] + my_quantile_upper[as.character(ii)]
    smoother_pts$interval_lower[current_ind] <-
      smoother_pts$pts[current_ind] + my_quantile_lower[as.character(ii)]

  }
}


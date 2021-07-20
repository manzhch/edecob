"diff_from_first_approx" = stats::median(bt_Y[window_ind]) -
  med_pts[as.logical(seq(min(date),
                         max(date) - med_win_size/2 + as.difftime(1, units = "days"),
                         as.difftime(1, units = "days"))
                     == win_beg_day + med_win_size / 2)]
# difference between initial median and bootstrapped median

# takes the bootstrapped smoother and checks what the ratio is of
# bootstrapped smoothers in the pointwise alpha_p quantiles
ratio_in_ci <- function(alpha_p,
                        bt_S,
                        med_pts, # bt_S/bt_tot_rep and med_pts should have the same size
                        date,
                        bt_tot_rep,
                        alpha = 0.05) {
  # initialize variables
  uniq_date <- unique(date)
  my_quantile_lower <- numeric(length(uniq_date))
  my_quantile_upper <- numeric(length(uniq_date))
  init_est <- bt_S - rep(med_pts, bt_tot_rep)
  names(init_est) <- bt_S$date

  # calculate quantiles
  for (ii in 1:length(uniq_date)) {
    my_quantile_lower[ii] <-
      quantile(init_est[names(init_est) == uniq_date[[ii]]], alpha_p)
    my_quantile_upper[ii] <-
      quantile(init_est[names(init_est) == uniq_date[[ii]]],
               1 - alpha_p)
  }
  names(my_quantile_lower) <- uniq_date
  names(my_quantile_upper) <- uniq_date

  # determine for each data point whether it is between the quantiles
  bt_S$in_interval <- logical(nrow(bt_S))
  bt_S$in_interval <-
    as.logical(pmin(init_est <= rep(my_quantile_upper, bt_tot_rep),
                    bt_S$value >= rep(my_quantile_lower, bt_tot_rep)))

  # determine for each bootstrap whether all data points are between the quantiles
  bt_S$bt_in_interval <- logical(nrow(bt_S))
  for (ii in 1:bt_tot_rep) {
    bt_S$bt_in_interval[bt_S$bt_no == ii] <-
      as.logical(min(bt_S$in_interval[bt_S$bt_no == ii]))
  }

  bt_in_interval_ratio <-
    sum(bt_S$bt_in_interval[bt_S$date == min(unlist(bt_S$date))]) /
    bt_tot_rep

  return(bt_in_interval_ratio - 1 + alpha)
}


# looking for the right alpha_p
root_search <- vector(mode = "list", length = length(patient_list))
i <- 1
for (my_patient in patient_list) {
  print(my_patient)
  S_star_one_patient <- S_star[S_star$USUBJID == my_patient, ]
  if (nrow(S_star_one_patient) > bt_tot_rep &&
      nrow(S_star_one_patient[!is.na(S_star_one_patient$value), ]) > 0 &&
      # nrow(FUTT[FUTT$USUBJID == my_patient,]) > 0 &&
      sign(pointwise_error_search(0)) != sign(pointwise_error_search(0.05))) {
    # print(1)
    root_search[[i]] <- stats::uniroot(pointwise_error_search, lower = 0, upper = 0.05)
  } else if (nrow(S_star_one_patient) > bt_tot_rep &&
             nrow(S_star_one_patient[!is.na(S_star_one_patient$value), ]) > 0 &&
             sign(pointwise_error_search(0)) == sign(pointwise_error_search(0.05))) {
    # print(2)
    root_search[[i]] <- NA
    print("sign same")
  } else if (nrow(S_star_one_patient) > bt_tot_rep &&
             nrow(S_star_one_patient[!is.na(S_star_one_patient$value), ]) == 0) {

    root_search[[i]] <- NA
    print("many NA")
  } else {
    # print(3)
    root_search[[i]] <- NA
  }
  i <- i + 1
}


# CI ----------

# calculate CI
for (my_patient in patient_list) {
  # profvis({
  print(my_patient)
  S_star_one_patient <- S_star[S_star$USUBJID == my_patient, ]
  S_star_one_patient$date <- as.Date(unlist(S_star_one_patient$date), format = "%Y-%m-%d", origin = "1970-01-01")

  # calculate pointwise quantile
  dates <- sort(unlist(unique(S_star_one_patient$date)))
  my_quantile_lower <- numeric(length(dates))
  my_quantile_upper <- numeric(length(dates))


  # calculate quantiles
  if (nrow(S_star_one_patient) > 0){
    both_quantiles <- sapply(split(S_star_one_patient,
                                   factor(unlist(S_star_one_patient$date))),
                             function(x){
                               if (sum(!is.na(x)) > 0 &&
                                   !is.na(root_search[[which(patient_list == my_patient)]]) &&
                                   !is.null(root_search[[which(patient_list == my_patient)]]$root)) {
                                 # print(unlist(x$date))
                                 # print(unlist(x$value))
                                 return(c(quantile(unlist(x$diff_from_first_approx),
                                                   root_search[[which(patient_list == my_patient)]]$root),
                                          quantile(unlist(x$diff_from_first_approx),
                                                   1 - root_search[[which(patient_list == my_patient)]]$root)))
                               } else {
                                 return(NA)
                               }
                             })
  } else {
    both_quantiles <- NA
  }
  # both_quantiles <- sapply(1:length(dates), function(x){
  #   S_star_diff_from_first_approx_one_patient <- S_star$diff_from_first_approx[
  #     as.logical((S_star$date == dates[x]) * (S_star$USUBJID == my_patient))]
  #   if (sum(!is.na(S_star_diff_from_first_approx_one_patient)) > 0 &&
  #       !is.na(root_search[[which(patient_list == my_patient)]]) &&
  #       !is.null(root_search[[which(patient_list == my_patient)]]$root)) {
  #     return(c(quantile(unlist(S_star_diff_from_first_approx_one_patient),
  #                       root_search[[which(patient_list == my_patient)]]$root),
  #              quantile(unlist(S_star_diff_from_first_approx_one_patient),
  #                       1 - root_search[[which(patient_list == my_patient)]]$root)))
  #   } else {
  #     return(NA)
  #   }
  # })

  if(!is.na(both_quantiles)) {
    my_quantile_lower <- both_quantiles[1, ]
    my_quantile_upper <- both_quantiles[2, ]
  }

  # print(both_quantiles)
  # for (i in 1:length(dates)) {
  #   S_star_diff_from_first_approx_one_patient <- S_star$diff_from_first_approx[
  #     as.logical((S_star$date == dates[i]) * (S_star$USUBJID == my_patient))]
  #   if (sum(!is.na(S_star_diff_from_first_approx_one_patient)) > 0 && !is.null(root_search[[which(patient_list == my_patient)]]$root)) {
  #     my_quantile_lower[i] <- quantile(S_star_diff_from_first_approx_one_patient,
  #                                      root_search[[which(patient_list == my_patient)]]$root)
  #     my_quantile_upper[i] <- quantile(S_star_diff_from_first_approx_one_patient,
  #                                      1 - root_search[[which(patient_list == my_patient)]]$root)
  #   }
  # }
  if(!is.null(dates)){
    names(my_quantile_lower) <- as.Date(dates, format = "%Y-%m-%d", origin = "1970-01-01")
    names(my_quantile_upper) <- as.Date(dates, format = "%Y-%m-%d", origin = "1970-01-01")

  }


  for(i in unique(as.character(sort(median_points$date[
    median_points$date %in% unique(S_star$date)])))){
    i <- as.Date(i, "%Y-%m-%d")
    current_ind <- as.logical((median_points$date == i) * (median_points$USUBJID == my_patient))
    median_points$interval_upper[current_ind] <-
      median_points$QRSRESN[current_ind] + my_quantile_upper[as.character(i)]
    median_points$interval_lower[current_ind] <-
      median_points$QRSRESN[current_ind] + my_quantile_lower[as.character(i)]

  }
} #)
#

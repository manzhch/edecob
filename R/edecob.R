# function to do everything

edecob <- function(data,
                   date,
                   subj_id = "subj1",
                   bt_tot_rep,
                   learn_dur,
                   basel_dur,
                   event_min_dur,
                   thresh_diff = 0.1,
                   smoother = "mov_med",
                   width = as.difftime(12, units = "weeks"),
                   alpha = 0.05,
                   plot = TRUE) {

  # check if date and data length match
  if (length(data) != length(date)) {
    stop("Length mismatch between data and date.")
  }

  # remove data for learning period
  data_non_learn <- data[date >= date[1] + learn_dur]
  date_non_learn <- date[date >= date[1] + learn_dur]


  # calculate the smoother
  if (smoother == "mov_med") {
    smoother_pts <- mov_med(data_non_learn, date_non_learn, subj_id, width)
  }

  # calculate residuals of the smoother
  smoother_resid <- mov_med_resid(data_non_learn, date_non_learn, smoother_pts)

  # bootstrap the errors of the AR model fitted on the residuals
  bt_smoother <- bt_smoother(smoother_resid, date, width, bt_tot_rep)

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, alpha)

  # detect events using confidence bands

}

# function to do everything

edecob <- function(data,
                   date,
                   subj_id = "subj1",
                   bt_tot_rep,
                   learn_dur,
                   basel_dur,
                   event_min_dur = as.difftime(12, units = "weeks"),
                   thresh_diff = 0.1,
                   smoother = "mov_med",
                   width = as.difftime(12, units = "weeks"),
                   alpha = 0.05,
                   plot = TRUE) {

  # check if date and data length match
  if (length(data) != length(date)) {
    warning("Length mismatch between data and date. Truncating the longer vector.")

    if (length(data) > length(date)) {
      data <- data[1:length(date)]
    } else if (length(date) > length(data)) {
      date <- date[1:length(data)]
    }
  }

  # sort data by date
  data <- data[order(date)]
  date <- date[order(date)]

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
  bt_smoother <- bt_smoother(smoother_resid, date_non_learn, width, bt_tot_rep)

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, alpha)

  # detect events using confidence bands
  event <- detect_event(data, date, conf_band, learn_dur, basel_dur, event_min_dur, thresh_diff)

  # plot
  if (plot == T) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package \"gglot2\" needed for plots.", call. = FALSE)
    } else {
      #plot
    }
  }

  return(event)
}

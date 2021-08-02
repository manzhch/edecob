# function to do everything

#' Detects events in the data.
#'
#' @param data A vector containing the data.
#' @param date A vector containing the dates at which the data points were collected
#' @param subj_id A string containing an identifier for the data.
#' @param bt_tot_rep The amount of bootstrap repetitions to be performed.
#' @param learn_dur The amount of days where the data is to be ignored.
#' @param basel_dur The width of the window from which the median is taken to use as baseline.
#' @param event_min_dur The minimal duration the worsening needs to be sustained to define an event.
#' @param thresh_diff The percentage difference between the baseline and threshold. A negative number indicates
#' @param smoother Which smoother is to be used. Currently only "mov_med" for moving median available.
#' @param width The width of the window for the moving median.
#' @param alpha The confidence level for the simultaneous confidence bands.
#' @param plot Whether a plot should be generated.
#'
#' @return A list of four values: `event_detected` gives whether an
#'   event was detected, `event_detection_date` gives the date at which the
#'   event was detected, `event_duration` gives the duration the event is
#'   sustained, and `event_censored` gives whether the detected event is censored.
#' @export
#'
#' @examples
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
  bt_smoother <- bt_smoother(smoother_pts, smoother_resid, date_non_learn, width, bt_tot_rep)

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

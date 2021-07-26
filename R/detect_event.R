
detect_event <- function(data,
                         date,
                         conf_band,
                         learn_dur,
                         basel_dur,
                         event_min_dur,
                         thresh_diff = 0.1) {

  # calculate baseline and threshold
  basel <- stats::median(data[as.logical(
    (date >= min(date) + learn_dur) *
      (date < min(date) + learn_dur + basel_dur))])
  thresh <- basel * (1 - thresh_diff)


  # initialize variables
  event_detected <- NA
  event_detect_date <- NA
  event_dur <- as.difftime(0, units = "days")
  conf_band$is_below_thresh <- conf_band$intvl_upper < thresh


  # find sequences of consecutive days where CI below threshold
  below_thresh_runs <- with(rle(conf_band$is_below_thresh), {
    ok <- values == TRUE
    ends <- cumsum(lengths)[ok]
    starts <- ends - lengths[ok] + 1
    dur <- ends - starts + 1
    event <- dur >= as.numeric(event_min_dur, units = "days")
    cbind(starts, ends, dur, event)
  })


  # if events found
  if (sum(below_thresh_runs[, "event"]) > 0) {
    event_detected <- TRUE
    event_detection_date <-
      conf_band$date[below_thresh_runs[, "starts"][
        which.max(below_thresh_runs[, "event"])]]
    event_duration <-
      as.difftime(below_thresh_runs[, "dur"]
                  [which.max(below_thresh_runs[, "event"])],
                  units = "days")
  } else {
    event_detected <- FALSE
    event_detection_date <- max(data$date) # censoring date
    event_duration <- max(below_thresh_runs[, "dur"]) # longest sequence below threshold
  }

  output <- list(event_detected, event_detection_date, event_duration)
  return(output)

}

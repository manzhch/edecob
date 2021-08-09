
#' Detect Events
#'
#' Detect events using the confidence bounds. In the case where the threshold is
#' below the baseline, an event is detected if all the points of the upper bound of the
#' confidence band are below the threshold for \code{event_min_dur}
#' consecutive days.
#'
#' @inheritParams edecob
#' @param conf_band A data frame containing the confidence bounds. Ideally the output of \code{\link{conf_band}}.
#'
#' @return A list of four values: \describe{ \item{\code{event_detected}}{gives
#'   whether an event was detected}
#'   \item{\code{event_detection_study_day}}{gives the study_day at which the
#'   event was detected} \item{\code{event_duration}}{gives the duration the
#'   event is sustained} \item{\code{event_censored}}{gives whether the detected
#'   event is censored}
#'   }
#' @export
#'
#' @examples
detect_event <- function(data,
                         study_day,
                         conf_band,
                         learn_dur,
                         basel_dur,
                         event_min_dur,
                         thresh_diff = -0.1) {

  # calculate baseline and threshold
  basel <- stats::median(data[as.logical(
    (study_day >= min(study_day) + learn_dur) *
      (study_day < min(study_day) + learn_dur + basel_dur))])
  thresh <- basel * (1 + thresh_diff)


  # check if threshold above or below baseline
  if (thresh_diff > 0) {
    conf_band$is_below_thresh <- conf_band$upper > thresh
  } else if (thresh_diff < 0) {
    conf_band$is_below_thresh <- conf_band$upper < thresh
  } else if (thresh_diff == 0) {
    warning("Threshold equal to baseline. Detecting events if data points below threshold.")
    conf_band$is_below_thresh <- conf_band$upper < thresh
  }

  # if there is at least one time point at which the upper interval is below the threshold
  if (sum(conf_band$is_below_thresh) > 0) {

    # find sequences of consecutive days where CI below threshold
    below_thresh_runs <- with(rle(conf_band$is_below_thresh), {
      ok <- values == TRUE
      ends <- cumsum(lengths)[ok]
      starts <- ends - lengths[ok] + 1
      dur <- ends - starts + 1
      event <- dur >= event_min_dur
      cbind(starts, ends, dur, event)
    })
    print(rle(conf_band$is_below_thresh))
    print(below_thresh_runs)

    # if events found
    if (sum(below_thresh_runs[, "event"]) > 0) {
      event_detected <- TRUE
      event_detection_study_day <-
        conf_band$study_day[below_thresh_runs[, "starts"][
          which.max(below_thresh_runs[, "event"])]]
      event_duration <-
        below_thresh_runs[, "dur"][which.max(below_thresh_runs[, "event"])]
      event_censored <- (event_detection_study_day + event_duration - 1) >= max(conf_band$study_day)

    } else {
      event_detected <- FALSE
      event_detection_study_day <- max(study_day) # censoring study_day
      event_duration <- max(below_thresh_runs[, "dur"]) # longest sequence below threshold
      event_censored <- FALSE
    }
  } else {
    event_detected <- FALSE
    event_detection_study_day <- max(study_day) # censoring study_day
    event_duration <- 0 # longest sequence below threshold
    event_censored <- FALSE
  }

  output <- list(event_detected = event_detected,
                 event_detection_study_day = event_detection_study_day,
                 event_duration = event_duration,
                 event_censored = event_censored)

  output[[3]] <- unname(output[[3]])
  output[[4]] <- unname(output[[4]])

  return(output)

}

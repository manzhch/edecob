
#' Detect Events
#'
#' Detect events using the confidence bounds. An event is detected if all the
#' points of the upper or lower bound of the confidence band are below or above
#' the threshold for \code{min_change_dur} consecutive days.
#'
#' @inheritParams edecob
#' @param conf_band A data frame containing the confidence bounds. Ideally the output of \code{\link{conf_band}}.
#' @param basel The baseline for the data.
#' @param thresh The threshold for the data.
#'
#' @return A list of four values: \describe{ \item{\code{event_detected}}{gives
#'   whether an event was detected}
#'   \item{\code{event_onset}}{gives the study_day at which the
#'   event was detected} \item{\code{event_duration}}{gives the duration the
#'   event is sustained} \item{\code{event_stop}}{gives whether the detected
#'   event is censored}
#'   }
#' @export
#'
#' @examples
detect_event <- function(conf_band,
                         basel,
                         thresh,
                         min_change_dur) {


  # check if threshold above or below baseline
  if (thresh > basel) {
    conf_band$is_below_thresh <- conf_band$upper > thresh
  } else if (thresh < basel) {
    conf_band$is_below_thresh <- conf_band$upper < thresh
  } else if (thresh == basel) {
    warning("Threshold equal to baseline. Detecting events if data points below threshold.")
    conf_band$is_below_thresh <- conf_band$upper < thresh
  }

  gap_below_thresh <- data.frame(
    "study_day" = min(conf_band$study_day):max(conf_band$study_day),
    "value" = FALSE)

   gap_below_thresh$value[gap_below_thresh$study_day %in% conf_band$study_day] <-
    conf_band$is_below_thresh

  # if there is at least one time point at which the upper interval is below the threshold
  if (sum(gap_below_thresh$value) > 0) {

    # find sequences of consecutive days where CI below threshold
    below_thresh_runs <- with(rle(gap_below_thresh$value), {
      ok <- values == TRUE
      ends <- cumsum(lengths)[ok]
      starts <- ends - lengths[ok] + 1
      dur <- ends - starts + 1
      event <- dur >= min_change_dur
      cbind(starts, ends, dur, event)
    })

    # if events found
    if (sum(below_thresh_runs[, "event"]) > 0) {
      event_detected <- TRUE
      event_onset <-
        gap_below_thresh$study_day[below_thresh_runs[, "starts"][
          which.max(below_thresh_runs[, "event"])]]
      event_duration <-
        below_thresh_runs[, "dur"][which.max(below_thresh_runs[, "event"])]
      event_stop <- (event_onset + event_duration - 1) >= max(conf_band$study_day)

    } else {
      event_detected <- FALSE
      event_onset <- max(conf_band$study_day) # censoring study_day
      event_duration <- max(below_thresh_runs[, "dur"]) # longest sequence below threshold
      event_stop <- FALSE
    }
  } else {
    event_detected <- FALSE
    event_onset <- max(conf_band$study_day) # censoring study_day
    event_duration <- 0 # longest sequence below threshold
    event_stop <- FALSE
  }

  output <- list(event_detected = event_detected,
                 event_onset = event_onset,
                 event_duration = event_duration,
                 event_stop = event_stop)

  output[[3]] <- unname(output[[3]])
  output[[4]] <- unname(output[[4]])

  return(output)

}


#' Detect Events
#'
#' Detect events using the confidence bounds. An event is detected if all the
#' points of the upper or lower bound of the confidence band are below or above
#' the threshold for \code{min_change_dur} consecutive days.
#'
#' @inheritParams edecob
#' @param conf_band A data frame containing the confidence bounds. Ideally the output of \code{\link{conf_band}}.
#' @param detec_lower The lower detection bound in the same units as the values in data.
#' @param detec_upper The upper detection bound in the same units as the values in data.
#'
#' @return A list of four values: \describe{ \item{\code{event_detected}}{gives
#'   whether an event was detected}
#'   \item{\code{event_onset}}{gives the time_point at which the
#'   event was detected} \item{\code{event_duration}}{gives the duration the
#'   event is sustained} \item{\code{event_stop}}{gives whether the detected
#'   event is censored}
#'   }
#' @export
#'
#'
detect_event <- function(conf_band,
                         detec_lower,
                         detec_upper,
                         min_change_dur) {

  if (nrow(conf_band) == 0) {
    event_detected <- FALSE
    event_onset <- max(conf_band$time_point) # censoring time_point
    event_duration <- 0 # longest sequence below threshold
    event_stop <- FALSE
  } else {

  # check if threshold above or below baseline
  conf_band$upper_inside <- conf_band$upper < detec_upper
  conf_band$lower_inside <- conf_band$lower > detec_lower

  # in case of gap set upper_inside and lower_inside to FALSE
  gap_upper_inside <- data.frame(
    "time_point" = seq(min(conf_band$time_point), max(conf_band$time_point)),
    "value" = FALSE)
  gap_lower_inside <- data.frame(
    "time_point" = seq(min(conf_band$time_point), max(conf_band$time_point)),
    "value" = FALSE)

  gap_upper_inside$value[gap_upper_inside$time_point %in% conf_band$time_point] <-
    conf_band$upper_inside
  gap_lower_inside$value[gap_lower_inside$time_point %in% conf_band$time_point] <-
    conf_band$lower_inside


  gap_both_inside <- data.frame(
    "time_point" = seq(min(conf_band$time_point), max(conf_band$time_point)),
    "value" = as.logical(gap_upper_inside$value * gap_lower_inside$value))

  # if there is at least one time point at which the both intervals are inside the detection interval
  if (sum(gap_both_inside) > 0) {

    # find sequences of consecutive days where both bound inside
    both_inside_runs <- with(rle(gap_both_inside$value), {
      ok <- values == TRUE
      ends <- cumsum(lengths)[ok]
      starts <- ends - lengths[ok] + 1
      dur <- ends - starts + 1
      event <- dur >= min_change_dur
      cbind(starts, ends, dur, event)
    })

    # if events found
    if (sum(both_inside_runs[, "event"]) > 0) {
      event_detected <- TRUE
      event_onset <-
        gap_both_inside$time_point[both_inside_runs[, "starts"][
          which.max(both_inside_runs[, "event"])]]
      event_duration <-
        both_inside_runs[, "dur"][which.max(both_inside_runs[, "event"])]
      event_stop <- (event_onset + event_duration - 1) >= max(conf_band$time_point)

    } else {

      event_detected <- FALSE
      event_onset <- max(conf_band$time_point) # censoring time_point
      event_duration <- max(c(both_inside_runs[, "dur"], 0)) # longest sequence below threshold
      event_stop <- FALSE
    }
  } else {
    event_detected <- FALSE
    event_onset <- max(conf_band$time_point) # censoring time_point
    event_duration <- 0 # longest sequence below threshold
    event_stop <- FALSE
  }
  }
  output <- list(event_detected = event_detected,
                 event_onset = event_onset,
                 event_duration = event_duration,
                 event_stop = event_stop)

  output[[3]] <- unname(output[[3]])
  output[[4]] <- unname(output[[4]])

  return(output)

}

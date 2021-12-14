
# summary <- function(object) {
#   UseMethod("summary", object)
# }

#' Summarizing Event Detection Results
#'
#' \code{summary} method for class "\code{edecob}", gives a summary of an edecob object
#'
#' @param object An object of class "\code{edecob}" for which the summary will be shown.
#' @param ... Other arguments
#'
#' @return A short summary whether an event was detected and the parameters of the \code{edecob} function call.
#' @export
#'
summary.edecob <- function(object, ...) {
  cat("edecob object\n\n")
  cat("Call: ")
  print(object$call)
  cat("\n")
  cat("Event detected:", object$event$event_detected, "\n")
  if (object$event$event_detected) {
    cat("Event onset ", object$time_unit, ": ", object$event$event_onset, "\n", sep = "")
    if (object$event$event_duration == 1) {
      cat("Event duration: ", object$event$event_duration, " ", object$time_unit, "\n", sep = "")
    } else {
      cat("Event duration: ", object$event$event_duration, " ", object$time_unit, "s\n", sep = "")
    }
    cat("Event sustained until end of observation:", object$event$event_stop, "\n")
  } else {
      if (object$event$event_duration == 1) {
        cat("Longest sustained change: ", object$event$event_duration, " ", object$time_unit, "\n", sep = "")
      } else {
        cat("Longest sustained change: ", object$event$event_duration, " ", object$time_unit, "s\n", sep = "")
      }
  }
  cat("---\n")
  cat("Source", ": ", object$source, "\n", sep = "")
  cat("Lower detection bound:", object$detec_lower, "\n")
  cat("Upper detection bound:", object$detec_upper, "\n")
  if (object$min_change_dur == 1) {
    cat("Minimal duration of change for event detection: ", object$min_change_dur, " ", object$time_unit, "\n", sep = "")
  } else {
    cat("Minimal duration of change for event detection: ", object$min_change_dur, " ", object$time_unit, "s\n", sep = "")
  }
  cat("Bootstrap repetitions:", object$bt_tot_rep, "\n")
  cat("Confidence band level: ", object$conf_band_lvl*100, "%", "\n", sep = "")
}

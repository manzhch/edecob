
# summary <- function(object) {
#   UseMethod("summary", object)
# }

#' Summarizing Event Detection Results
#'
#' \code{summary} method for class "\code{edecob}", gives a summary of an edecob object
#'
#' @param object An object of class "\code{edecob}" for which the summary will be shown.
#'
#' @return A short summary whether an event was detected and the parameters of the \code{edecob} function call.
#' @export
#'
summary.edecob <- function(object) {
  cat("edecob object\n\n")
  cat("Call: ")
  print(object$call)
  cat("\n")
  cat("Event detected:", object$event$event_detected, "\n")
  if (object$event$event_detected) {
    cat("Event ", object$colnames[2], ": ", object$event$event_onset, "\n", sep = "")
    cat("Event duration: ", object$event$event_duration, " ", object$colnames[2], "s\n", sep = "")
    cat("Event sustained until end of observation:", object$event$event_stop, "\n")
  } else {
    cat("Longest sustained change: ", object$event$event_duration, " ", object$colnames[2], "s\n", sep = "")
  }
  cat("---\n")
  cat(object$colnames[1], ": ", object$data$subj_id[1], "\n", sep = "")
  cat("Baseline:", object$baseline, "\n")
  cat("Threshold:", object$threshold, "\n")
  cat("Minimal duration of change for event detection: ", object$min_change_dur, " ", object$colnames[2], "s\n", sep = "")
  cat("Bootstrap repetitions:", object$bt_tot_rep, "\n")
  cat("Confidence band level: ", object$conf_band_lvl*100, "%", "\n", sep = "")
}

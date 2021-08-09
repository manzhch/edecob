#' edecob: Event Detection Using Confidence Bounds
#'
#' Detect sustained worsening/improvement in time-dependent data.
#'
#' @section Main functions:
#'   \describe{
#'     \item{\code{edecob}}{Detects event in the data provided}
#'     \item{\code{edecob_plot}}{Generates a ggplot to visualize the data}
#'   }
#'
#' @section Other functions:
#'   \describe{
#'     \item{\code{mov_med}}{Calculates the moving median}
#'     \item{\code{smoother_resid}}{Calculates the residuals of the moving median}
#'     \item{\code{bt_smoother}}{Bootstraps the smoother}
#'     \item{\code{conf_band}}{Calculates the confidence bounds}
#'     \item{\code{detect_event}}{Detects events using the confidence bounds}
#' }
#' @docType package
#' @name edecob-package
NULL
#> NULL

#' Event Detection using Confidence Bounds
#'
#' Calculate a smoother of the data and bootstrap it to form simultaneous
#' confidence bounds. Define an event if the confidence bound is below the
#' threshold for a predefined amount of time.
#'
#' Runs the following functions in this order: \enumerate{
#'   \item The smoother depending on choice: \code{\link{mov_med}}
#'   \item \code{\link{smoother_resid}}
#'   \item \code{\link{bt_smoother}}
#'   \item \code{\link{conf_band}}
#'   \item \code{\link{detect_event}}
#' }
#'
#' @param data A vector containing the data.
#' @param study_day A vector containing the study days at which the data points
#'   were collected
#' @param subj_id A string containing an identifier for the data.
#' @param bt_tot_rep The amount of bootstrap repetitions to be performed.
#' @param learn_dur The amount of days where the data is to be ignored.
#' @param basel_dur The width of the window in days from which the median is
#'   taken to use as baseline.
#' @param event_min_dur The minimal duration the worsening needs to be sustained
#'   to define an event in days.
#' @param thresh_diff The percentage difference between the baseline and
#'   threshold. A negative number indicates a threshold below the baseline and
#'   vice versa.
#' @param smoother Which smoother is to be used. Currently only `mov_med` for
#'   moving median available.
#' @param width The width of the window for the moving median in days.
#' @param alpha The confidence level for the simultaneous confidence bands.
#' @param plot Whether a plot should be generated.
#'
#' @return A list of four values: \describe{ \item{\code{event_detected}}{gives
#'   whether an event was detected}
#'   \item{\code{event_detection_study_day}}{gives the study_day at which the
#'   event was detected} \item{\code{event_duration}}{gives the duration the
#'   event is sustained} \item{\code{event_censored}}{gives whether the detected
#'   event is censored} }
#'
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#' @examples
edecob <- function(data,
                   study_day,
                   subj_id = "subj1",
                   bt_tot_rep,
                   learn_dur,
                   basel_dur,
                   event_min_dur = 12*7,
                   thresh_diff = -0.1,
                   smoother = "mov_med",
                   width = 12*7,
                   alpha = 0.05,
                   plot = TRUE) {

  # check if study_day and data length match
  if (length(data) != length(study_day)) {
    warning("Length mismatch between data and study_day. Truncating the longer vector.")

    if (length(data) > length(study_day)) {
      data <- data[1:length(study_day)]
    } else if (length(study_day) > length(data)) {
      study_day <- study_day[1:length(data)]
    }
  }

  # sort data by study_day
  data <- data[order(study_day)]
  study_day <- study_day[order(study_day)]

  # remove data for learning period
  data_non_learn <- data[study_day >= study_day[1] + learn_dur]
  study_day_non_learn <- study_day[study_day >= study_day[1] + learn_dur]


  # calculate the smoother
  if (smoother == "mov_med") {
    smoother_pts <- mov_med(data_non_learn, study_day_non_learn, subj_id, width)
  }

  # calculate residuals of the smoother
  smoother_resid <- smoother_resid(data_non_learn, study_day_non_learn, smoother_pts)

  # bootstrap the errors of the AR model fitted on the residuals
  bt_smoother <- bt_smoother(smoother_pts, smoother_resid, study_day_non_learn, width, bt_tot_rep)

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, alpha)

  # detect events using confidence bands
  event <- detect_event(data, study_day, conf_band, learn_dur, basel_dur, event_min_dur, thresh_diff)

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

#' edecob: Event Detection Using Confidence Bounds
#'
#' Detect sustained change in time-dependent data.
#'
#' @section Main function:
#'   \describe{
#'     \item{\code{edecob}}{Detects events in the data provided and visualizes the data}
#'   }
#'
#' @section Other functions:
#'   \describe{
#'     \item{\code{mov_med}}{Calculates the moving median}
#'     \item{\code{smoother_resid}}{Calculates the residuals of the moving median}
#'     \item{\code{bt_smoother}}{Bootstraps the smoother}
#'     \item{\code{conf_band}}{Calculates the confidence bounds}
#'     \item{\code{detect_event}}{Detects events using the confidence bounds}
#'     \item{\code{edecob_plot}}{Generates a ggplot to visualize the data}
#' }
#' @docType package
#' @name edecob-package
NULL
#> NULL



#' Event DEtection using COnfidence Bounds
#'
#' Calculate a smoother of the data and bootstrap it to form simultaneous
#' confidence bounds. Define an event if the confidence bound is below or above
#' the threshold for a predefined amount of time.
#'
#' @section The Model:
#' An autoregressive (AR) model is used for the residuals of the smoother:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sum^{p}_{j = 1} \phiⱼ \eta(t - j) + \epsilon}
#' where \eqn{t} is the point in time, \eqn{Y(t)} the data point,
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the residual of the smoother, \eqn{p}
#' the order of the AR model, \eqn{\phiⱼ} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model.
#'
#' @section Bootstrap:
#' The bootstrap procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(tᵢ) = Y(tᵢ) - S(tᵢ)}.
#'   \item Fit an AR(p) model to \eqn{\eta(tᵢ)} to obtain the coefficients
#'     \eqn{\phi₁, \dots, \phiₚ} and residuals \eqn{\epsilon(tᵢ) = \eta(tᵢ) -
#'     \sum^{p}_{j = 1} \phiⱼ \eta(tᵢ - tᵢ₋ⱼ)}.
#'   \item Resample \eqn{\epsilon(tᵢ)*} from \eqn{\epsilon(tₚ₊₁), \dots,
#'     \epsilon(tₙ)} to obtain \deqn{Y(tᵢ)* = S(tᵢ) + \eta(tᵢ)*,} where  \deqn{\eta(tᵢ)* = \sum^{p}_{j=1} \phiⱼ \eta(tᵢ₋ⱼ)*+ \epsilon(tᵢ₋ⱼ)*.}
#'   \item Compute \eqn{S(.)* = g(Y(t₁), \dots, Y(tₙ))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times.
#' }
#'
#' @section Calculation of the Confidence Bounds:
#' The confidence bounds are calculated as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ qₓ(tᵢ), q₁₋ₓ(tᵢ) i = 1,\dots, N}
#'
#'     where \deqn{qₓ(tᵢ) = inf {u; P*[S(tᵢ)*ᵦ - S(tᵢ) \le u] \ge x} } is a
#'     pointwise bootstrap quantile, \eqn{S(tᵢ)*ᵦ} the bootstrapped smoother,
#'     and \eqn{N} the number of measurements.
#'   \item We vary the pointwise error \eqn{x} until \deqn{P*[qₓ(tᵢ) \le S(tᵢ)*ᵦ - S(tᵢ) \le q₁₋ₓ(tᵢ) ⩝ i = 1,\dots, N] ≈ 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[qₓ(tᵢ), q₁₋ₓ(tᵢ)]} is approximately \eqn{1-\alpha}.
#'   \item We let
#'   \deqn{ Iₙ(tᵢ) = [S(tᵢ) +  qₓ(tᵢ), S(tᵢ) + q₁₋ₓ(tᵢ)] ⩝ i = 1, \dots, N}
#'   confidence bounds. Then \eqn{{Iₙ(tᵢ); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
#'
#'}
#'
#' @param data A data frame containing the data. The values should be in the
#'   first column, the study days in the second column, and the subject
#'   identifier in the third column. The study days are expected to be integes
#'   numbers. Study day zero may be the randomization date, the first day
#'   after baseline, or whatever is preferred. The subject identifier is expected
#'   to be a string.
#' @param bt_tot_rep The number of bootstrap repetitions to be performed.
#' @param min_change_dur The minimal duration the change needs to be sustained
#'   before an event is detected. In other words, the smallest amount of
#'   consecutive days the change needs to persist to define an event. Given in
#'   days. Defaults to 84, i.e. 12 weeks.
#' @param thresh_diff The percentage difference between the baseline and
#'   threshold. A negative number indicates a threshold below the baseline and
#'   similarly with a positive number. If it is zero, the algorithm will detect events below the
#'   baseline. Defaults to zero. To detect events above the baseline with threshold
#'   equal to baseline, use a very small number, i.e. 1e-10.
#' @param smoother Which smoother is to be used. Use \code{mov_med} for the
#'   moving median. When using the moving median, the parameter \code{width} must
#'   be given specifying the size of the window over which the moving median is
#'   to be taken. Note that the width is the total size of the window, meaning
#'   that for a day x the measurements from day x - \code{width}/2 to x +
#'   \code{width}/2 will be used for the calculation of the median. Defaults to the moving median.
#' @param conf_band_lvl The confidence level for the simultaneous confidence
#'   bands. Defaults to 0.95.
#' @param plot Whether a plot should be generated. Defaults to FALSE.
#' @param basel_start The study day of first data point to be used for the
#'   baseline. Defaults to the first study day at which a data point is present.
#' @param basel_end The study day of last data point to be used for the
#'   baseline.
#' @param basel_method How the baseline should be calculated. Can take any of
#'   the following values: \code{median, mean, quantile}. When using quantile,
#'   an additional argument \code{quantile} needs to be given. Defaults to
#'   median.
#' @param thresh_method Whether the threshold is a percentage or absolute
#'   difference from the baseline. Use \code{percentage} or \code{absolute}.
#'   Defaults to absolute.
#' @param ... Additional parameters to be given to the function.
#'
#' @return A list of four values: \describe{ \item{\code{event_detected}}{gives
#'   whether an event was detected}
#'   \item{\code{event_detection_study_day}}{gives the study_day at which the
#'   event was detected. It is defined as the first day the upper or lower bound
#'   of the confidence band stays below or above the threshold, and after which it
#'   stays below or above the threshold for at least \code{min_change_dur} days.}
#'   \item{\code{event_duration}}{gives the duration the
#'   event is sustained, i.e. the number of days the upper or lower bound
#'   of the confidence band stays below or above the threshold.}
#'   \item{\code{event_stop}}{gives whether the detected change is sustained until
#'   the last day at which we can calculate the confidence bound or not. } }
#'
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#' @examples
edecob <- function(data,
                   bt_tot_rep,
                   basel_start = min(data[, 2]),
                   basel_end,
                   basel_method = "median",
                   thresh_diff = 0,
                   thresh_method = "absolute",
                   min_change_dur = 12*7,
                   smoother = "mov_med",
                   conf_band_lvl = 0.95,
                   plot = FALSE,
                   ...) {

  # EXTRACT ... !!!!!!!! may contain quantile or width.
  colnames(data) <- c("value, study_day, subj_id")

  # check if study_day and data length match
  # if (length(data) != length(study_day)) {
  #   warning("Length mismatch between data and study_day. Truncating the longer vector.")
  #
  #   if (length(data) > length(study_day)) {
  #     data <- data[1:length(study_day)]
  #   } else if (length(study_day) > length(data)) {
  #     study_day <- study_day[1:length(data)]
  #   }
  # }

  # sort data by study_day
  data <- data[order(data$study_day), ]

  # remove data for learning period
  data_non_learn <- data[data$study_day >= basel_start, ]

  # calculate baseline
  if (basel_method == "median") {
    basel <- stats::median(data$value[as.logical(
      (data$study_day >= basel_start) * (data$study_day <= basel_end))])
  } else if (basel_method == "mean") {
    basel <- mean(data$value[as.logical(
      (data$study_day >= basel_start) * (data$study_day <= basel_end))])
  } else if (basel_method == "quantile") {
    basel <- stats::quantile(data$value[as.logical(
      (data$study_day >= basel_start) * (data$study_day <= basel_end))], quantile)
  } else {
    warning("Baseline method not recognized. Default to median.")
    basel <- stats::median(data$value[as.logical(
      (data$study_day >= basel_start) * (data$study_day <= basel_end))])
  }

  # calculate threshold
  if (thresh_method == "percentage") {
    thresh <- basel * (1 + thresh_diff)
  } else if (thresh_method == "absolute") {
    thresh <- basel + thresh_diff
  } else {
    warning("Threshold method not recognized. Default to absolute.")
    thresh <- basel + thresh_diff
  }


  # calculate the smoother
  if (smoother == "mov_med") {
    smoother_pts <- mov_med(data_non_learn, width)
  }

  # calculate residuals of the smoother
  smoother_resid <- smoother_resid(data_non_learn, smoother_pts)

  # bootstrap the errors of the AR model fitted on the residuals
  bt_smoother <- bt_smoother(data_non_learn, smoother_pts, smoother_resid, width, bt_tot_rep)

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)

  # detect events using confidence bands
  event <- detect_event(conf_band, basel, thresh, min_change_dur)

  output <- list(
    "data" = data,
    "study_day" = study_day,
    "subj_id" = subj_id,
    "bt_tot_rep" = bt_tot_rep,
    "basel_start" = basel_start,
    "basel_end" = basel_end,
    "min_change_dur" = min_change_dur,
    "thresh_diff" = thresh_diff,
    "smoother" = smoother,
    "width" = width,
    "conf_band_lvl" = conf_band_lvl,
    "smoother_pts" = smoother_pts,
    "conf_band" = conf_band,
    "event" = event
  )

  # plot
  if (plot == T) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package \"gglot2\" needed for plots.", call. = FALSE)
    } else {
      #plot
      output[["plot"]] <- edecob_plot(data,
                                      study_day,
                                      basel,
                                      thresh,
                                      smoother_pts,
                                      conf_band,
                                      event,
                                      subj_id,
                                      basel_start,
                                      basel_end,
                                      width,
                                      label = "Data")
    }
  }

  class(output) <- "edecob"

  return(output)
}

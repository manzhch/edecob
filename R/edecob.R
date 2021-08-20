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
#' Calculate a smoother of the data and bootstrap the errors of the autoregressive
#' model fitted on the smoother to form simultaneous
#' confidence bounds. Define an event if the simultaneous confidence bound is
#' below or above a threshold for a predefined amount of time.
#'
#' @section Moving Median:
#'   Consider a sample \eqn{X₁,\dots, Xₙ} of size \eqn{n} and the
#'   reordering \eqn{X₍₁₎,\dots, X₍ₙ₎} such
#'   that \eqn{X₍₁₎ \le X₍₂₎ \le \dots \le X₍ₙ₎}, commonly
#'   called the order statistic. Then for \eqn{n} even the median usually
#'   defined as \deqn{median(X₁,\dots, Xₙ) = X₍ₖ₎, where k = n/2.} In the
#'   case where \eqn{n} is odd the median is
#'   defined as \deqn{median(X₁,\dots, Xₙ) = 1/2(X₍ₖ₎ + X₍ₖ₊₁₎), where k = n/2.} Let the
#'   study days at which the measurements \eqn{X₁, \dots, Xₙ} were taken
#'   be \eqn{t₁, \dots, tₙ}.
#'   Let \eqn{T} a fixed positive amount of time. Then the moving median at time
#'   point \eqn{t} with window size \eqn{T} is defined as
#'   \deqn{S(t) = median({Xⱼ | t - T/2 \le tⱼ \le t + T/2}).}
#'
#' @section The Model:
#' An autoregressive (AR) model is used to model the residuals of the smoother \eqn{\eta}:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sumᵖⱼ₌₁ \phiⱼ \eta(t - j) + \epsilon} where
#' variable \eqn{t} is the study day, \eqn{Y(t)} the data point at study day \eqn{t},
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the difference between the smoother
#' and the measurement at study day \eqn{t}, \eqn{p}
#' the order of the AR model, \eqn{\phiⱼ} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model. The order is calculated using the
#' Akaike information criterion (AIC).
#'
#' @section Bootstrap:
#' The bootstrap procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(tᵢ) = Y(tᵢ) - S(tᵢ)}.
#'   \item Fit an AR(p) model to \eqn{\eta(tᵢ)} to obtain the coefficients
#'     \eqn{\phi₁,\dots, \phiₚ} and \eqn{\epsilon(tᵢ) = \eta(tᵢ) -
#'     \sumᵖⱼ₌₁ \phiⱼ \eta(tᵢ - tᵢ₋ⱼ)} the error of the AR model.
#'   \item Resample \eqn{\epsilon(tᵢ)*} with replacement from \eqn{\epsilon(tₚ₊₁),\dots,
#'     \epsilon(tₙ)} to obtain \deqn{Y(tᵢ)* = S(tᵢ) + \eta(tᵢ)*,} where  \deqn{\eta(tᵢ)* = \sumᵖⱼ₌₁ \phiⱼ \eta(tᵢ₋ⱼ)*+ \epsilon(tᵢ₋ⱼ)*} the bootstrapped residuals of the smoother.
#'   \item Compute \eqn{S(.)* = g(Y(t₁),\dots, Y(tₙ))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times to obtain \eqn{S(tᵢ)*ᵦ} for \eqn{β = 1,\dots,}
#'   \code{bt_tot_rep}.
#' }
#'
#' @section Calculation of the Confidence Bounds:
#' The confidence bounds are calculated as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ qₓ(tᵢ), q₁₋ₓ(tᵢ) i = 1,\dots, N}
#'
#'     where \deqn{qₓ(tᵢ) = inf{u; P*[S(tᵢ)*ᵦ - S(tᵢ) \le u] \ge x} } is a
#'     pointwise bootstrap quantile, \eqn{S(tᵢ)*ᵦ} the bootstrapped smoother,
#'     and \eqn{N} the number of measurements or rows in \code{data}, in our case the number of rows.
#'   \item We vary the pointwise error \eqn{x} until \deqn{P*[qₓ(tᵢ) \le S(tᵢ)*ᵦ - S(tᵢ) \le q₁₋ₓ(tᵢ) ⩝ i = 1,\dots, N] ≈ 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[qₓ(tᵢ), q₁₋ₓ(tᵢ)]} is approximately \eqn{1-\alpha}.
#'   \item We define
#'   \deqn{ Iₙ(tᵢ) = [S(tᵢ) +  qₓ(tᵢ), S(tᵢ) + q₁₋ₓ(tᵢ)] ⩝ i = 1,\dots, N}
#'   the confidence bounds. Then \eqn{{Iₙ(tᵢ); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
#'
#'}
#'
#' @section Notes:
#' For the moving median, the width is the total size of the window, meaning
#' that for the value corresponding to day x, the data points from day
#' x - \code{width}/2 to x + \code{width}/2 will be used for the calculation
#' of the median.
#'
#' @seealso \code{\link{summary.edecob}}, \code{\link{plot.edecob}}
#'
#' @param data A data frame or matrix containing the data. The data points should be in the
#'   first column, the study days in the second column, and the subject
#'   identifier in the third column. If a matrix is used, the subject identifier
#'   is to be specified in the row names. The study days are expected to be integer
#'   numbers. Study day zero may be the randomization date, the first day
#'   after baseline, or whatever is preferred. The subject identifier is expected
#'   to be a string.
#' @param smoother Which smoother is to be used. Use \code{mov_med} for the
#'   moving median. When using the moving median, the parameter \code{width} must
#'   be given to specify the size of the window over which the moving median is
#'   to be taken. Defaults to the moving median.
#' @param bt_tot_rep The number of times we perform the bootstrap. Defaults to 100.
#' @param min_change_dur The minimal number of days that the confidence bounds
#'   need to stay below or above the threshold in order for an event to be
#'   detected. Defaults to 84, i.e. 12 weeks.
#' @param thresh_diff The absolute or percentage (depending on \code{thresh_method})
#'   difference between the baseline and
#'   threshold. A negative number indicates a threshold below the baseline and
#'   similarly a positive number indicates a threshold above the baseline. If it
#'   is zero, the algorithm will detect events below the baseline.
#'   Defaults to -20%. To detect events above the baseline with threshold
#'   equal to baseline, use a very small number like 1e-10.
#' @param conf_band_lvl The confidence level for the simultaneous confidence
#'   bands. Defaults to 0.95.
#' @param basel_start The study day of first data point to be used for the
#'   baseline. Defaults to the smallest available study day.
#' @param basel_end The study day of last data point to be used for the
#'   baseline. Defaults to six weeks after \code{basel_start}.
#' @param basel_method How the baseline should be calculated. Can take any of
#'   the following values: \code{median, mean, quantile}. When using quantile,
#'   an additional argument \code{quantile} specifying the quantile needs to be
#'   given. Defaults to median.
#' @param thresh_method Whether the threshold is a percentage or absolute
#'   difference from the baseline. Use \code{percent} or \code{absolute}.
#'   Defaults to percent.
#' @param ... Additional parameters to be given to the function. When the moving
#'   median is used, a \code{width} parameter is expected. When using the quantile
#'   for baseline method, a \code{quantile} parameter is expected.
#'
#' @return A list of four variables: \describe{ \item{\code{event_detected}}{gives
#'   whether an event was detected}
#'   \item{\code{event_detection_study_day}}{gives the study_day at which the
#'   event was detected. It is defined as the first day the upper or lower bound
#'   of the confidence band is below or above the threshold, and after which it
#'   stays below or above the threshold for at least \code{min_change_dur} days.}
#'   \item{\code{event_duration}}{gives the number of days the upper or lower bound
#'   of the confidence band stays below or above the threshold.}
#'   \item{\code{event_stop}}{gives whether the the confidence bounds stay below or
#'   above the threshold until
#'   the last day at which we can calculate the confidence bound or not. } }
#'
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#'   Hogg, R., McKean, J. and Craig, A. (2014).
#'   \emph{Introduction to mathematical statistics.} Harlow: Pearson Education.
#'
#' @examples
edecob <- function(data,
                   smoother = "mov_med",
                   basel_start = min(data[, 2]),
                   basel_end = basel_start + 41,
                   basel_method = "median",
                   thresh_diff = -20,
                   thresh_method = "percent",
                   min_change_dur = 12*7,
                   conf_band_lvl = 0.95,
                   bt_tot_rep = 100,
                   ...) {

  # EXTRACT ... !!!!!!!! may contain quantile or width.

  # may be matrix!!!
  colnames(data) <- c("value", "study_day", "subj_id")
  # print(data)
  width <- 12*7
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
  # data <- data[order(data$study_day), ]

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
    warning("Baseline method not recognized. Defaulting to median.")
    basel <- stats::median(data$value[as.logical(
      (data$study_day >= basel_start) * (data$study_day <= basel_end))])
  }

  # calculate threshold
  if (thresh_method == "percent") {
    thresh <- basel * (1 + thresh_diff/100)
  } else if (thresh_method == "absolute") {
    thresh <- basel + thresh_diff
  } else {
    warning("Threshold method not recognized. Defaulting to percent.")
    thresh <- basel * (1 + thresh_diff/100)
  }


  # calculate the smoother
  if (smoother == "mov_med") {
    smoother_pts <- mov_med(data_non_learn, width)
  }

  # calculate residuals of the smoother
  smoother_resid <- smoother_resid(data_non_learn, smoother_pts)


  # bootstrap the errors of the AR model fitted on the residuals
  if (smoother == "mov_med") {

    bt_smoother <- bt_smoother(data_non_learn, smoother, smoother_pts, smoother_resid, bt_tot_rep, width)

  }

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)

  # detect events using confidence bands
  event <- detect_event(conf_band, basel, thresh, min_change_dur)

  output <- list(
    "event" = event,
    "conf_band" = conf_band,
    "smoother_pts" = smoother_pts,
    "data" = data,
    "baseline" = basel,
    "threshold" = thresh,
    "smoother" = smoother,
    "basel_start" = basel_start,
    "basel_end" = basel_end,
    "basel_method" = basel_method,
    "thresh_diff" = thresh_diff,
    "thresh_method" = thresh_method,
    "min_change_dur" = min_change_dur,
    "conf_band_lvl" = conf_band_lvl,
    "bt_tot_rep" = bt_tot_rep
  )

  # plot
  # if (plot == T) {
  #   if (!requireNamespace("ggplot2", quietly = TRUE)) {
  #     warning("Package \"gglot2\" needed for plots.", call. = FALSE)
  #   } else {
  #     #plot
  #     output[["plot"]] <- edecob_plot(data,
  #                                     basel,
  #                                     thresh,
  #                                     smoother_pts,
  #                                     conf_band,
  #                                     event,
  #                                     basel_start,
  #                                     basel_end,
  #                                     width,
  #                                     label = "Data")
  #   }
  # }

  class(output) <- "edecob"

  return(output)
}

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
#' The term study day is not to be taken too seriously; it only needs to be a
#' number specifying the time point at which the measurement was taken.
#'
#' For the moving median, the width is the total size of the window, meaning
#' that for the value corresponding to day x, the data points from day
#' x - \code{width}/2 to x + \code{width}/2 will be used for the calculation
#' of the median.
#'
#' @seealso \code{\link{summary.edecob}}, \code{\link{plot.edecob}}
#'
#' @param data A data frame in long format containing the data for which events
#'   is to be detected. The first column
#'   contains stings specifying the subject identifier, the second column contains
#'   numbers specifying the study day, and the third column contains the numerical
#'   values of the measurements.
#' @param smoother Which smoother is to be used. Use \code{mov_med} for the
#'   moving median. When using the moving median, the parameter \code{width} must
#'   be given to specify the size of the window over which the moving median is
#'   to be taken. Defaults to the moving median.
#' @param baseline A number specifying the baseline.
#' @param threshold A number specifying the threshold. If the threshold is smaller
#'   than the than the baseline, an event will be detected
#'   when the confidence band stays below the threshold for \code{min_change_dur}
#'   days. Similarly, if the threshold is larger than baseline, an event will be detected
#'   when the confidence band stays above the threshold for \code{min_change_dur}
#'   days. If the threshold is equal to the baseline, events below the threshold
#'   will be detected. For detection of events above the threshold with the
#'   threshold equal to the baseline, choose a threshold minimally larger than
#'   the baseline.
#' @param bt_tot_rep The number of times we perform the bootstrap. Because of
#'   run time, it is recommended to keep this number below 500, especially with
#'   large data sets. Defaults to 100.
#' @param min_change_dur The minimal number of days that the confidence bounds
#'   need to stay below or above the threshold in order for an event to be
#'   detected. Defaults to 84, i.e. 12 weeks.
#' @param conf_band_lvl The confidence level for the simultaneous confidence
#'   bands. Defaults to 0.95.
#' @param ... Additional parameters to be given to the function. When the moving
#'   median is used, a \code{width} parameter is expected.
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
#'   edecob(data.frame(data = LakeHuron, year = tsp(LakeHuron)[1]:tsp(LakeHuron)[2], subject = "LakeHuron"))
edecob <- function(data,
                   smoother = "mov_med",
                   baseline,
                   threshold,
                   min_change_dur = 12*7,
                   conf_band_lvl = 0.95,
                   bt_tot_rep = 100,
                   ...) {

  # EXTRACT ... !!!!!!!! may contain width.
  colnames(data) <- c("subj_id", "study_day", "value")
  width <- 12*7

  # calculate the smoother
  if (smoother == "mov_med") {
    smoother_pts <- mov_med(data, width)
  } else {
    warning("Smoother not recognized. Defaulting to moving median.")

  }

  # calculate residuals of the smoother
  smoother_resid <- smoother_resid(data, smoother_pts)

  # bootstrap the errors of the AR model fitted on the residuals
  if (smoother == "mov_med") {
    bt_smoother <- bt_smoother(data, smoother, smoother_pts, smoother_resid, bt_tot_rep, width)
  }

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)

  # detect events using confidence bands
  event <- detect_event(conf_band, baseline, threshold, min_change_dur)

  output <- list(
    "event" = event,
    "conf_band" = conf_band,
    "smoother_pts" = smoother_pts,
    "data" = data,
    "baseline" = baseline,
    "threshold" = threshold,
    "smoother" = smoother,
    "min_change_dur" = min_change_dur,
    "conf_band_lvl" = conf_band_lvl,
    "bt_tot_rep" = bt_tot_rep
  )

  class(output) <- "edecob"

  return(output)
}

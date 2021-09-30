#' edecob: Event Detection Using Confidence Bounds
#'
#' Detect sustained change (or "events") in time-dependent data while accounting for
#' heterogeneity in assessment frequency and noise. The data is first smoothed
#' (e.g. using the moving median) before an autoregressive model is fitted on
#' the residuals of the smoother.
#' The errors of this autoregressive model are then bootstrapped to find
#' simultaneous confidence bands for the smoother.
#'
#' @section Main function:
#'   \describe{
#'     \item{\code{edecob}}{Detects sustained change in the data provided}
#'   }
#'
#' @section Other functions:
#'   \describe{
#'     \item{\code{plot.edecob}}{Generates a ggplot to visualize the output of \code{edecob}}
#'     \item{\code{summary.edecob}}{Gives a summary of the sustained change detected}
#'     \item{\code{mov_med}}{Calculates the moving median for use as smoother}
#'     \item{\code{smoother_resid}}{Calculates the residuals of the smoother}
#'     \item{\code{bt_smoother}}{Bootstraps the errors of the autoregressive model fitted on the residuals of the smoother}
#'     \item{\code{conf_band}}{Calculates the confidence bounds which will be used to detect sustained change}
#'     \item{\code{detect_event}}{Detects sustained change using the confidence bounds
#' }
#' @docType package
#' @name edecob-package
NULL
#



#' Event DEtection using COnfidence Bounds
#'
#' Calculate a smoother of longitudinal data of the same measure and bootstrap the errors of the autoregressive
#' model fitted on the smoother to form simultaneous
#' confidence bounds of a certain level (mathematical details below).
#' Define an event if the simultaneous confidence bound is
#' within a chosen interval for a predefined amount of time. When data from
#' multiple sources is provided, the calculation will be done separately for
#' each source.
#'
#'
#' @param data A data frame in long format containing the data for which events
#'   is to be detected. This means that each measurement corresponds to a row
#'   and the columns are (in order): source (the device or person from which the
#'   data was collected), point in time,
#'   measurement value, lower detection bound, and upper detection bound.
#'
#'   The source is expected to
#'   be a string; the time point, measurements, and detection bounds are expected to be numerical.
#'   The detection bounds are in absolute value in the same unit as the
#'   values and each is expected to be identical for the same source. For details see below.
#' @param smoother Which smoother is to be used. Use \code{mov_med} for the
#'   moving median. When using the moving median, the parameter \code{width} must
#'   be given to specify the size of the window over which the moving median is
#'   to be taken. Defaults to the moving median.
#' @param bt_tot_rep The number of iterations for the bootstrap computation. Because of
#'   run time, it is recommended to keep this number below 500. Defaults to 100.
#' @param min_change_dur The minimal number of days that the confidence bounds
#'   need to stay inside the detection bounds in order for an event to be
#'   detected. Defaults to 84, i.e. 12 weeks.
#' @param conf_band_lvl The confidence level for the simultaneous confidence
#'   bands. Defaults to 0.95.
#' @param time_unit A string containing the unit of time used, in singular form.
#'   Defaults to day.
#' @param ... Additional parameters to be given to the function. Possible
#'   parameters for the model are \code{order} and \code{min_pts_in_win}. For
#'   the moving median, a \code{width} is required.
#'
#'   The parameter \code{min_pts_in_win}
#'   defines the minimal number of
#'   measurements required to be in the time window for the median to be calculated.
#'   Defaults to 1.
#'
#'   If the parameter \code{order} is given, that number will be the (maximal)
#'   order of the autoregressive model. If no \code{order} is given, it will be
#'   determined using the Akaike information criterion.
#'
#'   When the moving
#'   median is used as the smoother, \code{width} is expected. If no \code{width} is
#'   given, it will default to 84.
#'
#'
#' @details
#' In case detection is wanted for a one sided change (e.g. give an event if
#' the confidence bounds drop below a threshold) then the upper or lower detection
#' bound can be chosen to be Inf or -Inf respectively.
#'
#' For the moving median, the width is the total size of the window, meaning
#' that for the value corresponding to day x, the data points from day
#' x - \code{width}/2 to x + \code{width}/2 will be used for the calculation
#' of the median.
#'
#' If there is no data for two times \code{width} consecutive time units, there
#' will be time points at which no confidence bound can be calculated. In this
#' case, it will be assumed that the confidence bound is outside of the
#' detection interval when detecting sustained change.
#'
#' In case there are multiple instances where the algorithm would detect a
#' sustained change (i.e. if after the first sustained change the confidence
#' bounds leave the detection interval and then return into it
#' for longer than \code{min_change_dur} time units) then only the first
#' sustained change would be detected.
#'
#' Please note that the event onset could be on a date where there are no actual
#' measurements. This can happen when there is a gap in the data. In this case, the
#' confidence bounds will extend into the gap.
#' If the confidence bounds in this period are outside the detection interval and
#' remain outside for the next \code{min_change_duration} time units,
#' the event onset will be in this gap.
#'
#' The censoring date is based on the last date where the confidence bounds can be
#' calculated. We do not extend the confidence bounds to the last data point so
#' that the confidence bounds don't change in case we obtain new measurements
#' with time points later than the latest time point at which we have a measurement.
#'
#'
#' @return If \code{data} contains only a single source, the function returns
#'   a list of 13 variables: \describe{
#'   \item{\code{event}}{gives a list with four values: \code{event_detected},
#'     \code{event_onset}, \code{event_duration}, and \code{event_stop}.
#'     \describe{\item{\code{event_detected}}{gives
#'     whether an event was detected}
#'     \item{\code{event_onset}}{gives the first time point at which the upper or lower bound
#'     of the confidence band is inside the detection bounds, and after which it
#'     stays inside the detection bounds for at least \code{min_change_dur}
#'     consecutive time units} \item{\code{event_duration}}{gives the number of time units the upper or lower bound
#'     of the confidence band stays inside the detection bounds
#'     after \code{event_onset}} \item{\code{event_stop}}{gives whether the confidence
#'     bounds stay inside the detection bounds until
#'     the last time point at which we can calculate the confidence bound or not.}}
#'     }
#'   \item{\code{conf_band}}{gives a data frame containing the confidence bands.
#'     The columns are source, time point, lower bound, and upper
#'     bound of the confidence band.}
#'   \item{\code{smoother_pts}}{gives a data frame containing the smoother.
#'     The columns are source, time point, and the smoother}
#'   \item{\code{data}}{gives the data but with four additional columns:
#'     \code{event_detected}, \code{event_onset}, \code{event_duration}, and
#'    \code{event_stop}. They contain the same values as in \code{event}.}
#'   \item{\code{detec_lower}}{gives the lower detection bound.}
#'   \item{\code{detec_upper}}{gives the upper detection bound.}
#'   \item{\code{smoother}}{gives the smoother used.}
#'   \item{\code{min_change_dur}}{gives the smallest consecutive number of time units
#'     the confidence bounds must stay within the detection bounds in order for an event to be detected.}
#'   \item{\code{conf_band_lvl}}{gives the level of the simultaneous confidence band.}
#'   \item{\code{bt_tot_rep}}{gives the total amount of bootstrap repetitions performed.}
#'   \item{\code{call}}{gives the function call.}
#'   \item{\code{col_names}}{gives the original column names of the data.}
#'   \item{\code{time_unit}}{gives the unit of time used.}}
#'   If \code{data} contains more than one source, the output will be a list
#'   with one more element than the number of sources in \code{data}. Every
#'   element in this list will again be a list named after the corresponding sources
#'   with 13 items as described above except for the last one. The last element
#'   in the list is called \code{event_info} and is a data frame containing the
#'   information from \code{event} from each patient. \code{event_info} will thus
#'   have the following columns: \code{source}, \code{event_detected},
#'   \code{event_onset}, \code{event_duration}, and \code{event_stop}.
#'
#'
#' @section Mathematical background:
#' The mathematical background will be explained in the following sections.
#' @section Moving Median:
#'   Consider a sample \eqn{X₁,\dots, Xₙ} of size \eqn{n} and the
#'   reordering \eqn{X₍₁₎,\dots, X₍ₙ₎} such
#'   that \eqn{X₍₁₎ \le X₍₂₎ \le \dots \le X₍ₙ₎}, commonly
#'   called the order statistic. Then for \eqn{n} even the median usually
#'   defined as \deqn{median(X₁,\dots, Xₙ) = X₍ₖ₎, where k = n/2.} In the
#'   case where \eqn{n} is odd the median is
#'   defined as \deqn{median(X₁,\dots, Xₙ) = 1/2(X₍ₖ₎ + X₍ₖ₊₁₎), where k = n/2.} Let the
#'   time points at which the measurements \eqn{X₁, \dots, Xₙ} were taken
#'   be \eqn{t₁, \dots, tₙ}.
#'   Let \eqn{T} a fixed positive amount of time. Then
#'   \deqn{S(t) = median({Xⱼ | t - T/2 \le tⱼ \le t + T/2})}
#'   is defined as the moving median at time point \eqn{t} with window size \eqn{T}.
#' @section The Model:
#' An autoregressive (AR) model is used to model the residuals of the smoother \eqn{\eta}:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sumᵖⱼ₌₁ \phiⱼ \eta(t - j) + \epsilon} where
#' variable \eqn{t} is the time point, \eqn{Y(t)} the data point at time point \eqn{t},
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the difference between the smoother
#' and the measurement at time point \eqn{t}, \eqn{p}
#' the order of the AR model, \eqn{\phiⱼ} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model. The order is calculated using the
#' Akaike information criterion (AIC) if it was not given in the function call.
#' @section Bootstrap:
#' In the following, the star * denotes a bootstrapped value. The bootstrap
#' procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(tᵢ) = Y(tᵢ) - S(tᵢ)}.
#'   \item Fit an AR(p) model to \eqn{\eta(tᵢ)} to obtain the coefficients
#'     \eqn{\phi₁,\dots, \phiₚ} and \eqn{\epsilon(tᵢ) = \eta(tᵢ) -
#'     \sumᵖⱼ₌₁ \phiⱼ \eta(tᵢ - tᵢ₋ⱼ)} the error of the AR model.
#'   \item Randomly choose a \eqn{\epsilon(tᵢ)*} with replacement from \eqn{\epsilon(tₚ₊₁),\dots,
#'     \epsilon(tₙ)} to obtain \deqn{Y(tᵢ)* = S(tᵢ) + \eta(tᵢ)*,} where  \deqn{\eta(tᵢ)* = \sumᵖⱼ₌₁ \phiⱼ \eta(tᵢ₋ⱼ)*+ \epsilon(tᵢ₋ⱼ)*} the bootstrapped residuals of the smoother.
#'   \item Compute \eqn{S(.)* = g(Y(t₁),\dots, Y(tₙ))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times to obtain \eqn{S(tᵢ)*ᵦ} for \eqn{β = 1,\dots,}
#'   \code{bt_tot_rep}.
#' }
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
#'
#' @seealso \code{\link{summary.edecob}}, \code{\link{plot.edecob}}
#'
#'
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
#' # Let us examine the example_data dataset
#' head(example_data)
#' example_event <- edecob(example_data, width = 49)
#' names(example_event)
#'
#' # example_event contains the event data for each source
#' plot(example_event$`Subject 1`)
#'
#' # and a data frame containing the event information for all patients
#' example_event$event_info
#'
#' # Using this data frame, we can draw a survival plot
#' library("survival")
#' plot(survfit(Surv(time = event_onset, event = event_detected) ~ 1,
#'              data = example_event$event_info),
#'      conf.int = FALSE, xlim = c(0,350), ylim = c(0,1), mark.time = TRUE)
#'
#'
#'
edecob <- function(data,
                   smoother = "mov_med",
                   min_change_dur = 12*7,
                   conf_band_lvl = 0.95,
                   bt_tot_rep = 100,
                   time_unit = "day",
                   ...) {

  data_raw <- data

  if (!("col_names" %in% names(match.call()))) {
    col_names <- colnames(data)
    colnames(data) <- c("source", "time_point", "value", "detec_lower", "detec_upper")
  } else {
    col_names <- list(...)$col_names
  }

  data <- data.frame("source" = unlist(data$source),
                     "time_point" = unlist(data$time_point),
                     "value" = unlist(data$value),
                     "detec_lower" = unlist(data$detec_lower),
                     "detec_upper" = unlist(data$detec_upper))

  stopifnot(
    "Data not a data frame" = is.data.frame(data),
    "Data empty" = nrow(data) > 0,
    # "source not character" = is.character(data[,1]),
    "Time points not numeric" = is.numeric(data[,2]),
    "Measurements not numeric" = is.numeric(data[,3]),
    "Upper bound of detection interval not numeric" = is.numeric(data[,5]),
    "Lower bound of detection interval not numeric" = is.numeric(data[,4]),
    "Upper bound of detection interval not all equal for at least one source" = {
      all(do.call(c, lapply(unique(data$source), function(x){
        return(length(unique(data$detec_upper[data$source == x])) == 1)
      })))},
    "Lower bound of detection interval not all equal for at least one source" = {
      all(do.call(c, lapply(unique(data$source), function(x){
        return(length(unique(data$detec_lower[data$source == x])) == 1)
      })))},
    "Upper bound of detection interval contains NA values" = sum(is.na(data[,5])) == 0,
    "Lower bound of detection interval contains NA values" = sum(is.na(data[,4])) == 0,
    "Lower bound of detection interval is larger than upper bound of detection interval for at least one source" = sum(data[,4] > data[,5]) == 0,
    "Lower bound of detection interval is equal to upper bound of detection interval for at least one source" = sum(data[,4] == data[,5]) == 0
  )



  if (sum(is.na(data$value)) > 1) {
    warning("Removing rows where value is NA", immediate. = TRUE)
    data <- data[!is.na(data$value), ]
  }
  if (sum(is.na(data$time_point)) > 1) {
    warning("Removing rows where time point is NA", immediate. = TRUE)
    data <- data[!is.na(data$time_point), ]
  }



  # making sure that width is not repeated in future function calls
  no_width_given <- FALSE
  if (smoother == "mov_med" && !("width" %in% names(list(...)))) {
    warning("Parameter width not given after choosing the moving median as smoother. Defaulting to 84.", immediate. = TRUE)
    width <- 84
    no_width_given <- TRUE

  } else if (smoother == "mov_med") {
    width <- list(...)$width
  }

  # multiple patients
  if (length(unique(data$source)) > 1) {
    patients_event_data <-
      lapply(split(data, factor(data$source)), edecob,
             smoother, min_change_dur,conf_band_lvl, bt_tot_rep, time_unit,
             "col_names" = col_names, ...)
    patients_event_data$event_info <- as.data.frame(do.call(rbind,
      lapply(patients_event_data, function(x) {
        return(list(x$event$source, x$event$event_detected, x$event$event_onset, x$event$event_duration, x$event$event_stop))
      })))
    colnames(patients_event_data$event_info) <-
      c("source", "event_detected", "event_onset", "event_duration", "event_stop")
    patients_event_data$event_info$source <- unlist(patients_event_data$event_info$source)
    patients_event_data$event_info$event_detected <- unlist(patients_event_data$event_info$event_detected)
    patients_event_data$event_info$event_onset <- unlist(patients_event_data$event_info$event_onset)
    patients_event_data$event_info$event_duration  <- unlist(patients_event_data$event_info$event_duration)
    patients_event_data$event_info$event_stop  <- unlist(patients_event_data$event_info$event_stop)
    return(patients_event_data)
  }

  print(data[1,1])

  # calculate the smoother
  if (smoother == "mov_med") {
    if ("min_pts_in_win" %in% names(match.call)) {
      smoother_pts <- mov_med(data, width, list(...)$min_pts_in_win)
    }
    smoother_pts <- mov_med(data, width)
  } else {
    warning("Smoother not recognized. Defaulting to moving median.")

    if (!("width" %in% names(match.call()))) {
      warning("Parameter width not given after choosing the moving median as smoother. Defaulting to 84.")
      width <- 84
      no_width_given <- TRUE
    }
    smoother_pts <- mov_med(data, width)
  }

  # calculate residuals of the smoother
  smoother_resid <- smoother_resid(data, smoother_pts)

  # bootstrap the errors of the AR model fitted on the residuals
  if (no_width_given) {
    bt_smoother <- bt_smoother(data, smoother, smoother_pts, smoother_resid, bt_tot_rep, "width" = width, ...)
  } else {
    bt_smoother <- bt_smoother(data, smoother, smoother_pts, smoother_resid, bt_tot_rep, ...)
  }

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)

  # detect events using confidence bands
  event <- detect_event(conf_band, data$detec_upper[1], data$detec_lower[1], min_change_dur)

  # add columns with event information to data
  data_raw$event <- event$event_detected
  data_raw$event_onset <- event$event_onset
  data_raw$event_duration <- event$event_duration
  data_raw$event_stop <- event$event_stop

  # compose output
  output <- list(
    "source" = data$source[1],
    "event" = event,
    "conf_band" = conf_band,
    "smoother_pts" = smoother_pts,
    "data" = data_raw,
    "smoother" = smoother,
    "detec_lower" = data$detec_lower[1],
    "detec_upper" = data$detec_upper[1],
    "min_change_dur" = min_change_dur,
    "conf_band_lvl" = conf_band_lvl,
    "bt_tot_rep" = bt_tot_rep,
    "call" = match.call(),
    "col_names" = col_names,
    "time_unit" = time_unit
  )
  class(output) <- "edecob"

  return(output)
}


#' edecob: Event Detection Using Confidence Bounds
#'
#' Detect sustained change (or events) in time-dependent data while accounting for
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
#'
NULL



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
#'   data was collected), point in time, and
#'   measurement value. If custom detection bounds are chosen, the folloing two
#'   columns must be added: lower detection bound, and upper detection bound.
#'
#'   The source is expected to
#'   be a string; the time point are integers; measurements, and detection bounds are expected to be numerical.
#'   The detection bounds are in absolute value in the same unit as the
#'   values and each is expected to be identical for the same source.
#'
#'   In case detection is wanted for a one sided change (e.g. give an event if
#'   the confidence bounds drop below a threshold) then the upper or lower detection
#'   bound can be chosen to be Inf or -Inf respectively.
#' @param smoother A string specifying which smoother is to be used. Use \code{mov_med} for the
#'   moving median. When using the moving median, the parameter \code{med_win} must
#'   be given to specify the size of the window over which the moving median is
#'   to be taken. Defaults to the moving median.
#' @param resample_method A string that determines how to resample the errors of the
#'   autoregression for the bootstrap. Default is \code{all}, meaning that the epsilon of a
#'   certain time point are resampled from all time points. \code{past} only
#'   considers epsilon corresponding to a time point prior to the one being
#'   resampled. \code{window} resamples the epsilon from the window from which
#'   the moving median is taken.
#' @param bt_tot_rep The number of iterations for the bootstrap computation. Because of
#'   run time, it is recommended to keep this number below 500. Defaults to 100.
#' @param min_change_dur The minimal number of days that the confidence bounds
#'   need to stay inside the detection bounds in order for an event to be
#'   detected. Defaults to 84, i.e. 12 weeks.
#' @param conf_band_lvl The confidence level for the simultaneous confidence
#'   bands. Defaults to 0.95. When detection of events using only the smoother
#'   is desired, \code{conf_band_lvl} can be chosen to be 0.
#' @param time_unit A string containing the unit of time used, in singular form.
#'   Defaults to day.
#' @param detect A string specifying how the detection bounds are to be chosen.
#'   \code{below}, \code{above}, and \code{custom} can be chosen. \code{below}
#'   detects decreases in value, \code{above} detects increases in value, and
#'   \code{custom} can be used to manually add detection bounds for each subject.
#'   When \code{above} or \code{below} are used, the detection bound will be x percent
#'   above or below the median of the first y days, where x is \code{detect_factor}
#'   and y is \code{detect period}.
#' @param detect_factor A number specifying the factor by which the median of
#'   the fist \code{bline_period} days is to be multiplied to obtain the
#'   detection bounds. E.g. 0.9 sets the detection bound 10 percent below the
#'   said median.
#' @param bline_period A number specifying the number of time units from which
#'   data should be taken to calculate the median to obtain the detection bounds.
#' @param ... Additional parameters to be given to the function. Possible
#'   parameters for the model are \code{order} and \code{min_pts_in_win}. For
#'   the moving median, a \code{med_win} is required. When resampling from
#'   window, a \code{resample_win} may be given.
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
#'   median is used as the smoother, \code{med_win} is expected. If no \code{med_win} is
#'   given, it will default to \code{c(-42, 42)}.
#'
#'   When resampling from window, one can choose the window size for the
#'   resampling window with \code{resample_win} by giving a window like e.g. \code{c(-14,14)}..
#'
#'
#' @details
#'
#' For the moving median, the med_win is the total size of the window, meaning
#' that for the value corresponding to day x, the data points from day
#' x + \code{med_win[1]} to x + \code{med_win[2]} will be used for the calculation
#' of the median.
#'
#' If there is no data for two times \code{med_win[2]-med_win[1]} consecutive time units, there
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
#' The \code{edecob} function runs the functions \code{mov_med}, \code{smoother_resid},
#' \code{bt_smoother}, \code{conf_band}, and \code{detect_event} in this order
#' for all subjects given. If desired, the functions can also manually be
#' applied for the data to obtain e.g. the confidence bands. Note that in order
#' to run one of these functions, the output of the previous functions are needed.
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
#'   \item{\code{resample_method}}{gives the resampling method used for the bootstrap.}
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
#'   Consider a sample \eqn{X_1,\dots, X_n} of size \eqn{n} and the
#'   reordering \eqn{X_{(1)},\dots, X_{(n)}} such
#'   that \eqn{X_{(1)} \le X_{(2)} \le \dots \le X_{(n)}}, commonly
#'   called the order statistic. Then for \eqn{n} even the median usually
#'   defined as \deqn{median(X_1,\dots, X_n) = X_{(k)}, \mathrm{where} \; k = n/2.} In the
#'   case where \eqn{n} is odd the median is
#'   defined as \deqn{median(X_1,\dots, X_n) = 1/2(X_{(k)} + X_{(k+1)}), \mathrm{where} \; k = n/2.} Let the
#'   study days at which the measurements \eqn{X_1, \dots, X_n} were taken
#'   be \eqn{t_1, \dots, t_n}.
#'   Let \eqn{T} a fixed positive amount of time. Then the moving median at time
#'   point \eqn{t} with window size \eqn{T} is defined as
#'   \deqn{S(t) = median({X_j | t - T/2 \le t_j \le t + T/2}).}
#'
#' @section The Model:
#' An autoregressive (AR) model is used to model the residuals of the smoother \eqn{\eta}:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sum^p_{j = 1} \phi_j \eta(t - j) + \epsilon} where
#' variable \eqn{t} is the study day, \eqn{Y(t)} the data point at study day \eqn{t},
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the difference between the smoother
#' and the measurement at study day \eqn{t}, \eqn{p}
#' the order of the AR model, \eqn{\phi_j} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model. The order is calculated using the
#' Akaike information criterion (AIC) if it was not given in the function call.
#'
#' @section Bootstrap:
#' In the following, the star * denotes a bootstrapped value. The bootstrap
#' procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(t_i) = Y(t_i) - S(t_i)}.
#'   \item Fit an AR(p) model to \eqn{\eta(t_i)} to obtain the coefficients
#'     \eqn{\phi_1,\dots, \phi_p} and \eqn{\epsilon(t_i) = \eta(t_i) -
#'     \sum^p_{j = 1} \phi_j \eta(t_i - t_{i-j})} the error of the AR model.
#'   \item Randomly choose a \eqn{\epsilon(t_i)^*} with replacement from \eqn{\epsilon(t_{p+1}),\dots,
#'     \epsilon(t_n)} to obtain \deqn{Y(t_i)^* = S(t_i) + \eta(t_i)^*,} where  \deqn{\eta(t_i)^* = \sum^p_{j = 1} \phi_j \eta(t_{i-j})^*+ \epsilon(t_{i-j})^*} the bootstrapped residuals of the smoother.
#'   \item Compute \eqn{S(.)^* = g(Y(t_1),\dots, Y(t_n))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times to obtain \eqn{S(t_i)^*_b} for \eqn{\beta = 1,\dots,}
#'   \code{bt_tot_rep}.
#' }
#'
#' @section Calculation of the Confidence Bounds:
#' The confidence bounds are calculated as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ q_x(t_i), q_{1-x}(t_i) i = 1,\dots, N}
#'     where \deqn{q_x(t_i) = inf\left\{u; P^*[S(t_i)^*_b - S(t_i) \le u] \ge x\right\} } is a
#'     pointwise bootstrap quantile, \eqn{S(t_i)^*_b} the bootstrapped smoother,
#'     and \eqn{N} the number of measurements or rows in \code{data}, in our case the number of rows.
#'   \item We vary the pointwise error \eqn{x} until \deqn{P^*[q_x(t_i) \le S(t_i)^*_b - S(t_i) \le q_{1-x}(t_i) \forall i = 1,\dots, N] \approx 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[q_x(t_i), q_{1-x}(t_i)]} is approximately \eqn{1-\alpha}.
#'   \item We define
#'   \deqn{ I_n(t_i) = [S(t_i) +  q_x(t_i), S(t_i) + q_{1-x}(t_i)] \forall i = 1,\dots, N}
#'   the confidence bounds. Then \eqn{{I_n(t_i); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
#'
#'}
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
#' library(edecob)
#'
#' # Let us examine the example_data dataset
#' head(example_data, 3)
#' #>     subject study_day jump_height detect_lower detect_upper
#' #> 1 Subject 1         1    55.60844         -Inf     54.41227
#' #> 2 Subject 1         4    57.77688         -Inf     54.41227
#' #> 3 Subject 1         7    57.59584         -Inf     54.41227
#'
#' # We apply the main fuction of the package onto our example_data
#' example_event <- edecob(example_data, med_win = c(-21,21), bt_tot_rep = 20,
#'                         min_change_dur = 50)
#' #> Warning in edecob(example_data, med_win = c(-21, 21), bt_tot_rep = 20,
#' #> min_change_dur = 50) :
#' #>   Removing rows where value is NA
#' names(example_event)
#' #> [1] "Subject 1"  "Subject 2"  "Subject 3"  "event_info"
#'
#' # example_event contains the event data for each source
#' plot(example_event$`Subject 1`)
#' plot(example_event$`Subject 2`)
#' plot(example_event$`Subject 3`)
#'
#' # example_event also contains a data frame containing the event information for all patients
#' example_event$event_info
#' #>           event_detected event_onset event_duration event_stop
#' #> Subject 1           TRUE         169             87       TRUE
#' #> Subject 2           TRUE         205             51       TRUE
#' #> Subject 3          FALSE         306             38      FALSE
#'
#' # Using this data frame, we can draw a survival plot
#' library("survival")
#' plot(survfit(Surv(time = event_onset, event = event_detected) ~ 1,
#'              data = example_event$event_info),
#'      conf.int = FALSE, xlim = c(0,350), ylim = c(0,1), mark.time = TRUE,
#'      xlab = "Study Day", ylab = "Survival Probability", main = "Survival plot")
edecob <- function(data,
                   smoother = "mov_med",
                   resample_method = "all",
                   min_change_dur = 84,
                   conf_band_lvl = 0.95,
                   bt_tot_rep = 100,
                   time_unit = "day",
                   detect = "below",
                   detect_factor = 1,
                   bline_period = 14,
                   ...) {

  data_raw <- data

  if (!("col_names" %in% names(match.call()))) {
    col_names <- colnames(data)
    if (detect == "custom"){
      colnames(data) <- c("source", "time_point", "value", "detec_lower", "detec_upper")
    } else {
      colnames(data) <- c("source", "time_point", "value")
    }

  } else {
    col_names <- list(...)$col_names
  }

  if (detect == "custom"){
    data <- data.frame("source" = unlist(data$source),
                       "time_point" = unlist(data$time_point),
                       "value" = unlist(data$value),
                       "detec_lower" = unlist(data$detec_lower),
                       "detec_upper" = unlist(data$detec_upper))
  } else {
    data <- data.frame("source" = unlist(data$source),
                       "time_point" = unlist(data$time_point),
                       "value" = unlist(data$value))
  }

  stopifnot(
    "Data not a data frame" = is.data.frame(data),
    "Data empty" = nrow(data) > 0,
    "Time points not numeric" = is.numeric(data[,2]),
    "Measurements not numeric" = is.numeric(data[,3])
  )

  if (detect == "custom") {
    stopifnot(
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
  }

  if ("med_win" %in% names(list(...))) {
    stopifnot("Window of the moving median does not contain two numbers" = (length(list(...)$med_win) == 2),
               "Lower bound of the window for the moving median is not smaller than the upper bound" = (list(...)$med_win[1] <= list(...)$med_win[2]))
  }


  if (sum(is.na(data$value)) > 1) {
    warning("Removing rows where value is NA", immediate. = TRUE)
    data <- data[!is.na(data$value), ]
  }
  if (sum(is.na(data$time_point)) > 1) {
    warning("Removing rows where time point is NA", immediate. = TRUE)
    data <- data[!is.na(data$time_point), ]
  }
  if (!all(data[,2] == floor(data[,2]))) {
    warning("Time points not integer. Rounding time points down to closest integer.", immediate. = TRUE)
    data[,2] <- floor(data[,2])
  }



  # making sure that width is not repeated in future function calls
  no_width_given <- FALSE
  if (!("med_win" %in% names(list(...)))) { # (smoother == "mov_med" || smoother == "mov_mean") &&
    warning("Parameter med_win not given after choosing the moving median as smoother. Defaulting to c(-42,42)", immediate. = TRUE)
    med_win <- c(-42, 42)
    no_width_given <- TRUE

  } else if (smoother == "mov_med") {
    med_win <- list(...)$med_win
  }

  # multiple patients
  if (length(unique(data$source)) > 1) {
    patients_event_data <-
      lapply(split(data, factor(data$source)), edecob,
             smoother, resample_method, min_change_dur,conf_band_lvl,
             bt_tot_rep, time_unit, detect, detect_factor, bline_period,
             "col_names" = col_names, "med_win" = med_win, ...)
    patients_event_data$event_info <- as.data.frame(do.call(rbind,
      lapply(patients_event_data, function(x) {
        return(list(x$event$source, x$event$event_detected, x$event$event_onset, x$event$event_duration, x$event$event_stop))
      })))
    patients_event_data$event_info$source <- row.names(patients_event_data$event_info)
    patients_event_data$event_info <- patients_event_data$event_info[,c(6, 2:5)]
    colnames(patients_event_data$event_info) <-
      c("source", "event_detected", "event_onset", "event_duration", "event_stop")
    patients_event_data$event_info$source <- unlist(patients_event_data$event_info$source)
    patients_event_data$event_info$event_detected <- unlist(patients_event_data$event_info$event_detected)
    patients_event_data$event_info$event_onset <- unlist(patients_event_data$event_info$event_onset)
    patients_event_data$event_info$event_duration  <- unlist(patients_event_data$event_info$event_duration)
    patients_event_data$event_info$event_stop  <- unlist(patients_event_data$event_info$event_stop)
    return(patients_event_data)
  }

  data <- data[order(data$time_point), ]

  # calculate the upper and lower detection bounds if needed
  if (detect == "custom") {
    detection_bound_lower <- data$detec_lower[1]
    detection_bound_upper <- data$detec_upper[1]
  } else if (detect == "below") {
    detection_bound_lower <- -Inf
    bound_data <- data$value[data$time_point <= data$time_point[1] + bline_period]
    detection_bound_upper <- median(bound_data)*detect_factor
  } else if (detect == "above") {
    detection_bound_upper <- Inf
    bound_data <- data$value[data$time_point <= data$time_point[1] + bline_period]
    detection_bound_lower <- median(bound_data)*detect_factor
  }


  # calculate the smoother
  if (smoother == "mov_med") {
    if ("min_pts_in_win" %in% names(match.call)) {
      smoother_pts <- mov_med(data, med_win, list(...)$min_pts_in_win)
    }
    smoother_pts <- mov_med(data, med_win)
  } else if (smoother == "mov_mean") {
    if ("min_pts_in_win" %in% names(match.call)) {
      smoother_pts <- mov_mean(data, med_win, list(...)$min_pts_in_win)
    }
    smoother_pts <- mov_mean(data, med_win)
  } else {
    warning("Smoother not recognized. Defaulting to moving median.")

    if (!("med_win" %in% names(match.call()))) {
      warning("Parameter med_win not given after choosing the moving median as smoother. Defaulting to c(-42,42).")
      med_win <- c(-42,42)
      no_width_given <- TRUE
    }
    smoother_pts <- mov_med(data, med_win)
  }

  # calculate residuals of the smoother
  smoother_resid <- smoother_resid(data, smoother_pts)

  # bootstrap the errors of the AR model fitted on the residuals
  if (no_width_given) {
    bt_smoother <- bt_smoother(data, smoother, resample_method, smoother_pts, smoother_resid, bt_tot_rep, "med_win" = med_win, ...)
  } else {
    bt_smoother <- bt_smoother(data, smoother, resample_method, smoother_pts, smoother_resid, bt_tot_rep, ...)
  }

  # calculate the confidence bands
  conf_band <- conf_band(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)

  # detect events using confidence bands
  event <- detect_event(conf_band, detection_bound_lower, detection_bound_upper, min_change_dur)

  # add columns with event information to data
  colnames(data_raw) <- colnames(data) #<- c("source", "time_point", "value", "detec_lower", "detec_upper")
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
    "resample_method" = resample_method,
    "detec_lower" = detection_bound_lower,
    "detec_upper" = detection_bound_upper,
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

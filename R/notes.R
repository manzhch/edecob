# legacy ----

# edecob
#' @param thresh_diff The absolute or percentage (depending on \code{thresh_method})
#'   difference between the baseline and
#'   threshold. A negative number indicates a threshold below the baseline and
#'   similarly a positive number indicates a threshold above the baseline. If it
#'   is zero, the algorithm will detect events below the baseline.
#'   Defaults to -20%. To detect events above the baseline with threshold
#'   equal to baseline, use a very small number like 1e-10.
#' @param thresh_method Whether the threshold is a percentage or absolute
#'   difference from the baseline. Use \code{percent} or \code{absolute}.
#'   Defaults to percent.
#' @param basel_start The study day of first data point to be used for the
#'   baseline. Defaults to the smallest available study day.
#' @param basel_end The study day of last data point to be used for the
#'   baseline. Defaults to six weeks after \code{basel_start}.
#' @param basel_method How the baseline should be calculated. Can take any of
#'   the following values: \code{median, mean, quantile}. When using quantile,
#'   an additional argument \code{quantile} specifying the quantile needs to be
#'   given. Defaults to median.
#' @param ... Additional parameters to be given to the function. When the moving
#'   median is used, a \code{width} parameter is expected. When using the quantile
#'   for baseline method, a \code{quantile} parameter is expected.
#'


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
# data_non_learn <- data[data$study_day >= basel_start, ]

# calculate baseline
# if (basel_method == "median") {
#   basel <- stats::median(data$value[as.logical(
#     (data$study_day >= basel_start) * (data$study_day <= basel_end))])
# } else if (basel_method == "mean") {
#   basel <- mean(data$value[as.logical(
#     (data$study_day >= basel_start) * (data$study_day <= basel_end))])
# } else if (basel_method == "quantile") {
#   basel <- stats::quantile(data$value[as.logical(
#     (data$study_day >= basel_start) * (data$study_day <= basel_end))], quantile)
# } else {
#   warning("Baseline method not recognized. Defaulting to median.")
#   basel <- stats::median(data$value[as.logical(
#     (data$study_day >= basel_start) * (data$study_day <= basel_end))])
# }


# calculate threshold
# if (thresh_method == "percent") {
#   thresh <- basel * (1 + thresh_diff/100)
# } else if (thresh_method == "absolute") {
#   thresh <- basel + thresh_diff
# } else {
#   warning("Threshold method not recognized. Defaulting to percent.")
#   thresh <- basel * (1 + thresh_diff/100)
# }

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



# plot
#' @param ... If data from more than one subject is included in
#'   \code{event_data}, specify \code{subj_id}. Defaults to the first subject
#'   found in \code{event_data$data}.


# summary
# cat("Baseline of", object$basel, "calculated using the", object$basel_method,
#     "using data from study day", object$basel_start,
#     "to", object$basel_end, "\n")
# if (object$thresh_method == "percentage") {
#   cat("Threshold of", object$thresh, "is", object$thresh*10, "%", "difference from baseline", "\n")
# } else if (object$thresh_method == "absolute") {
#   cat("Threshold of", object$thresh, "is", object$thresh, object$thresh_method, "difference from baseline", "\n")
# }



# doc edecob --------



#' Event DEtection using COnfidence Bounds
#'
#' Calculate a smoother of the data and bootstrap the errors of the autoregressive
#' model fitted on the smoother to form simultaneous
#' confidence bounds. Define an event if the simultaneous confidence bound is
#' below or above a threshold for a predefined amount of time.
#'
#' @section Moving Median:
#'   Consider a sample \eqn{X_1,\dots, X_n} of size \eqn{n} and the
#'   reordering \eqn{X_{(1)},\dots, X_{(n)}} such
#'   that \eqn{X_{(1)} \le X_{(2)} \le \dots \le X_{(n)}}, commonly
#'   called the order statistic. Then for \eqn{n} even the median usually
#'   defined as \deqn{median(X_1,\dots, X_n) = X_{(k)}, where k = n/2.} In the
#'   case where \eqn{n} is odd the median is
#'   defined as \deqn{median(X_1,\dots, X_n) = 1/2(X_{(k)} + X_{(k+1)}), where k = n/2.} Let the
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
#'   \item Randomly choose a \eqn{\epsilon(t_i)*} with replacement from \eqn{\epsilon(t_{p+1}),\dots,
#'     \epsilon(t_n)} to obtain \deqn{Y(t_i)* = S(t_i) + \eta(t_i)*,} where  \deqn{\eta(t_i)* = \sum^p_{j = 1} \phi_j \eta(t_{i-j})*+ \epsilon(t_{i-j})*} the bootstrapped residuals of the smoother.
#'   \item Compute \eqn{S(.)* = g(Y(t_1),\dots, Y(t_n))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times to obtain \eqn{S(t_i)*_b} for \eqn{\beta = 1,\dots,}
#'   \code{bt_tot_rep}.
#' }
#'
#' @section Calculation of the Confidence Bounds:
#' The confidence bounds are calculated as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ q_x(t_i), q_{1-x}(t_i) i = 1,\dots, N}
#'
#'     where \deqn{q_x(t_i) = inf{u; P*[S(t_i)*_b - S(t_i) \le u] \ge x} } is a
#'     pointwise bootstrap quantile, \eqn{S(t_i)*_b} the bootstrapped smoother,
#'     and \eqn{N} the number of measurements or rows in \code{data}, in our case the number of rows.
#'   \item We vary the pointwise error \eqn{x} until \deqn{P*[q_x(t_i) \le S(t_i)*_b - S(t_i) \le q_{1-x}(t_i) \forall i = 1,\dots, N] \approx 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[q_x(t_i), q_{1-x}(t_i)]} is approximately \eqn{1-\alpha}.
#'   \item We define
#'   \deqn{ I_n(t_i) = [S(t_i) +  q_x(t_i), S(t_i) + q_{1-x}(t_i)] \forall i = 1,\dots, N}
#'   the confidence bounds. Then \eqn{{I_n(t_i); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
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
#' If the parameter \code{order} it not given, the order of the autoregressive
#' model will be determined using the Akaike information criterion (AIC).
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
#'   median is used, a \code{width} parameter is expected. If no \code{width} is
#'   given, it will default to 84. If the parameter
#'   \code{order} is given, that number will be the (maximal) order of the autoregressive
#'   model.
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
#' # We look at the level of Lake Huron
#' LakeHuron_event <-
#'   edecob(data.frame(subject = "LakeHuron",
#'                     year = tsp(LakeHuron)[1]:tsp(LakeHuron)[2],
#'                     data = LakeHuron),
#'          baseline = 580.5, threshold = 579.7, width = 10, min_change_dur = 20)
#' summary(LakeHuron_event)
#' # Notice in the plot that the event onset does not happen the first time the
#' # confidence bands drop below the threshold as the change is not sutained for long enough
#' plot(LakeHuron_event)
#'
#' # Let us see what happens when we introduce a gap into the data
#' LakeHuron_event2 <-
#'   edecob(data.frame(subject = "LakeHuron",
#'                     year = tsp(LakeHuron)[1]:tsp(LakeHuron)[2],
#'                     data = c(LakeHuron[1:25], rep(NA, 20), LakeHuron[46:98])),
#'          baseline = 580.5, threshold = 579.7, width = 10, min_change_dur = 10)
#' summary(LakeHuron_event2)
#' # Notice that the confidence bounds become wider around the gap
#' plot(LakeHuron_event2)


# doc conf_band


#' Confidence Bounds of the Smoother
#'
#' Calculate the confidence bounds of the smoother function using the bootstrap.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ q_x(t_i), q_{1-x}(t_i) i = 1,\dots, N}
#'
#'     where \deqn{q_x(t_i) = inf {u; P*[S(t_i)*_b - S(t_i) \le u] \ge x} } is a
#'     pointwise bootstrap quantile, \eqn{S(t_i)*_b the bootstrapped smoother},
#'     and \eqn{N} the number of measurements.
#'   \item We vary the pointwise error \eqn{2x} until \deqn{P*[q_x(t_i) \le S(t_i)*_b - S(t_i) \le q_{1-x}(t_i) \forall i = 1,\dots, N] \approx 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[q_x(t_i), q_{1-x}(t_i)]} is approximately \eqn{1-\alpha}.
#'   \item We let
#'   \deqn{ I_n(t_i) = [S(t_i) +  q_x(t_i), S(t_i) + q_{1-x}(t_i)] \forall i = 1, \dots, N}
#'   confidence bounds. Then \eqn{{I_n(t_i); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
#'
#'}
#'
#'
#' @inheritParams edecob
#' @param smoother_pts A data frame containing the smoother. Preferably the
#'   output of one of the smoother functions included in this package.
#' @param bt_smoother A data frame containing the bootstrapped smoother. Use the output of \code{bt_smoother}.
#'
#' @return A data frame containing the upper confidence bound, the lower confidence bound,
#'   and the study day corresponding to the bounds.
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'

# doc bt_smoother -------

#' Bootstrap the Smoother
#'
#' Fit an autoregressive model on the residuals of the smoother and bootstrap
#' the error of the autoregressive model. Reconstruct the smoother by adding the
#' bootstrapped error, the fitted autoregressive model, and the smoother to form
#' the bootstrapped smoother.
#'
#' An autoregressive (AR) model is used for the residuals of the smoother:
#' \deqn{Y(t) = S(t) + \eta(t)}
#' \deqn{\eta(t) = \sum^{p}_{j = 1} \phi_j \eta(t - j) + \epsilon}
#' where \eqn{t} is the point in time, \eqn{Y(t)} the data point,
#' \eqn{S(t)} a smoother, \eqn{\eta(t)} the residual of the smoother, \eqn{p}
#' the order of the AR model, \eqn{\phi_j} the coefficients of the AR model, and
#' \eqn{\epsilon} the error of the AR model.
#'
#' The bootstrap procedure is as follows:
#' \enumerate{
#'   \item Compute the smoother \eqn{S(t)}.
#'   \item Compute the residuals \eqn{\eta(t_i) = Y(t_i) - S(t_i)}.
#'   \item Fit an AR(p) model to \eqn{\eta(t_i)} to obtain the coefficients
#'     \eqn{\phi_1, \dots, \phi_p} and residuals \eqn{\epsilon(t_i) = \eta(t_i) -
#'     \sum^{p}_{j = 1} \phi_j \eta(t_i - t_{i-j})}.
#'   \item Resample \eqn{\epsilon(t_i)*} from \eqn{\epsilon(t_{p+1}), \dots,
#'     \epsilon(t_n)} to obtain \deqn{Y(t_i)* = S(t_i) + \eta(t_i)*,} where  \deqn{\eta(t_i)* = \sum^{p}_{j=1} \phi_j \eta(t_{i-j})*+ \epsilon(t_{i-j})*.}
#'   \item Compute \eqn{S(.)* = g(Y(t_1), \dots, Y(t_n))} where \eqn{g} is the
#'     function with which the smoother is calculated.
#'   \item Repeat steps 4 and 5 \code{bt_tot_rep} times.
#' }
#'
#'
#' @param smoother_pts A data frame containing the smoother. Preferably the
#'   output of one of the smoother functions included in this package.
#' @param resid A vector containing the residuals of the smoother to the data
#'   points.
#' @inheritParams edecob
#'
#' @return A data frame containing the bootstrap repetitions of the smoother.
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'




# doc mov_med -----


#' Moving Median over a Time Window
#'
#' Calculates the moving median over a time window around a time point for the
#' all time points between the first and last study day provided.
#'
#' Consider a sample \eqn{X_1,\dots, X_n} of size \eqn{n} and the
#' reordering \eqn{X_{(1)},\dots, X_{(n)}} such
#' that \eqn{X_{(1)} \le X_{(2)} \le \dots \le X_{(n)}}, commonly
#' called the order statistic. Then for \eqn{n} even the median usually
#' defined as \deqn{median(X_1,\dots, X_n) = X_{(k)}, where k = n/2.} In the
#' case where \eqn{n} is odd the median is
#' defined as \deqn{median(X_1,\dots, X_n) = 1/2(X_{(k)} + X_{(k+1)}), where k = n/2.} Let the
#' study days at which the measurements \eqn{X_1, \dots, X_n} were taken
#' be \eqn{t_1, \dots, t_n}.
#' Let \eqn{T} a fixed positive amount of time. Then the moving median at time
#' point \eqn{t} with window size \eqn{T} is defined as
#' \deqn{S(t) = median({X_j | t - T/2 \le t_j \le t + T/2}).}
#'
#' For the initial time points where the time difference between the first data
#' point and the time point for which we are calculating the median is less than
#' half the \code{width}, we do not have enough data points to form a window
#' which has the same size to both sides of the time point. In this case fewer
#' data points are used to calculate the median and the time window is not
#' symmetric around the time point for which we are calculating the median.
#'
#' No median is calculated if the time difference between the last data point
#' and the current time point for which we are calculating the median is less
#' than half the \code{width}. We do not calculate the median using a smaller
#' time window so that the values do not change upon receiving new data.
#'
#' @inheritParams edecob
#' @param width The width of the window over which the moving median is taken in
#'   number of days.
#'
#' @return A data frame containing the values of the moving median, the study
#'   day to which it corresponds, the time window from which it was calculated,
#'   and the subject id corresponding to the data.
#' @export
#'


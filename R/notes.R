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

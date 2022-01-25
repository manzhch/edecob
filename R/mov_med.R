
#' Moving Median over a Time Window
#'
#' Calculates the moving median over a time window around a time point for the
#' all time points between the first and last time point provided.
#'
#' Consider a sample \eqn{X_1,\dots, X_n} of size \eqn{n} and the
#' reordering \eqn{X_{(1)},\dots, X_{(n)}} such
#' that \eqn{X_{(1)} \le X_{(2)} \le \dots \le X_{(n)}}, commonly
#' called the order statistic. Then for \eqn{n} even the median usually
#' defined as \deqn{median(X_1,\dots, X_n) = X_{(k)}, \mathrm{where } \;  k = n/2.} In the
#' case where \eqn{n} is odd the median is
#' defined as \deqn{median(X_1,\dots, X_n) = 1/2(X_{(k)} + X_{(k+1)}), \mathrm{where } \;  k = n/2.} Let the
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
#' time window so that the values do not change upon receiving new data with
#' time points newer than that of the old data.
#'
#' @inheritParams edecob
#' @param med_win A vector containing two numbers specifying the window over
#'   which the moving median is to be taken. More specifically, when given a
#'   certain time point, the numbers specify the maximum number of time units
#'   difference that a time point can have such that that data point will be
#'   considered for the moving median. Note that the first number must be smaller
#'   than the second number. The numbers may be negative.
#' @param min_pts_in_win The minimal number of measurements required to be in
#'   the time window in order for the median to be calculated.
#'
#' @return A data frame containing the values of the moving median, the study
#'   day to which it corresponds, the time window from which it was calculated,
#'   and the subject id corresponding to the data.
#' @export
#'
mov_med <- function(data,
                    med_win= c(-42, 42),
                    min_pts_in_win = 1) {

  source <- data$source[1]
  # the number of days for which the moving median will be calculated
  dur <- max(data$time_point) - min(data$time_point) + 1 - med_win[2]

  if (dur < 0) {
    dur <- 0
  }

  med_pts_med <- numeric(dur)
  med_pts_time_point <- numeric(dur)
  med_pts_win_beg <- numeric(dur)
  med_pts_win_end <- numeric(dur)
  med_pts_source <- character(dur)
  ll <- 1 # index for the vectors above



  # calculate moving medians of the first 6 weeks by taking smaller window
  # by taking whatever we have from baseline until current data$time_point + 6 weeks
  # (basically creating an unbalanced window)

  first_time_point <- min(data$time_point)
  last_time_point <- max(data$time_point)

  med_time_point <- first_time_point
  while (med_time_point <= min(last_time_point, last_time_point - med_win[2])) {
    win_ind <- as.logical((data$time_point >= max(first_time_point, med_time_point + med_win[1])) *
                          (data$time_point <= med_time_point + med_win[2]))

    win <- data$value[win_ind]

    if (sum(win_ind) >= min_pts_in_win) {

      # saving the values for the window
      med_pts_med[ll] <- stats::median(win)
      med_pts_time_point[ll] <- med_time_point
      med_pts_win_beg[ll] <- med_time_point + med_win[1]
      med_pts_win_end[ll] <- med_time_point + med_win[2]
      med_pts_source[ll] <- source

      ll <- ll + 1
    }
    med_time_point <- med_time_point + 1

  }


  # compile median point data into dataframe
  med_pts <- data.frame(
    "source" = med_pts_source,
    "time_point" = med_pts_time_point,
    "value" = med_pts_med,
    "win_beg" = med_pts_win_beg,
    "win_end" = med_pts_win_end,
    stringsAsFactors = FALSE
  )
  med_pts <- med_pts[med_pts$source != "", ]

  return(med_pts)
}





mov_mean <- function(data,
                     mean_win= c(-42, 42),
                     min_pts_in_win = 1) {

  source <- data$source[1]
  # the number of days for which the moving median will be calculated
  dur <- max(data$time_point) - min(data$time_point) + 1 - mean_win[2]

  if (dur < 0) {
    dur <- 0
  }

  mean_pts_med <- numeric(dur)
  mean_pts_time_point <- numeric(dur)
  mean_pts_win_beg <- numeric(dur)
  mean_pts_win_end <- numeric(dur)
  mean_pts_source <- character(dur)
  ll <- 1 # index for the vectors above



  # calculate moving medians of the first 6 weeks by taking smaller window
  # by taking whatever we have from baseline until current data$time_point + 6 weeks
  # (basically creating an unbalanced window)

  first_time_point <- min(data$time_point)
  last_time_point <- max(data$time_point)

  mean_time_point <- first_time_point
  while (mean_time_point <= min(last_time_point, last_time_point - mean_win[2])) {
    win_ind <- as.logical((data$time_point >= max(first_time_point, mean_time_point + mean_win[1])) *
                            (data$time_point <= mean_time_point + mean_win[2]))

    win <- data$value[win_ind]

    if (sum(win_ind) >= min_pts_in_win) {

      # saving the values for the window
      mean_pts_med[ll] <- mean(win)
      mean_pts_time_point[ll] <- mean_time_point
      mean_pts_win_beg[ll] <- mean_time_point + mean_win[1]
      mean_pts_win_end[ll] <- mean_time_point + mean_win[2]
      mean_pts_source[ll] <- source

      ll <- ll + 1
    }
    mean_time_point <- mean_time_point + 1

  }


  # compile median point data into dataframe
  mean_pts <- data.frame(
    "source" = mean_pts_source,
    "time_point" = mean_pts_time_point,
    "value" = mean_pts_med,
    "win_beg" = mean_pts_win_beg,
    "win_end" = mean_pts_win_end,
    stringsAsFactors = FALSE
  )
  mean_pts <- mean_pts[mean_pts$source != "", ]

  return(mean_pts)
}

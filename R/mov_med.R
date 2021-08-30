
#' Moving Median over a Time Window
#'
#' Calculates the moving median over a time window around a time point for the
#' all time points between the first and last study day provided.
#'
#' Consider a sample \eqn{X₁,\dots, Xₙ} of size \eqn{n} and the
#' reordering \eqn{X₍₁₎,\dots, X₍ₙ₎} such
#' that \eqn{X₍₁₎ \le X₍₂₎ \le \dots \le X₍ₙ₎}, commonly
#' called the order statistic. Then for \eqn{n} even the median usually
#' defined as \deqn{median(X₁,\dots, Xₙ) = X₍ₖ₎, where k = n/2.} In the
#' case where \eqn{n} is odd the median is
#' defined as \deqn{median(X₁,\dots, Xₙ) = 1/2(X₍ₖ₎ + X₍ₖ₊₁₎), where k = n/2.} Let the
#' study days at which the measurements \eqn{X₁, \dots, Xₙ} were taken
#' be \eqn{t₁, \dots, tₙ}.
#' Let \eqn{T} a fixed positive amount of time. Then the moving median at time
#' point \eqn{t} with window size \eqn{T} is defined as
#' \deqn{S(t) = median({Xⱼ | t - T/2 \le tⱼ \le t + T/2}).}
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
mov_med <- function(data,
                    width = 12*7) {


  subj_id <- data$subj_id[1]

  # the number of days for which the moving median will be calculated
  dur <- max(data$study_day) - min(data$study_day) + 1 - width/2

  med_pts_med <- numeric(dur)
  med_pts_study_day <- numeric(dur)
  med_pts_win_beg <- numeric(dur)
  med_pts_win_end <- numeric(dur)
  med_pts_subj_id <- character(dur)
  ll <- 1 # index for the vectors above



  # calculate moving medians of the first 6 weeks by taking smaller window
  # by taking whatever we have from baseline until current data$study_day + 6 weeks
  # (basically creating an unbalanced window)

  first_study_day <- min(data$study_day)
  win_beg_day <- first_study_day

  while (win_beg_day < first_study_day + width / 2) {

    # determine window
    win_ind <- as.logical((data$study_day >= first_study_day) *
      (data$study_day < win_beg_day + width / 2))

    win <- data$value[win_ind]

    if (sum(win_ind) > 0) {

      # saving the values for the window
      med_pts_med[ll] <- stats::median(win)
      med_pts_study_day[ll] <- win_beg_day
      med_pts_win_beg[ll] <- win_beg_day
      med_pts_win_end[ll] <- win_beg_day + width - 1
      med_pts_subj_id[ll] <- subj_id

      ll <- ll + 1
    }
    win_beg_day <- win_beg_day + 1
  }

  # calculate moving median for the rest of the data
  win_beg_day <- first_study_day
  last_study_day <- max(data$study_day)
  while (win_beg_day < last_study_day - width) {

    # determining winow
    win_ind <- as.logical(
      (data$study_day >= win_beg_day) * (data$study_day < win_beg_day + width))
    win <- data$value[win_ind]

    if (sum(win_ind) > 0 &&
      win_beg_day + width / 2 >= first_study_day) {

      # saving the values for the win
      med_pts_med[ll] <- stats::median(win)
      med_pts_study_day[ll] <- ceiling(win_beg_day + width / 2)
      med_pts_win_beg[ll] <- win_beg_day
      med_pts_win_end[ll] <- win_beg_day + width - 1
      med_pts_subj_id[ll] <- subj_id

      ll <- ll + 1
    }

    win_beg_day <- win_beg_day + 1
  }


  # compile median point data into dataframe
  med_pts <- data.frame(
    "subj_id" = med_pts_subj_id,
    "study_day" = med_pts_study_day,
    "value" = med_pts_med,
    "win_beg" = med_pts_win_beg,
    "win_end" = med_pts_win_end,
    stringsAsFactors = FALSE
  )
  med_pts <- med_pts[med_pts$subj_id != "", ]
}

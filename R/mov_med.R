
mov_med <- function(data,
                    date,
                    subj_id = "subj1",
                    width = as.difftime(12, units = "weeks")) {

  # check if data is sorted by time
  # check if data and date same length
  if (length(data) != length(date)) {
    warning("Data and date length mismatch")
  }

  dur <- as.numeric(
    max(date)
    - min(date)
    + as.difftime(1, units = "days")
    - as.difftime(6, units = "weeks"),
    units = "days")
  med_pts_med <- numeric(dur)
  med_pts_date <- as.Date(character(dur), "%Y-%m-%d")
  med_pts_win_beg <- as.Date(character(dur), "%Y-%m-%d")
  med_pts_win_end <- as.Date(character(dur), "%Y-%m-%d")
  med_pts_subj_id <- character(dur)
  ll <- 1 # index for the lists above



  # calculate moving medians of the first 6 weeks by taking smaller window
  # by taking whatever we have from baseline until current date + 6 weeks
  # (basically creating an unbalanced window)

  first_date <- min(date)
  win_beg_day <- first_date

  while (win_beg_day < first_date + width / 2) {

    # determine window
    win_ind <- as.logical((date >= first_date) *
      (date < win_beg_day + width / 2))
    win <- data[win_ind]

    if (sum(win_ind) > 0) {

      # saving the values for the window
      med_pts_med[ll] <- stats::median(win)
      med_pts_date[ll] <- win_beg_day
      med_pts_win_beg[ll] <- win_beg_day
      med_pts_win_end[ll] <- win_beg_day + 2 * width - 1
      med_pts_subj_id[ll] <- subj_id

      ll <- ll + 1
    }
    win_beg_day <- win_beg_day + 1
  }

  # calculate moving median for the rest of the data
  win_beg_day <- first_date
  last_date <- max(date)
  while (win_beg_day < last_date - width) {

    # determining winow
    win_ind <- as.logical(
      (date >= win_beg_day) * (date < win_beg_day + width))
    win <- data[win_ind]

    if (sum(win_ind) > 0 &&
      win_beg_day + width / 2 >= first_date) {

      # saving the values for the win
      med_pts_med[ll] <- stats::median(win)
      med_pts_date[ll] <- win_beg_day + width / 2
      med_pts_win_beg[ll] <- win_beg_day
      med_pts_win_end[ll] <- win_beg_day + 2 * width - 1
      med_pts_subj_id[ll] <- subj_id

      ll <- ll + 1
    }

    win_beg_day <- win_beg_day + 1
  }


  # compile median point data into dataframe
  med_pts <- data.frame(
    "pts" = med_pts_med,
    "win_beg" = med_pts_win_beg,
    "win_end" = med_pts_win_end,
    "date" = med_pts_date,
    "subj_id" = med_pts_subj_id,
    stringsAsFactors = FALSE
  )
  med_pts <- med_pts[!is.na(med_pts$date), ]
}

mov_med_resid <- function(data, date, smoother_pts) {

  resid <- numeric(length(data))

  if (length(data) > 0 && nrow(smoother_pts) > 0) {

    # calculate residuals for those dates where there is both a
    # smoother and a data point
    for (ii in 1:nrow(smoother_pts)) {
      dataset_ind <- as.logical(date == smoother_pts$date[ii])
      resid[dataset_ind] <- data[dataset_ind] - smoother_pts$pts[ii]
    }

    return(resid)

  } else {
    warning("Cannot calculate residuals (no data or smoother points)")
  }
}

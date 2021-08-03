
mov_med <- function(data,
                    study_day,
                    subj_id = "subj1",
                    width = 12*7) {

  # check if data is sorted by time
  # check if data and study_day same length
  if (length(data) != length(study_day)) {
    warning("Data and study_day length mismatch")
  }

  dur <- max(study_day) - min(study_day) + 1 - width/2

  med_pts_med <- numeric(dur)
  med_pts_study_day <- numeric(dur)
  med_pts_win_beg <- numeric(dur)
  med_pts_win_end <- numeric(dur)
  med_pts_subj_id <- character(dur)
  ll <- 1 # index for the lists above



  # calculate moving medians of the first 6 weeks by taking smaller window
  # by taking whatever we have from baseline until current study_day + 6 weeks
  # (basically creating an unbalanced window)

  first_study_day <- min(study_day)
  win_beg_day <- first_study_day

  while (win_beg_day < first_study_day + width / 2) {

    # determine window
    win_ind <- as.logical((study_day >= first_study_day) *
      (study_day < win_beg_day + width / 2))
    win <- data[win_ind]

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
  last_study_day <- max(study_day)
  while (win_beg_day < last_study_day - width) {

    # determining winow
    win_ind <- as.logical(
      (study_day >= win_beg_day) * (study_day < win_beg_day + width))
    win <- data[win_ind]

    if (sum(win_ind) > 0 &&
      win_beg_day + width / 2 >= first_study_day) {

      # saving the values for the win
      med_pts_med[ll] <- stats::median(win)
      med_pts_study_day[ll] <- win_beg_day + width / 2
      med_pts_win_beg[ll] <- win_beg_day
      med_pts_win_end[ll] <- win_beg_day + width - 1
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
    "study_day" = med_pts_study_day,
    "subj_id" = med_pts_subj_id,
    stringsAsFactors = FALSE
  )
  med_pts <- med_pts[med_pts$subj_id != "", ]
}

mov_med_resid <- function(data, study_day, smoother_pts) {

  resid <- numeric(length(data))

  if (length(data) > 0 && nrow(smoother_pts) > 0) {

    # calculate residuals for those study_days where there is both a
    # smoother and a data point
    for (ii in 1:nrow(smoother_pts)) {
      dataset_ind <- as.logical(study_day == smoother_pts$study_day[ii])
      resid[dataset_ind] <- data[dataset_ind] - smoother_pts$pts[ii]
    }

    return(resid)

  } else {
    warning("Cannot calculate residuals (no data or smoother points)")
  }
}

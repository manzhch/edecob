
mov_med <- function(dataa,
                    date,
                    pat_id,
                    width = as.difftime(12, units = "weeks")) {

  # check if data is sorted by time
  # check if data and date same length

  dur <- as.numeric(
    max(data)
    - min(data)
    + as.difftime(1, units = "days")
    - as.difftime(6, units = "weeks"),
    units = "days")
  med_pts_med <- numeric(dur)
  med_pts_date <- as.Date(character(dur), "%Y-%m-%d")
  med_pts_win_beg <- as.Date(character(dur), "%Y-%m-%d")
  med_pts_win_end <- as.Date(character(dur), "%Y-%m-%d")
  med_pts_USUBJID <- character(dur)
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
      med_pts_USUBJID[ll] <- pat_id

      ll <- ll + 1
    }
    win_beg_day <- win_beg_day + 1
  }

  # calculate moving median for the rest of the data
  win_beg_day <- first_date
  last_date <- max(date)
  while (win_beg_day < last_date - width) {

    # determining win
    win_ind <- as.logical((date >= win_beg_day) *
      (date < win_beg_day + width))
    win_one <- one_pat[win_one_ind, ]

    if (sum(win_one_ind) > 0 &&
      win_beg_day + width / 2 >= first_date) {

      # saving the values for the win
      med_pts_QRSRESN[l] <- med(win_one$QRSRESN)
      # empty_row$date <-  med(win_one$date)
      med_pts_date[l] <- win_beg_day + width / 2
      med_pts_win_beg[l] <- win_beg_day
      med_pts_win_end[l] <- win_beg_day + 2 * width - 1
      med_pts_USUBJID[l] <- my_pat
      # if(length(intersect(which(FUTT$USUBJID == my_pat), which(FUTT$date == win_beg_day + 2*width-1))) > 0){
      #   med_pts_residuals[l] <- FUTT$QRSRESN[
      #     intersect(which(FUTT$USUBJID == my_pat),
      #               which(FUTT$date == win_beg_day + 2*width-1))] -
      #     med_pts_QRSRESN[l]
      # }else{
      #   med_pts_residuals[l] <- NA
      # }


      l <- l + 1
    }

    win_beg_day <- win_beg_day + 1
  }


  # compile med point data into dataframe
  med_pts <- data.frame(
    "QRSRESN" = med_pts_QRSRESN,
    "win_beg" = med_pts_win_beg,
    "win_end" = med_pts_win_end,
    "date" = med_pts_date,
    "pat_id" = med_pts_pat_id,
    # "residuals" = med_pts_residuals,
    stringsAsFactors = FALSE
  )
  med_pts <- med_pts[!is.na(med_pts$date), ]
}

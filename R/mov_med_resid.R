
mov_med_resid <- function(data, date, med_pts) {

  resid <- numeric(length(data))

  if (nrow(data) > 0 && nrow(med_pts) > 0) {
    for (ii in 1:nrow(med_pts)) {
      dataset_ind <- as.logical(date == med_pts$date[ii])
      resid[dataset_ind] <- data[dataset_ind] - med_pts$med[ii]
    }
  }
}

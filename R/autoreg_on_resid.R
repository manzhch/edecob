# fits autoregression model on residuals
#
# resid: the residuals of the smoother
#
# output: the autoregression model
autoreg_on_resid <- function(resid){

  # fit AR model
  data_ind <- !is.na(resid)
  if (sum(data_ind) > 1 && stats::var(resid[data_ind]) != 0) {
    ar_resid <- stats::ar(resid[data_ind])
    return(ar_resid)
  }else{
    return(NA)
  }
}

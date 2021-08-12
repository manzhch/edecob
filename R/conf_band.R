
# takes the bootstrapped smoother and checks what the ratio is of
# bootstrapped smoothers in the pointwise alpha_p quantiles
ratio_in_ci <- function(alpha_p,
                        bt_smoother,
                        smoother_pts, # bt_smoother/bt_tot_rep and smoother_pts should have the same size
                        bt_tot_rep,
                        conf_band_lvl) {
  # initialize variables
  uniq_study_day <- unique(bt_smoother$study_day)
  uniq_study_day <- uniq_study_day[!is.na(uniq_study_day)]
  my_quantile_lower <- numeric(length(uniq_study_day))
  my_quantile_upper <- numeric(length(uniq_study_day))
  bt_smoother$init_est <- numeric(nrow(bt_smoother))

  # calculate difference of bootstrapped smoother vs non-bootstrapped smoother
  bt_smoother$init_est <- bt_smoother$value -
    rep(smoother_pts$pts[smoother_pts$study_day %in% uniq_study_day], bt_tot_rep)

  # calculate quantiles
  for (ii in 1:length(uniq_study_day)) {
    # print(bt_smoother$init_est[bt_smoother$study_day == uniq_study_day[[ii]]])
    my_quantile_lower[ii] <-
      stats::quantile(bt_smoother$init_est[bt_smoother$study_day == uniq_study_day[ii]], alpha_p)
    my_quantile_upper[ii] <-
      stats::quantile(bt_smoother$init_est[bt_smoother$study_day == uniq_study_day[ii]],
               1 - alpha_p)
  }
  names(my_quantile_lower) <- uniq_study_day
  names(my_quantile_upper) <- uniq_study_day

  # determine for each data point whether it is between the quantiles
  bt_smoother$in_intvl <- logical(nrow(bt_smoother))
  bt_smoother$in_intvl <-
    as.logical(pmin(bt_smoother$init_est <= rep(my_quantile_upper, bt_tot_rep),
                    bt_smoother$value >= rep(my_quantile_lower, bt_tot_rep)))

  # determine for each bootstrap whether all data points are between the quantiles
  bt_smoother$bt_in_intvl <- logical(nrow(bt_smoother))
  for (ii in 1:bt_tot_rep) {
    bt_smoother$bt_in_intvl[bt_smoother$bt_rep == ii] <-
      as.logical(min(bt_smoother$in_intvl[bt_smoother$bt_rep == ii]))
  }

  bt_in_intvl_ratio <-
    sum(bt_smoother$bt_in_intvl[bt_smoother$study_day == min(bt_smoother$study_day)]) /
    bt_tot_rep


  return(bt_in_intvl_ratio - conf_band_lvl)
}

# looking for the right alpha_p such that 95% of the bootstrap curves
# are within CI
find_ptw_conf_band_lvl <- function(bt_smoother,
                           smoother_pts,
                           bt_tot_rep,
                           conf_band_lvl){
  # print(ratio_in_ci(0, bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl))
  # print(ratio_in_ci(0.05, bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl))
  if (nrow(bt_smoother) > bt_tot_rep &&
      nrow(bt_smoother[!is.na(bt_smoother$value), ]) > 0 &&
      sign(ratio_in_ci(0, bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)) !=
      sign(ratio_in_ci(0.05, bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl))) {

    return(stats::uniroot(ratio_in_ci, c(0,0.05),#lower = 0, upper = 0.05,
                          bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl))

  } else if (nrow(bt_smoother) > bt_tot_rep &&
             nrow(bt_smoother[!is.na(bt_smoother$value), ]) > 0 &&
             sign(ratio_in_ci(0, bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)) ==
             sign(ratio_in_ci(0.05, bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl))) {
    # print(2)
    return(NA)
    # print("sign same")
  } else if (nrow(bt_smoother) > bt_tot_rep &&
             nrow(bt_smoother[!is.na(bt_smoother$value), ]) == 0) {
    return(NA)
    # print("many NA")
  } else {
    # print(3)
    return(NA)
  }
}



#' Confidence Bounds of the Smoother
#'
#' Calculate the confidence bounds of the smoother function using the bootstrap.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ qₓ(tᵢ), q₁₋ₓ(tᵢ) i = 1,\dots, N}
#'
#'     where \deqn{qₓ(tᵢ) = inf {u; P*[S(tᵢ)*ᵦ - S(tᵢ) \le u] \ge x} } is a
#'     pointwise bootstrap quantile, \eqn{S(tᵢ)*ᵦ the bootstrapped smoother},
#'     and \eqn{N} the number of measurements.
#'   \item We vary the pointwise error \eqn{2x} until \deqn{P*[qₓ(tᵢ) \le S(tᵢ)*ᵦ - S(tᵢ) \le q₁₋ₓ(tᵢ) ⩝ i = 1,\dots, N] ≈ 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[qₓ(tᵢ), q₁₋ₓ(tᵢ)]} is approximately \eqn{1-\alpha}.
#'   \item We let
#'   \deqn{ Iₙ(tᵢ) = [S(tᵢ) +  qₓ(tᵢ), S(tᵢ) + q₁₋ₓ(tᵢ)] ⩝ i = 1, \dots, N}
#'   confidence bounds. Then \eqn{{Iₙ(tᵢ); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
#'
#'}
#'
#'
#' @inheritParams edecob
#' @param bt_smoother The bootstrapped smoother. Use the output of \code{bt_smoother()}.
#'
#' @return A data frame containing the upper confidence bound, the lower confidence bound,
#'   and the study day corresponding to the bounds.
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#' @examples
conf_band <- function(bt_smoother,
                      smoother_pts,
                      bt_tot_rep,
                      conf_band_lvl){

  # calculate pointwise quantile
  study_days <- sort(unique(bt_smoother$study_day))
  my_quantile_lower <- numeric(length(study_days))
  my_quantile_upper <- numeric(length(study_days))
  width <- smoother_pts$win_end[1] - smoother_pts$win_beg[1] + 1

  ptw_conf_band_lvl <- find_ptw_conf_band_lvl(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)
  init_est <- bt_smoother$value -
    rep(smoother_pts$pts[smoother_pts$study_day %in% bt_smoother$study_day], bt_tot_rep)
  names(init_est) <- bt_smoother$study_day

  # calculate quantiles
  if (nrow(bt_smoother) > 0) {
    both_quantiles <-
      sapply(split(init_est, factor(names(init_est))),
             function(x){
               if (sum(!is.na(x)) > 0 &&
                   !is.na(ptw_conf_band_lvl) &&
                   !is.null(ptw_conf_band_lvl$root)) {

                 return(c(stats::quantile(x, ptw_conf_band_lvl$root),
                          stats::quantile(x, 1 - ptw_conf_band_lvl$root)))
               } else {
                  return(NA)
               }
             })
  } else {
    both_quantiles <- NA
  }
  suppressWarnings(
  if (!is.na(both_quantiles)) {
    my_quantile_lower <- both_quantiles[1, ]
    my_quantile_upper <- both_quantiles[2, ]
  })

  if (!is.null(study_days)) {
    names(my_quantile_lower) <- study_days
    names(my_quantile_upper) <- study_days

  }

  study_days_matching <- unique(sort(smoother_pts$study_day[
    smoother_pts$study_day %in% unique(bt_smoother$study_day)]))
  study_days_matching <- study_days_matching[study_days_matching <= max(smoother_pts$study_day) - width/2]

  for (ii in study_days_matching) {
    current_ind <- smoother_pts$study_day == ii
    smoother_pts$intvl_upper[current_ind] <-
      smoother_pts$pts[current_ind] + my_quantile_upper[as.character(ii)]
    smoother_pts$intvl_lower[current_ind] <-
      smoother_pts$pts[current_ind] + my_quantile_lower[as.character(ii)]
  }

  conf_band <- data.frame(
    upper = smoother_pts$intvl_upper[!is.na(smoother_pts$intvl_upper)],
    lower = smoother_pts$intvl_lower[!is.na(smoother_pts$intvl_lower)],
    study_day = study_days_matching)

  return(conf_band)
}



# takes the bootstrapped smoother and checks what the ratio is of
# bootstrapped smoothers in the pointwise alpha_p quantiles
ratio_in_ci <- function(alpha_p,
                        bt_smoother,
                        smoother_pts, # bt_smoother/bt_tot_rep and smoother_pts should have the same size
                        bt_tot_rep,
                        conf_band_lvl) {
  # initialize variables
  uniq_time_point <- unique(bt_smoother$time_point)
  uniq_time_point <- uniq_time_point[!is.na(uniq_time_point)]
  my_quantile_lower <- numeric(length(uniq_time_point))
  my_quantile_upper <- numeric(length(uniq_time_point))
  bt_smoother$init_est <- numeric(nrow(bt_smoother))

  # calculate difference of bootstrapped smoother vs non-bootstrapped smoother
  bt_smoother$init_est <- bt_smoother$value -
    rep(smoother_pts$value[smoother_pts$time_point %in% uniq_time_point], bt_tot_rep)

  # calculate quantiles
  for (ii in 1:length(uniq_time_point)) {
    my_quantile_lower[ii] <-
      stats::quantile(bt_smoother$init_est[bt_smoother$time_point == uniq_time_point[ii]], alpha_p)
    my_quantile_upper[ii] <-
      stats::quantile(bt_smoother$init_est[bt_smoother$time_point == uniq_time_point[ii]],
               1 - alpha_p)
  }
  names(my_quantile_lower) <- uniq_time_point
  names(my_quantile_upper) <- uniq_time_point

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
    sum(bt_smoother$bt_in_intvl[bt_smoother$time_point == min(bt_smoother$time_point)]) /
    bt_tot_rep


  return(bt_in_intvl_ratio - conf_band_lvl)
}

# looking for the right alpha_p such that 95% of the bootstrap curves
# are within CI
find_ptw_conf_band_lvl <- function(bt_smoother,
                           smoother_pts,
                           bt_tot_rep,
                           conf_band_lvl){

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
    return(NA)
  } else if (nrow(bt_smoother) > bt_tot_rep &&
             nrow(bt_smoother[!is.na(bt_smoother$value), ]) == 0) {
    return(NA)
  } else {
    return(NA)
  }
}



#' Confidence Bounds of the Smoother
#'
#' Calculate the confidence bounds of the smoother function using the bootstrap.
#'
#' The procedure is as follows:
#' \enumerate{
#'   \item We compute the quantiles \deqn{ q_x(t_i), q_{1-x}(t_i) i = 1,\dots, N}
#'     where \deqn{q_x(t_i) = inf\left\{u; P^*[S(t_i)^*_b - S(t_i) \le u] \ge x\right\} } is a
#'     pointwise bootstrap quantile, \eqn{S(t_i)^*_b} the bootstrapped smoother,
#'     and \eqn{N} the number of measurements or rows in \code{data}, in our case the number of rows.
#'   \item We vary the pointwise error \eqn{x} until \deqn{P^*[q_x(t_i) \le S(t_i)^*_b - S(t_i) \le q_{1-x}(t_i) \forall i = 1,\dots, N] \approx 1-\alpha.}
#'     In other words, until the ratio of bootstrap curves that have all their points within
#'     \eqn{[q_x(t_i), q_{1-x}(t_i)]} is approximately \eqn{1-\alpha}.
#'   \item We define
#'   \deqn{ I_n(t_i) = [S(t_i) +  q_x(t_i), S(t_i) + q_{1-x}(t_i)] \forall i = 1,\dots, N}
#'   the confidence bounds. Then \eqn{{I_n(t_i); i = 1,\dots, N}} is a consistent simultaneous confidence band of level \eqn{1-\alpha}.
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
#'   and the time point corresponding to the bounds.
#' @export
#'
#' @references Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in
#'   Nonstationary Time Series. \emph{The Annals of Statistics}, 26(1), 48-83.
#'
#'
conf_band <- function(bt_smoother,
                      smoother_pts,
                      bt_tot_rep,
                      conf_band_lvl){

  # calculate pointwise quantile
  time_points <- sort(unique(bt_smoother$time_point))
  my_quantile_lower <- numeric(length(time_points))
  my_quantile_upper <- numeric(length(time_points))
  width <- smoother_pts$win_end[1] - smoother_pts$win_beg[1] + 1

  ptw_conf_band_lvl <- find_ptw_conf_band_lvl(bt_smoother, smoother_pts, bt_tot_rep, conf_band_lvl)
  init_est <- bt_smoother$value -
    rep(smoother_pts$value[smoother_pts$time_point %in% bt_smoother$time_point], bt_tot_rep)
  names(init_est) <- bt_smoother$time_point

  # calculate quantiles
  if (nrow(bt_smoother) > 0) {
    both_quantiles <-
      sapply(split(init_est, factor(unique(names(init_est)), unique(names(init_est)))),
             function(x){
               if (sum(!is.na(x)) > 0 &&
                   !is.na(ptw_conf_band_lvl)[1] &&
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
  if (!is.na(both_quantiles)[1]) {
    my_quantile_lower <- both_quantiles[1, ]
    my_quantile_upper <- both_quantiles[2, ]
  })

  if (!is.null(time_points)) {
    names(my_quantile_lower) <- time_points
    names(my_quantile_upper) <- time_points

  }

  time_points_matching <- unique(sort(smoother_pts$time_point[
    smoother_pts$time_point %in% unique(bt_smoother$time_point)]))
  time_points_matching <- time_points_matching[time_points_matching <= max(smoother_pts$time_point) - width/2]

  for (ii in time_points_matching) {
    current_ind <- smoother_pts$time_point == ii
    smoother_pts$intvl_upper[current_ind] <-
      smoother_pts$value[current_ind] + my_quantile_upper[as.character(ii)]
    smoother_pts$intvl_lower[current_ind] <-
      smoother_pts$value[current_ind] + my_quantile_lower[as.character(ii)]
  }

  conf_band <- data.frame(
    source = rep(smoother_pts$source[1], length(time_points_matching)),
    time_point = time_points_matching,
    lower = smoother_pts$intvl_lower[!is.na(smoother_pts$intvl_lower)],
    upper = smoother_pts$intvl_upper[!is.na(smoother_pts$intvl_upper)]
    )

  return(conf_band)
}


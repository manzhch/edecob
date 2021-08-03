#' Title
#'
#' @param data
#' @param date
#' @param smoother_pts
#' @param conf_band
#' @param event
#' @param subj_id
#' @param learn_dur
#' @param basel_dur
#' @param thresh_diff
#'
#' @return A `ggplot2` that visualizes the data.
#' @export
#'
#' @examples
#'
#' @importFrom rlang .data
edecob_plot <- function(data,
                        date,
                        smoother_pts,
                        conf_band,
                        event,
                        subj_id = "subj1",
                        learn_dur,
                        basel_dur,
                        thresh_diff = 0.1,
                        width) {


  subj_data <- data.frame(data = data, date = date)

  # calculate baseline and threshold
  basel <- stats::median(data[as.logical(
    (date >= min(date) + learn_dur) *
      (date < min(date) + learn_dur + basel_dur))])
  thresh <- basel * (1 - thresh_diff)


  # plotting data points, baseline, and threshold
  plot_colors <- character()
  plot_shape <- numeric()
  plot_fill <- character()
  plot_size <- numeric()
  patient_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = subj_data,
      ggplot2::aes(x = .data$date, y = .data$data),
      fill = "black",
      show.legend = FALSE,
      shape = 21,
      color = "grey",
      size = 4
    ) +
    ggplot2::geom_point(
      data = subj_data[which(subj_data$date < min(subj_data$date) + learn_dur), ],
      ggplot2::aes(x = .data$date, y = .data$data),
      color = "grey70",
      size = 3.8
    ) +
    ggplot2::labs(x = "date", y = "5UTT average turn speed (rad/s)") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = basel, color = "basel")) +
    # geom_hline(aes(yintercept = threshold, color = "threshold")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = basel*(1 - thresh_diff), color = "threshold")) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = min(subj_data$date) + learn_dur,
                   linetype = "basel period"),
               color = "blue") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = min(subj_data$date) + learn_dur + width - 1, linetype = "basel period"),
               color = "blue") +
    # ggplot2::geom_vline(ggplot2::aes(xintercept = my_patient_devices$RETDT),
    #            color = "grey75") +
    ggplot2::scale_color_manual(
      "",
      aesthetics  = "color",
      values = c("basel" = "black", "threshold" = "red")
    ) +
    ggplot2::scale_color_manual("",
                       aesthetics = "linetype",
                       values = c("basel period" = "dashed")) +
    ggplot2::scale_color_manual(
      "",
      aesthetics = "fill",
      breaks = c("TRUE", "FALSE"),
      values = c("red", "black")
    ) +
    ggplot2::ggtitle(subj_id) +
    ggplot2::theme(legend.title = ggplot2::element_blank())


  # plotting median points for the moving median
  if (!is.null(nrow(smoother_pts)) &&
      nrow(smoother_pts) > 0) {
    patient_plot <- patient_plot +
      ggplot2::geom_point(
        data = smoother_pts,
        ggplot2::aes(x = .data$date, y = .data$pts, pch = "moving median"),
        color = "orange",
        size = 3
      )

    plot_colors <- append(plot_colors, "orange")
    plot_shape <- append(plot_shape, 16)
    plot_fill <- append(plot_fill, "orange")
    plot_size <- append(plot_size, 5)

  }

  # draw event point
  if (event[[1]]) {
    patient_plot <- patient_plot +
      ggplot2::geom_point(
        data = data.frame(date = event[[2]], data1 = smoother_pts$pts[which(smoother_pts$date == event[[2]])]),
        ggplot2::aes(x = .data$date, y = .data$data1, shape = "moving median event onset"),
        color = "red",
        size = 9
      )
    plot_colors <- append(plot_colors, "red")
    plot_shape <- append(plot_shape, 17)
    plot_fill <- append(plot_fill, "red")
    plot_size <- append(plot_size, 8)

  }


  # redraw basel and threshold so that they are on top
  patient_plot <- patient_plot +
    ggplot2::scale_shape_manual("", values = plot_shape,
                       guide = ggplot2::guide_legend(
                         override.aes = list(color = plot_colors, size = plot_size, fill = plot_fill)
                       )) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = basel*(1 - thresh_diff)), color = "red") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = basel), color = "black") +
    ggplot2::theme(text = ggplot2::element_text(size = 25),
          axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 20)),
          axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))


  # CI
  patient_plot <- patient_plot +
    ggplot2::geom_ribbon(data = conf_band, ggplot2::aes(x = .data$date, ymin = .data$lower, ymax = .data$upper),
                alpha = 0.4, color = "transparent", fill = "blue")



  return(patient_plot)
}


# plot <- function(object) {
#   UseMethod("plot", object)
# }

#' Plot Event Data
#'
#' Generates a ggplot2 visualizing the data and event if any detected.
#'
#' The data points are plotted on the x-axis with the time point on the y-axis.
#' The data from the learning period are gray. The baseline period is outlined
#' by two vertical blue lines. The smoother is plotted in orange. The confidence
#' bound  is the blue area. If an event is detected, a red triangle will mark
#' the detection day.
#'
#' Since the output is a ggplot, it can be manipulated using the usual functions
#' from the ggplot2 package like \code{xlim}, \code{ylim}, \code{scale_x_continuous},
#' and many more.
#'
#' @param x The output of the \code{edecob} function for one subject. It is an object
#'   of class \code{edecob} containing the data and the event information.
#' @param ... Other arguments like \code{title}, \code{xlab}, or \code{ylab}.
#'
#' @return A `ggplot2` object that visualizes the data.
#' @export
#'
#'
#' @importFrom rlang .data
#' @importFrom utils capture.output
#' @importFrom graphics plot.new
plot.edecob <- function(x, ...) {

  stopifnot("Cannot plot data for multiple sources. Please choose one source to plot using \"$\"." = !("event_info" %in% names(x)))

  event_data <- x

  # if ggplot2 was not imported
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"gglot2\" needed for plots.", call. = FALSE)
  } else {

    if ("title" %in% names(list(...))) {
      title <- list(...)$title
    } else {
      title <- event_data$data$source[1]
    }

    if ("xlab" %in% names(list(...))) {
      xlab <- list(...)$xlab
    } else {
      xlab <- event_data$col_names[2]
    }

    if (!("ylab" %in% names(list(...)))) {
      ylab <- list(...)$ylab
    } else {
      ylab <- event_data$col_names[3]
    }

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mar = c(1,1,1,1))
    graphics::plot.new()

    # initialize variables
    dot_size <- 1.8
    smo_size <- 1
    eve_size <- 3
    txt_size <- 12

    data <- event_data$data[!is.na(event_data$data$value), ]
    detec_upper <- event_data$detec_upper
    detec_lower <- event_data$detec_lower
    smoother_pts <- event_data$smoother_pts
    conf_band <- event_data$conf_band
    event <- event_data$event
    source <- event_data$data$source[1]

    subj_data <- data.frame(
      time_point = data[data$source == source, 2],
      data = data[data$source == source, 3]
    )

    # plotting data points, baseline, and threshold
    plot_colors <- character()
    plot_shape <- numeric()
    plot_fill <- character()
    plot_size <- numeric()
    patient_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = subj_data,
        ggplot2::aes(x = .data$time_point, y = .data$data),
        fill = "black",
        show.legend = FALSE,
        shape = 21,
        color = "transparent",
        size = dot_size,
        alpha = 0.8
      ) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = detec_lower, color = "Detection Bounds")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = detec_upper, color = "Detection Bounds")) +
      ggplot2::scale_color_manual(
        "",
        aesthetics  = "color",
        values = c(#"Lower Detection Bound" = "red",
                   "Detection Bounds" = "red",
                   "Smoother" = "orange")
      ) +
      ggplot2::scale_color_manual("",
        aesthetics = "fill",
        values = c("Confidence Band" = "blue")
      ) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.title = ggplot2::element_blank())


    # plotting median points for the moving median
    if (!is.null(nrow(smoother_pts)) && nrow(smoother_pts) > 0) {

      # create new data frame but containing NA values when no smoother point
      smoother_na <- data.frame(
        pts = rep(NA, max(smoother_pts$time_point) - min(smoother_pts$time_point) + 1),
        time_point = min(smoother_pts$time_point):max(smoother_pts$time_point))

      smoother_na$pts[smoother_na$time_point %in% smoother_pts$time_point] <-
        smoother_pts$value

      patient_plot <- patient_plot +
        ggplot2::geom_line(
          data = smoother_na,
          ggplot2::aes(x = .data$time_point,
                       y = .data$pts,
                       color = "Smoother"),
          size = smo_size, linetype = "solid"
        )
    }

    # draw event point
    if (event[[1]]) {
      patient_plot <- patient_plot +
        ggplot2::geom_point(
          data = data.frame(time_point = event[[2]], data1 = smoother_pts$value[which(smoother_pts$time_point == event[[2]])]),
          ggplot2::aes(x = .data$time_point, y = .data$data1, shape = "Event Onset"),
          color = "red",
          size = eve_size
        )
      plot_colors <- append(plot_colors, "red")
      plot_shape <- append(plot_shape, 17)
      plot_fill <- append(plot_fill, "red")
      plot_size <- append(plot_size, eve_size)

    }


    # redraw basel and threshold so that they are on top
    patient_plot <- patient_plot +
      ggplot2::scale_shape_manual("", values = plot_shape,
                                  guide = ggplot2::guide_legend(
                                    override.aes = list(color = plot_colors, size = plot_size, fill = plot_fill)
                                  )) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = detec_upper), color = "red") #+


    # CI

    if (nrow(conf_band) > 0) {
      # create new data frame but containing NA values when no CI points
      conf_band_na <- data.frame(
        upper = rep(NA, max(conf_band$time_point) - min(conf_band$time_point) + 1),
        lower = rep(NA, max(conf_band$time_point) - min(conf_band$time_point) + 1),
        time_point = min(conf_band$time_point):max(conf_band$time_point))

      conf_band_na$upper[conf_band_na$time_point %in% conf_band$time_point] <-
        conf_band$upper
      conf_band_na$lower[conf_band_na$time_point %in% conf_band$time_point] <-
        conf_band$lower

      patient_plot <- patient_plot +
        ggplot2::geom_ribbon(data = conf_band_na, ggplot2::aes(x = .data$time_point, ymin = .data$lower, ymax = .data$upper, fill = "Confidence Band"),
                             alpha = 0.45, color = "transparent")
    }


    if (event[[1]]) {
      patient_plot <- patient_plot +
        ggplot2::geom_point(
          data = data.frame(time_point = event[[2]], data1 = smoother_pts$value[which(smoother_pts$time_point == event[[2]])]),
          ggplot2::aes(x = .data$time_point, y = .data$data1, shape = "Event Onset"),
          color = "red",
          size = eve_size
        )
    }


    # annotation below plot
    detec_finite <- min(data$value)
    if (!is.infinite(event_data$detec_lower) && !is.infinite(event_data$detec_upper)) {
      detec_finite <- c(event_data$detec_lower, event_data$detec_upper)
    } else if (!is.infinite(event_data$detec_lower) && is.infinite(event_data$detec_upper)) {
      detec_finite <- event_data$detec_lower
    } else if (is.infinite(event_data$detec_lower) && !is.infinite(event_data$detec_upper)) {
      detec_finite <- event_data$detec_upper
    }


    if (event$event_detected) {
      if (event$event_duration == 1) {
        text_below_plot_event <-
          paste("Event detected at ", event_data$time_unit, " ",
                event$event_onset, ", sustained for ", event$event_duration,
                " ", event_data$time_unit, ".", sep = "")
      } else {
        text_below_plot_event <-
          paste("Event detected at ", event_data$time_unit, " ",
                event$event_onset, ", sustained for ", event$event_duration,
                " ", event_data$time_unit, "s.", sep = "")
      }
    } else {
        if (event$event_duration == 1) {
          text_below_plot_event <-
            paste("Event detection censored at ", event_data$time_unit, " ", event$event_onset ,". Longest sustained change is ",
                  event$event_duration, " ", event_data$time_unit, ".", sep = "")

        } else {
          text_below_plot_event <-
            paste("Event detection censored at ", event_data$time_unit, " ", event$event_onset ,". Longest sustained change is ",
                  event$event_duration, " ", event_data$time_unit, "s.", sep = "")

        }
    }

    if (event_data$min_change_dur == 1) {
      text_below_plot_min_dur <-
        paste("Minimal duration of change for event detection: ",
              event_data$min_change_dur, " ", event_data$time_unit, sep = "")
    } else {
      text_below_plot_min_dur <-
        paste("Minimal duration of change for event detection: ",
              event_data$min_change_dur, " ", event_data$time_unit, "s", sep = "")
    }

    text_below_plot_CI <- paste("Level of confidence band: ", 100*event_data$conf_band_lvl, "%", sep = "")

    text_below_plot <- grid::textGrob(
      paste(text_below_plot_event, text_below_plot_min_dur, text_below_plot_CI, sep = "\n"),
      gp = grid::gpar(fontsize = 10),
      just = c("left", "top"))

    patient_plot <-  patient_plot +
      ggplot2::theme(plot.margin = ggplot2::unit(c(0.03,0.03,0.18,0.03), "npc")) +
      ggplot2::annotation_custom(grob = text_below_plot,
                                 xmin = min(data$time_point) - 0.15*(max(data$time_point) - min(data$time_point)),
                                 xmax = min(data$time_point) - 0.15*(max(data$time_point) - min(data$time_point)),
                                 ymin = min(c(data$value, detec_finite, conf_band_na$lower), na.rm = TRUE) - 0.29*(max(c(data$value, detec_finite)) - min(c(data$value, detec_finite))),
                                 ymax = min(c(data$value, detec_finite, conf_band_na$lower), na.rm = TRUE) - 0.29*(max(c(data$value, detec_finite)) - min(c(data$value, detec_finite))))


    # text if lower detection bound at -Inf
    if (is.infinite(detec_lower)) {

      inf_text <- grid::textGrob(
        paste("-", rawToChar(as.raw(c(226, 136, 158))), sep = ""),
        gp = grid::gpar(fontsize = 12, col = "red"),
        just = c("right", "top"))

      patient_plot <- patient_plot +
        ggplot2::annotation_custom(grob = inf_text,
                                   xmin = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   xmax = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   ymin = min(c(data$value, detec_finite, conf_band_na$lower), na.rm = TRUE) - 0.035*(max(data$value) - min(data$value)),
                                   ymax = min(c(data$value, detec_finite, conf_band_na$lower), na.rm = TRUE) - 0.035*(max(data$value) - min(data$value)))
    }

    # text if lower detection bound at Inf
    if (is.infinite(detec_upper)) {

      inf_text <- grid::textGrob(
        rawToChar(as.raw(c(226, 136, 158))),
        gp = grid::gpar(fontsize = 12, col = "red"),
        just = c("right", "bottom"))

      patient_plot <- patient_plot +
        ggplot2::annotation_custom(grob = inf_text,
                                   xmin = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   xmax = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   ymin = max(c(data$value, detec_finite, conf_band_na$upper), na.rm = TRUE) + 0.045*(max(data$value) - min(data$value)),
                                   ymax = max(c(data$value, detec_finite, conf_band_na$upper), na.rm = TRUE) + 0.045*(max(data$value) - min(data$value)))
    }

    # draw moving mean
    # patient_plot <- patient_plot +
    #   ggplot2::geom_line(
    #     data = mov_mean(event_data$data),
    #     ggplot2::aes(x = .data$time_point,
    #                  y = .data$value),
    #                 color = "green",
    #     size = 1, linetype = "solid"
    # )


    # Code to override clipping
    patient_plot <- ggplot2::ggplotGrob(patient_plot)
    patient_plot$layout$clip[patient_plot$layout$name == "panel"] <- "off"

    grid::grid.draw(patient_plot)
    graphics::par(mar = c(5, 4, 4, 2) + 0.1)

    invisible(capture.output(return(patient_plot)))
  }
}

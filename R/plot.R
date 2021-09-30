
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
plot.edecob <- function(x,
                        ...) {

  event_data <- x

  # if ggplot2 was not imported
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"gglot2\" needed for plots.", call. = FALSE)
  } else {

    if (!("title" %in% names(list(...)))) {
      title <- event_data$data$source[1]
    }

    if (!("xlab" %in% names(list(...)))) {
      xlab <- event_data$col_names[2]
    }

    if (!("ylab" %in% names(list(...)))) {
      ylab <- event_data$col_names[3]
    }

    par(mar = c(1,1,1,1))
    plot.new()

    # initialize variables
    dot_size <- 1.8
    smo_size <- 1
    eve_size <- 3
    txt_size <- 12

    data <- event_data$data
    detec_upper <- event_data$detec_upper
    detec_lower <- event_data$detec_lower
    smoother_pts <- event_data$smoother_pts
    conf_band <- event_data$conf_band
    event <- event_data$event
    source <- event_data$data$source[1]
    # basel_start <- event_data$basel_start
    # basel_end <- event_data$basel_end
    # width <- event_data$width

    subj_data <- data.frame(
      time_point = data$time_point[data$source == source],
      data = data$value[data$source == source]
    )

    # conf_band_text <- paste("Confidence Band (", event_data$conf_band_lvl*100, "%)", sep = "")
    # conf_band_text <- "Confidence Band"

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
      # ggplot2::geom_point(
      #   data = subj_data[which(subj_data$time_point < basel_start), ],
      #   ggplot2::aes(x = .data$time_point, y = .data$data),
      #   color = "grey70",
      #   size = dot_size*0.9
      # ) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = detec_lower, color = "Detection Bounds")) +
      # geom_hline(aes(yintercept = threshold, color = "threshold")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = detec_upper, color = "Detection Bounds")) +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = min(subj_data$time_point) + learn_dur,
      #                                  color = "Baseline Period"),
      #                     linetype = "dashed", key_glyph = "path") +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = basel_start,
      #                                  linetype = "Baseline Period"),
      #                     color = "blue") +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = basel_end
      #                                  ),
      #                     color = "blue", show.legend = F, linetype = "dashed") +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = my_patient_devices$RETDT),
      #            color = "grey75") +
      ggplot2::scale_color_manual(
        "",
        aesthetics  = "color",
        values = c(#"Lower Detection Bound" = "red",
                   "Detection Bounds" = "red",
                   "Smoother" = "orange")
      ) +
      # ggplot2::scale_color_manual("",
      #                             aesthetics = "linetype",
      #                             values = c("Baseline Period" = "dashed")) +
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

      # plot_colors <- append(plot_colors, "orange")
      # plot_shape <- append(plot_shape, 16)
      # plot_fill <- append(plot_fill, "orange")
      # plot_size <- append(plot_size, smo_size + 1)

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
      # ggplot2::geom_hline(ggplot2::aes(yintercept = basel), color = "black") #+
      # ggplot2::theme(text = ggplot2::element_text(size = txt_size),
      #                axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 20)),
      #                axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))


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
                                 ymin = min(c(data$value, detec_finite)) - 0.29*(max(c(data$value, detec_finite)) - min(c(data$value, detec_finite))),
                                 ymax = min(c(data$value, detec_finite)) - 0.29*(max(c(data$value, detec_finite)) - min(c(data$value, detec_finite))))

    # text if lower detection bound at -Inf
    if (is.infinite(detec_lower)) {

      inf_text <- grid::textGrob(
        "-∞",
        gp = grid::gpar(fontsize = 12, col = "red"),
        just = c("right", "top"))

      patient_plot <- patient_plot +
        ggplot2::annotation_custom(grob = inf_text,
                                   xmin = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   xmax = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   ymin = min(c(data$value, detec_finite)) - 0.035*(max(data$value) - min(data$value)),
                                   ymax = min(c(data$value, detec_finite)) - 0.035*(max(data$value) - min(data$value)))
    }

    # text if lower detection bound at Inf
    if (is.infinite(detec_upper)) {

      inf_text <- grid::textGrob(
        "∞",
        gp = grid::gpar(fontsize = 12, col = "red"),
        just = c("right", "bottom"))

      patient_plot <- patient_plot +
        ggplot2::annotation_custom(grob = inf_text,
                                   xmin = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   xmax = min(data$time_point) - 0.06*(max(data$time_point) - min(data$time_point)),
                                   ymin = max(c(data$value, detec_finite)) + 0.045*(max(data$value) - min(data$value)),
                                   ymax = max(c(data$value, detec_finite)) + 0.045*(max(data$value) - min(data$value)))
    }


    # print(patient_plot)
    #   ggplot2::labs(tag = paste("Minimal duration of change for event detection:", event_data$min_change_dur)) +
    #   ggplot2::theme(plot.tag.position = c(min(data$time_point), min(data$value) + 0.05*(max(data$value) - min(data$value))))
    # grid::grid.text((paste("Minimal duration of change for event detection:", event_data$min_change_dur)),
    #           x = grid::unit(0, "npc"), y = grid::unit(0, "npc"), just = c("left", "bottom"),
    #           gp = grid::gpar(fontsize = 12, col = "black"))
    # return(patient_plot)

    # Code to override clipping
    patient_plot <- ggplot2::ggplotGrob(patient_plot)
    patient_plot$layout$clip[patient_plot$layout$name == "panel"] <- "off"

    grid::grid.draw(patient_plot)

    invisible(capture.output(return(patient_plot)))
  }
}

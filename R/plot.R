
# plot <- function(object) {
#   UseMethod("plot", object)
# }

#' Plot Event Data
#'
#' Generates a ggplot2 visualizing the data and event if any detected.
#'
#' The data points are plotted on the x-axis with the study day on the y-axis.
#' The data from the learning period are gray. The baseline period is outlined
#' by two vertical blue lines. The smoother is plotted in orange. The confidence
#' bound  is the blue area. If an event is detected, a red triangle will mark
#' the detection day.
#'
#' Since the output is a ggplot, it can be manipulated using the usual functions
#' from the ggplot2 package like \code{xlim}, \code{ylim}, \code{scale_x_continuous},
#' and many more.
#'
#' @param event_data The output of the \code{edecob} function for one subject. It is an object
#'   of class \code{edecob} containing the data and the event information.
#' @param title The title of the plot. Defaults to the subject identifier for the
#'   subject to which \code{event_data} corresponds to.
#' @param xlab The label for the x-axis. Defaults to the name of the second
#'   column of the data when it was first entered into the \code{edecob} function.
#' @param ylab The label for the y-axis. Defaults to the name of the third
#'   column of the data when it was first entered into the \code{edecob} function.
#' @param ... Other arguments.
#'
#' @return A `ggplot2` object that visualizes the data.
#' @export
#'
#'
#' @importFrom rlang .data
#' @importFrom utils capture.output
#' @importFrom graphics plot.new
plot.edecob <- function(event_data,
                        ...,
                        title = event_data$data$subj_id[1],
                        xlab = event_data$col_names[2],
                        ylab = event_data$col_names[3]) {

  # if ggplot2 was not imported
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"gglot2\" needed for plots.", call. = FALSE)
  } else {


    plot.new()

    # initialize variables
    dot_size <- 1.8
    smo_size <- 1
    eve_size <- 3
    txt_size <- 12

    data <- event_data$data
    basel <- event_data$baseline
    thresh <- event_data$threshold
    smoother_pts <- event_data$smoother_pts
    conf_band <- event_data$conf_band
    event <- event_data$event
    subj_id <- event_data$data$subj_id[1]
    # basel_start <- event_data$basel_start
    # basel_end <- event_data$basel_end
    # width <- event_data$width

    subj_data <- data.frame(
      data = data$value[data$subj_id == subj_id],
      study_day = data$study_day[data$subj_id == subj_id])

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
        ggplot2::aes(x = .data$study_day, y = .data$data),
        fill = "black",
        show.legend = FALSE,
        shape = 21,
        color = "transparent",
        size = dot_size,
        alpha = 0.8
      ) +
      # ggplot2::geom_point(
      #   data = subj_data[which(subj_data$study_day < basel_start), ],
      #   ggplot2::aes(x = .data$study_day, y = .data$data),
      #   color = "grey70",
      #   size = dot_size*0.9
      # ) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = basel, color = "Baseline")) +
      # geom_hline(aes(yintercept = threshold, color = "threshold")) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = thresh, color = "Threshold")) +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = min(subj_data$study_day) + learn_dur,
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
        values = c("Baseline" = "black", "Threshold" = "red",
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
        pts = rep(NA, max(smoother_pts$study_day) - min(smoother_pts$study_day) + 1),
        study_day = min(smoother_pts$study_day):max(smoother_pts$study_day))

      smoother_na$pts[smoother_na$study_day %in% smoother_pts$study_day] <-
        smoother_pts$value

      patient_plot <- patient_plot +
        ggplot2::geom_line(
          data = smoother_na,
          ggplot2::aes(x = .data$study_day,
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
          data = data.frame(study_day = event[[2]], data1 = smoother_pts$value[which(smoother_pts$study_day == event[[2]])]),
          ggplot2::aes(x = .data$study_day, y = .data$data1, shape = "Event Onset"),
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
      ggplot2::geom_hline(ggplot2::aes(yintercept = thresh), color = "red") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = basel), color = "black") +
      ggplot2::theme(text = ggplot2::element_text(size = txt_size),
                     axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 20)),
                     axis.title.y.left = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))


    # CI

    # create new data frame but containing NA values when no CI points
    conf_band_na <- data.frame(
      upper = rep(NA, max(conf_band$study_day) - min(conf_band$study_day) + 1),
      lower = rep(NA, max(conf_band$study_day) - min(conf_band$study_day) + 1),
      study_day = min(conf_band$study_day):max(conf_band$study_day))

    conf_band_na$upper[conf_band_na$study_day %in% conf_band$study_day] <-
      conf_band$upper
    conf_band_na$lower[conf_band_na$study_day %in% conf_band$study_day] <-
      conf_band$lower

    patient_plot <- patient_plot +
      ggplot2::geom_ribbon(data = conf_band_na, ggplot2::aes(x = .data$study_day, ymin = .data$lower, ymax = .data$upper, fill = "Confidence Band"),
                           alpha = 0.45, color = "transparent")


    if (event[[1]]) {
      patient_plot <- patient_plot +
        ggplot2::geom_point(
          data = data.frame(study_day = event[[2]], data1 = smoother_pts$value[which(smoother_pts$study_day == event[[2]])]),
          ggplot2::aes(x = .data$study_day, y = .data$data1, shape = "Event Onset"),
          color = "red",
          size = eve_size
        )
    }


    # annotation below plot

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
                                 xmin = min(data$study_day) - 0.15*(max(data$study_day) - min(data$study_day)),
                                 xmax = min(data$study_day) - 0.15*(max(data$study_day) - min(data$study_day)),
                                 ymin = min(data$value) - 0.28*(max(data$value) - min(data$value)),
                                 ymax = min(data$value) - 0.28*(max(data$value) - min(data$value)))
    # print(patient_plot)
    #   ggplot2::labs(tag = paste("Minimal duration of change for event detection:", event_data$min_change_dur)) +
    #   ggplot2::theme(plot.tag.position = c(min(data$study_day), min(data$value) + 0.05*(max(data$value) - min(data$value))))
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

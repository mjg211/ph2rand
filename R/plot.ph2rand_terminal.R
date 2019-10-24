#' @export
plot.ph2rand_terminal <- function(x, k, output = F, summary = F, ...) {

  ##### Check inputs ###########################################################

  #check_terminal(x, "any")
  k  <- check_k(k, x$des)
  check_logical(output, "output")
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    #summary_plot_terminal(x)
  }

  ##### Perform main computations ##############################################

  x_internal                           <- x
  levels_decision                      <-
    sort(which(c("Continue to stage 2", "Do not reject", "Reject") %in%
                 levels(x_internal$terminal$decision)))
  if (any(x_internal$des$J == 1, x_internal$des$type != "fisher",
          all(x_internal$des$type == "fisher", x_internal$k == 1))) {
    x_internal$terminal$k              <- factor(x_internal$terminal$k,
                                                 labels = paste("italic(k) ==",
                                                                x_internal$k))
    if (x_internal$des$type == "fisher") {
      x_internal$terminal$xC           <- x_internal$terminal$xC1
      x_internal$terminal$xE           <- x_internal$terminal$xE1
    }
    plot                               <-
      ggplot2::ggplot(x_internal$terminal,
                      ggplot2::aes(x      = .data$xC,
                                   y      = .data$xE,
                                   colour = .data$decision,
                                   shape  = .data$decision)) +
      ggplot2::geom_point(size = 0.75) +
      ggplot2::facet_grid(.~k,
                          labeller = ggplot2::label_parsed) +
      ggplot2::xlab(expression(italic(x[C]))) +
      ggplot2::ylab(expression(italic(x[E]))) +
      ggplot2::scale_colour_manual(values = c("gray55", "firebrick2",
                                              "forestgreen")[levels_decision]) +
      ggplot2::scale_shape_manual(values = c(16, 4, 3)[levels_decision]) +
      ggplot2::coord_fixed(ratio = 1) +
      theme_ph2rand()
  } else {
    which_2                            <- which(x_internal$terminal$k == 2)
    x_internal$terminal$xC             <- x_internal$terminal$xC1
    x_internal$terminal$xC[which_2]    <- x_internal$terminal$xC1[which_2] +
                                            x_internal$terminal$xC2[which_2]
    x_internal$terminal$xE             <- x_internal$terminal$xE1
    x_internal$terminal$xE[which_2]    <- x_internal$terminal$xE1[which_2] +
                                            x_internal$terminal$xE2[which_2]
    x_internal$terminal$state          <- paste("k =", x_internal$terminal$k)
    x_internal$terminal$state[which_2] <-
      paste("k = ", x_internal$terminal$k[which_2], ": z1 = ",
            x_internal$terminal$z1[which_2], sep = "")
    levels_state                       <-
      paste("k = 2: z1 =", sort(unique(x_internal$terminal$z1[which_2])))
    if (1 %in% x$k) {
      levels_state                     <- c("k = 1", levels_state)
    }
    x_internal$terminal$state          <- factor(x_internal$terminal$state,
                                                 levels = levels_state)
    plot                               <-
      ggplot2::ggplot(x_internal$terminal,
                      ggplot2::aes(x      = .data$xC,
                                   y      = .data$xE,
                                   colour = .data$decision,
                                   shape  = .data$decision)) +
      ggplot2::geom_point(size = 0.75) +
      ggplot2::xlab(expression(italic(x[C]))) +
      ggplot2::ylab(expression(italic(x[E]))) +
      ggplot2::scale_colour_manual(values = c("gray55", "firebrick2",
                                              "forestgreen")[levels_decision]) +
      ggplot2::scale_shape_manual(values = c(16, 4, 3)[levels_decision]) +
      ggplot2::coord_fixed(ratio = 1) +
      theme_ph2rand() +
      ggplot2::labs(title = '{current_frame}') +
      gganimate::transition_manual(.data$state) +
      gganimate::ease_aes('linear')
    plot                               <-
      gganimate::animate(plot, fps = 2,
                         renderer = gganimate::gifski_renderer(loop = T))
  }
  print(plot)

  ##### Outputting #############################################################

  if (output) {
    return(list(k       = k,
                output  = output,
                plot    = plot,
                summary = summary,
                x       = x))
  }

}
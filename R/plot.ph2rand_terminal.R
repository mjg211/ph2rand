#' Plot the terminal points of a two-arm randomised clinical trial design for a
#' binary primary outcome variable
#'
#' \code{plot.ph2rand_terminal} plots the terminal points of a design returned
#' by \code{\link{terminal}}.
#' 
#' @param x An object of class \code{ph2rand_terminal}, as returned by
#' \code{\link{terminal}}.
#' @param output A \code{\link{logical}} variable indicating whether outputs
#' should be returned by the function.
#' @return If \code{output = T}, a \code{\link{list}} containing each of the
#' input parameters along with a plot in the slot \code{$plot}, which gives the
#' produced plot of the terminal points.
#' @examples
#' # The default two-stage design
#' des  <- des_two_stage()
#' # Its terminal points across stages 1 and 2
#' term <- terminal(des)
#' # The plot of them
#' plot(term)
#' # Its terminal points from stage 2 only
#' term <- terminal(des, 2)
#' # The plot of them
#' plot(term)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{terminal}}, \code{\link{plot.ph2rand_des}}.
#' @export
plot.ph2rand_terminal <- function(x, output = F, ...) {

  ##### Check inputs ###########################################################

  check_ph2rand_terminal(x)
  check_logical(output, "output")

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
    return(list(output = output,
                plot   = plot,
                x      = x))
  }

}
#' Plot the operating characteristics of a randomised clinical trial design that
#' assumes a Bernoulli distributed primary outcome variable
#'
#' \code{plot.ph2rand_des} plots the operating characteristics of a design
#' returned by \code{\link{des_one_stage}} or \code{\link{des_two_stage}}, under
#' a range of key response rate scenarios. For convenience, it also calls
#' \code{\link{plot.ph2rand_terminal}} to plot the terminal points of the
#' design.
#' 
#' @param x An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' @param k A \code{\link{numeric}} \code{\link{vector}} indicating which stages
#' to consider in determining the probability mass function. That is, it will
#' condition the calculations on the trial ending in the stages given in
#' \code{k}. Defaults to \code{1:des$J} (i.e., to all stages of the given
#' design).
#' @param output A \code{\link{logical}} variable indicating whether available
#' outputs should be returned by the function.
#' @param ... Not currently used.
#' @return If \code{output = TRUE}, a \code{\link{list}} containing each of the
#' input parameters along with a \code{\link{list}} in the slot \code{$plots},
#' which gives all of the available produced plots.
#' @examples
#' \donttest{
#' # The default two-stage design
#' des   <- des_two_stage()
#' # Print several key plots
#' plot(des)
#' # Determine and store all available plots
#' plots <- plot(des, output = TRUE)
#' }
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{plot.ph2rand_terminal}}.
#' @method plot ph2rand_des
#' @export
plot.ph2rand_des <- function(x, k = 1:x$J, output = FALSE, ...) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(x, "any", "x")
  check_k(k, x)
  check_logical(output, "output")

  ##### Perform main computations ##############################################

  plots                  <- list()
  terminal               <- terminal(x, k)
  plots$terminal         <- plot.ph2rand_terminal(terminal, output = TRUE)$plot
  seq_grid               <- seq(0, 1, 0.05)
  opchar_grid            <- opchar(x, as.matrix(expand.grid(seq_grid,
                                                            seq_grid)),
                                   k)$opchar
  plots$grid             <- list()
  plots$grid$`P(pi)`     <-
    ggplot2::ggplot(opchar_grid,
                    ggplot2::aes(x    = .data$piC,
                                 y    = .data$piE,
                                 z    = .data$`P(pi)`,
                                 fill = .data$`P(pi)`)) +
    ggplot2::xlab(expression(italic(pi[C]))) +
    ggplot2::ylab(expression(italic(pi[E]))) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::geom_raster() +
    theme_ph2rand() +
    ggplot2::labs(fill = expression(paste(italic(P), "(", pi, ")", sep = ""))) +
    ggplot2::theme(legend.position = "right",
                   legend.title    = ggplot2::element_text()) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::geom_contour(breaks = x$alpha, linetype = 2,
                          colour = "firebrick") +
    ggplot2::geom_contour(breaks = 1 - x$beta, linetype = 2,
                          colour = "forestgreen")
  x_internal             <- x
  if (length(x$Pi0) > 1) {
    plots$grid$`P(pi)`   <- plots$grid$`P(pi)` +
      ggplot2::geom_segment(x = x_internal$Pi0[1], y = x_internal$Pi0[1],
                            xend = x_internal$Pi0[2], yend = x_internal$Pi0[2],
                            colour = "firebrick") +
      ggplot2::geom_point(ggplot2::aes(x = x_internal$Pi0[1],
                                       y = x_internal$Pi0[1]),
                          colour = "firebrick") +
      ggplot2::geom_point(ggplot2::aes(x = x_internal$Pi0[2],
                                       y = x_internal$Pi0[2]),
                          colour = "firebrick")
  } else {
    plots$grid$`P(pi)`   <- plots$grid$`P(pi)` +
      ggplot2::geom_point(ggplot2::aes(x = x_internal$Pi0[1],
                                       y = x_internal$Pi0[1]),
                          colour = "firebrick")
  }
  if (length(x$Pi1) > 1) {
    plots$grid$`P(pi)`   <- plots$grid$`P(pi)` +
      ggplot2::geom_segment(x = x_internal$Pi1[1],
                            y = x_internal$Pi1[1] + x_internal$delta,
                            xend = x_internal$Pi0[2],
                            yend = x_internal$Pi0[2] + x_internal$delta,
                            colour = "forestgreen") +
      ggplot2::geom_point(ggplot2::aes(x = x_internal$Pi1[1],
                              y = x_internal$Pi1[1] + x_internal$delta),
                          colour = "forestgreen") +
      ggplot2::geom_point(ggplot2::aes(x = x_internal$Pi1[2],
                                       y = x_internal$Pi1[2] +
                                         x_internal$delta),
                          colour = "forestgreen")
  } else {
    plots$grid$`P(pi)`   <- plots$grid$`P(pi)` +
      ggplot2::geom_point(ggplot2::aes(x = x_internal$Pi1[1],
                              y = x_internal$Pi1[1] + x_internal$delta),
                          colour = "forestgreen")
  }
  if (x$J > 1) {
    plots$grid$`ESS(pi)` <-
      ggplot2::ggplot(opchar_grid,
                      ggplot2::aes(x    = .data$piC,
                                   y    = .data$piE,
                                   fill = .data$`ESS(pi)`)) +
      ggplot2::xlab(expression(italic(pi[C]))) +
      ggplot2::ylab(expression(italic(pi[E]))) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::geom_raster() +
      theme_ph2rand() +
      ggplot2::labs(fill = expression(paste(italic(ESS), "(", pi, ")",
                                            sep = ""))) +
      ggplot2::theme(legend.position = "right",
                     legend.title    = ggplot2::element_text()) +
      ggplot2::scale_fill_viridis_c()
    opchar_grid_rej      <- tidyr::gather(opchar_grid, "key", "Probability",
                                          7:10)
    opchar_grid_rej$key  <-
      factor(opchar_grid_rej$key,
             labels = paste0(rep(c(paste0(expression(italic(E)), "["),
                                   paste0(expression(italic(F)), "[")),
                                 each = 2), rep(1:2, 2), "]"))
    plots$grid$rejection <-
      ggplot2::ggplot(opchar_grid_rej,
                      ggplot2::aes(x    = .data$piC,
                                   y    = .data$piE,
                                   fill = .data$Probability)) +
      ggplot2::geom_raster() +
      ggplot2::facet_wrap(~key,
                          labeller = ggplot2::label_parsed) +
      ggplot2::xlab(expression(italic(pi[C]))) +
      ggplot2::ylab(expression(italic(pi[E]))) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::scale_fill_viridis_c() +
      theme_ph2rand() +
      ggplot2::theme(legend.position = "right",
                     legend.title    = ggplot2::element_text())
    opchar_grid_stop     <- tidyr::gather(opchar_grid, "key", "Probability",
                                          11:12)
    opchar_grid_stop$key <- factor(opchar_grid_stop$key,
                                   labels = c("italic(S)[1]", "italic(S)[2]"))
    plots$grid$stopping  <-
      ggplot2::ggplot(opchar_grid_stop,
                      ggplot2::aes(x    = .data$piC,
                                   y    = .data$piE,
                                   z    = .data$Probability,
                                   fill = .data$Probability)) +
      ggplot2::geom_raster() +
      ggplot2::facet_wrap(~key,
                          labeller = ggplot2::label_parsed) +
      ggplot2::xlab(expression(italic(pi[C]))) +
      ggplot2::ylab(expression(italic(pi[E]))) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::scale_fill_viridis_c() +
      theme_ph2rand() +
      ggplot2::theme(legend.position = "right",
                     legend.title    = ggplot2::element_text()) +
      ggplot2::geom_contour(breaks = 0.5, colour = "white", linetype = 2)
  }
  seq_null               <- seq(0, 1, 0.005)
  opchar_null            <- opchar(x, cbind(seq_null, seq_null), k)$opchar
  plots$null$`P(pi)`     <-
    ggplot2::ggplot(opchar_null,
                    ggplot2::aes(x = .data$piC,
                                 y = .data$`P(pi)`)) +
    ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
    ggplot2::ylab(expression(paste(italic(P), "(", pi, ")", sep = ""))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = x$alpha,
                        linetype   = 2) +
    ggplot2::geom_vline(xintercept = x$Pi0,
                        linetype   = 2) +
    theme_ph2rand()
  if (x$J > 1) {
    plots$null$`ESS(pi)` <-
      ggplot2::ggplot(opchar_null,
                      ggplot2::aes(x = .data$piC,
                                   y = .data$`ESS(pi)`)) +
      ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
      ggplot2::ylab(expression(paste(italic(ESS), "(", pi, ")", sep = ""))) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = x$Pi0,
                          linetype   = 2) +
      theme_ph2rand()
    opchar_null_rej      <- tidyr::gather(opchar_null, "key", "value", 7:10)
    opchar_null_rej      <-
      dplyr::mutate(opchar_null_rej,
                    key2 = rep(paste0(rep(c(paste0(expression(italic(E)), "["),
                                            paste0(expression(italic(F)), "[")),
                                          each = 2), rep(1:2, 2), "]"),
                               each = nrow(opchar_null)))
    plots$null$rejection <-
      ggplot2::ggplot(opchar_null_rej,
                      ggplot2::aes(x    = .data$piC,
                                   y    = .data$value,
                                   fill = .data$key2)) +
      ggplot2::geom_area() +
      ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
      ggplot2::ylab("Probability") +
      theme_ph2rand() +
      ggplot2::scale_fill_viridis_d(labels = scales::parse_format()) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(0, 1)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
    opchar_null_stop     <- tidyr::gather(opchar_null, "key", "value", 11:12)
    opchar_null_stop     <-
      dplyr::mutate(opchar_null_stop,
                    key2 = rep(paste0(rep(paste0(expression(italic(S)), "["), 2),
                                      1:2, "]"), each = nrow(opchar_null)))
    plots$null$stopping  <-
      ggplot2::ggplot(opchar_null_stop,
                      ggplot2::aes(x    = .data$piC,
                                   y    = .data$value,
                                   fill = .data$key2)) +
      ggplot2::geom_area() +
      ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
      ggplot2::ylab("Probability") +
      theme_ph2rand() +
      ggplot2::scale_fill_viridis_d(labels = scales::parse_format()) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(0, 1)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  }
  seq_alt               <- seq(0, 1 - x$delta, 0.005)
  opchar_alt            <- opchar(x, cbind(seq_alt, seq_alt + x$delta),
                                  k)$opchar
  plots$alt$`P(pi)`     <-
    ggplot2::ggplot(opchar_alt,
                    ggplot2::aes(x = .data$piC,
                                 y = .data$`P(pi)`)) +
    ggplot2::xlab(bquote(paste(italic(pi[C]), " = ", italic(pi[E]), " - ",
                               delta, " = ", italic(pi[E]), " - ",
                               .(x$delta)))) +
    ggplot2::ylab(expression(paste(italic(P), "(", pi, ")", sep = ""))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1 - x$beta,
                        linetype   = 2) +
    ggplot2::geom_vline(xintercept = x$Pi1,
                        linetype   = 2) +
    theme_ph2rand()
  if (x$J > 1) {
    plots$alt$`ESS(pi)` <-
      ggplot2::ggplot(opchar_alt,
                      ggplot2::aes(x = .data$piC,
                                   y = .data$`ESS(pi)`)) +
      ggplot2::xlab(bquote(paste(italic(pi[C]), " = ", italic(pi[E]), " - ",
                                 delta, " = ", italic(pi[E]), " - ",
                                 .(x$delta)))) +
      ggplot2::ylab(expression(paste(italic(ESS), "(", pi, ")", sep = ""))) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = x$Pi1,
                          linetype   = 2) +
      theme_ph2rand()
    opchar_alt_rej      <- tidyr::gather(opchar_alt, "key", "value", 7:10)
    opchar_alt_rej      <-
      dplyr::mutate(opchar_alt_rej,
                    key2 = rep(paste0(rep(c(paste0(expression(italic(E)), "["),
                                            paste0(expression(italic(F)), "[")),
                                          each = 2), rep(1:2, 2), "]"),
                               each = nrow(opchar_alt)))
    plots$alt$rejection <-
      ggplot2::ggplot(opchar_alt_rej,
                      ggplot2::aes(x    = .data$piC,
                                   y    = .data$value,
                                   fill = .data$key2)) +
      ggplot2::geom_area() +
      ggplot2::xlab(bquote(paste(italic(pi[C]), " = ", italic(pi[E]), " - ",
                                 delta, " = ", italic(pi[E]), " - ",
                                 .(x$delta)))) +
      ggplot2::ylab("Probability") +
      theme_ph2rand() +
      ggplot2::scale_fill_viridis_d(labels = scales::parse_format()) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(0, 1 - x$delta)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
    opchar_alt_stop     <- tidyr::gather(opchar_alt, "key", "value", 11:12)
    opchar_alt_stop     <-
      dplyr::mutate(opchar_alt_stop,
                    key2 = rep(paste0(rep(paste0(expression(italic(S)), "["), 2),
                                      1:2, "]"), each = nrow(opchar_alt)))
    plots$alt$stopping   <- ggplot2::ggplot(opchar_alt_stop,
                                            ggplot2::aes(x    = .data$piC,
                                                         y    = .data$value,
                                                         fill = .data$key2)) +
      ggplot2::geom_area() +
      ggplot2::xlab(bquote(paste(italic(pi[C]), " = ", italic(pi[E]), " - ",
                                 delta, " = ", italic(pi[E]), " - ",
                                 .(x$delta)))) +
      ggplot2::ylab("Probability") +
      theme_ph2rand() +
      ggplot2::scale_fill_viridis_d(labels = scales::parse_format()) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  limits = c(0, 1 - x$delta)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
    opchar              <- rbind(opchar_grid, opchar_null, opchar_alt)
    opchar              <- dplyr::arrange(opchar[!duplicated(opchar), ],
                                          .data$piC, .data$piE)
  }
  if (!output) {
    print(plots$null$`P(pi)`)
    print(plots$alt$`P(pi)`)
  }

  ##### Outputting #############################################################

  if (output) {
    return(list(k      = k,
                opchar = opchar,
                output = output,
                plots  = plots,
                x      = x))
  }

}
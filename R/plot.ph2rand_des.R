#' @export
plot.ph2rand_des <- function(x, k, output = F, summary = F, ...) {

  ##### Check inputs ###########################################################

  check_des(x, "any", "x")
  k  <- check_k(k, x)
  check_logical(output, "output")
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    #summary_plot_des(x, k)
  }

  ##### Perform main computations ##############################################

  plots                <- list()
  terminal             <- terminal(x, k)
  plots$terminal       <- plot(terminal, output = T)$plot
  seq_grid             <- seq(0, 1, 0.05)
  opchar_grid          <- opchar(x, as.matrix(expand.grid(seq_grid, seq_grid)),
                                 k)$opchar
  plots$grid           <- list()
  plots$grid$`P(pi)`   <-
    ggplot2::ggplot(opchar_grid,
                    ggplot2::aes(x    = .data$piC,
                                 y    = .data$piE,
                                 fill = .data$`P(pi)`)) +
    ggplot2::xlab(expression(italic(pi[C]))) +
    ggplot2::ylab(expression(italic(pi[E]))) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::geom_raster() +
    theme_ph2rand() +
    ggplot2::labs(fill = expression(paste(italic(P), "(", pi, ")", sep = ""))) +
    ggplot2::theme(legend.position = "right",
                   legend.title    = ggplot2::element_text()) +
    ggplot2::scale_fill_viridis_c()
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
  opchar_grid_rej      <- tidyr::gather(opchar_grid, "key", "Probability", 7:10)
  opchar_grid_rej$key  <-
    factor(opchar_grid_rej$key,
           labels = paste0(rep(c(paste0(expression(italic(E)), "["),
                                 paste0(expression(italic(F)), "[")), each = 2),
                           rep(1:2, 2), "]"))
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
  seq_null             <- seq(0, 1, 0.005)
  opchar_null          <- opchar(x, cbind(seq_null, seq_null), k)$opchar
  plots$null$`P(pi)`   <-
    ggplot2::ggplot(opchar_null,
                    ggplot2::aes(x = piC,
                                 y = `P(pi)`)) +
    ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
    ggplot2::ylab(expression(paste(italic(P), "(", pi, ")", sep = ""))) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = x$alpha,
                        linetype   = 2) +
    ggplot2::geom_vline(xintercept = x$Pi0,
                        linetype   = 2) +
    theme_ph2rand()
  plots$null$`ESS(pi)` <-
    ggplot2::ggplot(opchar_null,
                    ggplot2::aes(x = piC,
                                 y = `ESS(pi)`)) +
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
                    ggplot2::aes(x    = piC,
                                 y    = value,
                                 fill = key2)) +
    ggplot2::geom_area() +
    ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
    ggplot2::ylab("Probability") +
    theme_ph2rand() +
    ggplot2::scale_fill_viridis_d(labels = scales::parse_format()) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                limits = c(0, 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  opchar_null_stop      <- tidyr::gather(opchar_null, "key", "value", 11:12)
  opchar_null_stop      <-
    dplyr::mutate(opchar_null_stop,
                  key2 = rep(paste0(rep(paste0(expression(italic(S)), "["), 2),
                                    1:2, "]"), each = nrow(opchar_null)))
  plots$null$stopping   <-
    ggplot2::ggplot(opchar_null_stop,
                    ggplot2::aes(x    = piC,
                                 y    = value,
                                 fill = key2)) +
    ggplot2::geom_area() +
    ggplot2::xlab(expression(paste(italic(pi[C]), " = ", italic(pi[E])))) +
    ggplot2::ylab("Probability") +
    theme_ph2rand() +
    ggplot2::scale_fill_viridis_d(labels = scales::parse_format()) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                limits = c(0, 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  seq_alt              <- seq(0, 1 - x$delta, 0.005)
  opchar_alt           <- opchar(x, cbind(seq_alt, seq_alt + x$delta), k)$opchar
  plots$alt$`P(pi)`    <-
    ggplot2::ggplot(opchar_alt,
                    ggplot2::aes(x = piC,
                                 y = `P(pi)`)) +
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
  plots$alt$`ESS(pi)`  <-
    ggplot2::ggplot(opchar_alt,
                    ggplot2::aes(x = piC,
                                 y = `ESS(pi)`)) +
    ggplot2::xlab(bquote(paste(italic(pi[C]), " = ", italic(pi[E]), " - ",
                               delta, " = ", italic(pi[E]), " - ",
                               .(x$delta)))) +
    ggplot2::ylab(expression(paste(italic(ESS), "(", pi, ")", sep = ""))) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = x$Pi1,
                        linetype   = 2) +
    theme_ph2rand()
  opchar_alt_rej       <- tidyr::gather(opchar_alt, "key", "value", 7:10)
  opchar_alt_rej       <-
    dplyr::mutate(opchar_alt_rej,
                  key2 = rep(paste0(rep(c(paste0(expression(italic(E)), "["),
                                          paste0(expression(italic(F)), "[")),
                                        each = 2), rep(1:2, 2), "]"),
                             each = nrow(opchar_alt)))
  plots$alt$rejection  <-
    ggplot2::ggplot(opchar_alt_rej,
                    ggplot2::aes(x    = piC,
                                 y    = value,
                                 fill = key2)) +
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
  opchar_alt_stop      <- tidyr::gather(opchar_alt, "key", "value", 11:12)
  opchar_alt_stop      <-
    dplyr::mutate(opchar_alt_stop,
                  key2 = rep(paste0(rep(paste0(expression(italic(S)), "["), 2),
                                    1:2, "]"), each = nrow(opchar_alt)))
  plots$alt$stopping   <- ggplot2::ggplot(opchar_alt_stop,
                                          ggplot2::aes(x    = piC,
                                                       y    = value,
                                                       fill = key2)) +
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
  opchar               <- rbind(opchar_grid, opchar_null, opchar_alt)
  opchar               <- dplyr::arrange(opchar[!duplicated(opchar), ], piC,
                                         piE)
  print(plots$null$`P(pi)`)
  print(plots$alt$`P(pi)`)

  ##### Outputting #############################################################

  if (output) {
    return(list(k       = k,
                opchar  = opchar,
                output  = output,
                plots   = plots,
                summary = summary,
                x       = x))
  }

}
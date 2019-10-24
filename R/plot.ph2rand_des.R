plot.ph2rand_des <- function(x, k, ...) {
  opchar_grid <- opchar(x, as.matrix(expand.grid(seq(0, 1, 0.05),
                                                 seq(0, 1, 0.05))), k)
  opchar_null <- opchar(x, cbind(seq(0, 1, 0.01), seq(0, 1, 0.01)), k)
  opchar_alt  <- opchar(x, cbind(seq(0, 1 - x$delta, 0.01),
                                 seq(x$delta, 1, 0.01)), k)
  
}

plot.ph2rand_terminal <- function(x, ...) {
  
  internal_x              <- x
  if (any(internal_x$des$J == 1, internal_x$des$type != "fisher",
          all(internal_x$des$type == "fisher", internal_x$k == 1))) {
    internal_x$terminal$k <- factor(internal_x$terminal$k,
                                    labels = paste("italic(k) ==",
                                                   internal_x$k))
    levels_decision       <- sort(which(c("Continue to stage 2",
                                          "Do not reject", "Reject") %in%
                                          levels(internal_x$terminal$decision)))
    if (internal_x$des$type == "fisher") {
      internal_x$terminal$xC <- internal_x$terminal$xC1
      internal_x$terminal$xE <- internal_x$terminal$xE1
    }
    plot                  <-
      ggplot2::ggplot(internal_x$terminal,
                      ggplot2::aes(x      = .data$xC,
                                   y      = .data$xE,
                                   colour = .data$decision,
                                   shape  = .data$decision)) +
      ggplot2::geom_point(size = 0.75) +
      ggplot2::facet_grid(.~k, labeller = ggplot2::label_parsed) +
      ggplot2::xlab(expression(italic(x[C]))) +
      ggplot2::ylab(expression(italic(x[E]))) +
      ggplot2::scale_colour_manual(values = c("gray55", "firebrick2",
                                              "forestgreen")[levels_decision]) +
      ggplot2::scale_shape_manual(values = c(16, 4, 3)[levels_decision]) +
      ggplot2::coord_fixed(ratio = 1) +
      ph2rand:::theme_ph2rand()
  } else {
    levels_decision       <- sort(which(c("Continue to stage 2",
                                          "Do not reject", "Reject") %in%
                                          levels(internal_x$terminal$decision)))
    internal_x$terminal$xC <- internal_x$terminal$xC1
    internal_x$terminal$xC[which(internal_x$terminal$k == 2)] <- internal_x$terminal$xC1[which(internal_x$terminal$k == 2)] + internal_x$terminal$xC2[which(internal_x$terminal$k == 2)]
    internal_x$terminal$xE <- internal_x$terminal$xE1
    internal_x$terminal$xE[which(internal_x$terminal$k == 2)] <- internal_x$terminal$xE1[which(internal_x$terminal$k == 2)] + internal_x$terminal$xE2[which(internal_x$terminal$k == 2)]
    internal_x$terminal$state <- paste("k =", internal_x$terminal$k)
    internal_x$terminal$state[which(internal_x$terminal$k == 2)] <-
      paste("k = ", internal_x$terminal$k[which(internal_x$terminal$k == 2)],
            ": z1 = ", internal_x$terminal$z1[which(internal_x$terminal$k == 2)],
            sep = "")
    levels_state <- paste("k = 2: z1 =", sort(unique(internal_x$terminal$z1[which(internal_x$terminal$k == 2)])))
    if (1 %in% x$k) {
      levels_state <- c("k = 1", levels_state)
    }
    internal_x$terminal$state <- factor(internal_x$terminal$state,
                                        levels = levels_state)
    plot <- ggplot(internal_x$terminal, aes(xC, xE, colour = decision, shape = decision)) +
      geom_point(size = 0.75) +
      ggplot2::xlab(expression(italic(x[C]))) +
      ggplot2::ylab(expression(italic(x[E]))) +
      ggplot2::scale_colour_manual(values = c("gray55", "firebrick2",
                                              "forestgreen")[levels_decision]) +
      ggplot2::scale_shape_manual(values = c(16, 4, 3)[levels_decision]) +
      ggplot2::coord_fixed(ratio = 1) +
      ph2rand:::theme_ph2rand() +
      #facet_wrap(.~state) +
      # Here comes the gganimate specific bits
      labs(title = '{current_frame}') +
      transition_manual(state) +
      ease_aes('linear')
    animate(plot, fps = 2, renderer = gifski_renderer(loop = T))
    
    
    
    
  }
}
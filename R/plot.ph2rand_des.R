plot.ph2rand_des <- function(x, ...) {
  
}

plot.ph2rand_terminal <- function(x, ...) {
  
  if (x$des$J == 1) {
    
  } else {
    if (x$des$type != "fisher") {
      plot <- ggplot2::ggplot(x$terminal,
                              ggplot2::aes(x = xC,
                                           y = xE,
                                           colour = decision)) +
        ggplot2::geom_point() +
        ggplot2::facet_grid(.~k)
    }
  }
}
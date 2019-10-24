plot.ph2rand_pmf <- function(x, output = F, summary = F, ...) {
  
  ##### Check inputs ###########################################################
  
  #check_pmf(x, "any")
  check_logical(output, "output")
  check_logical(summary, "summary")
  
  ##### Print summary ##########################################################
  
  if (summary) {
    #summary_plot_pmf(x)
  }
  
  ##### Perform main computations ##############################################
  
  x_internal                   <- x
  if (all(x$des$type == "fisher", x$des$J == 2)) {
    which_2                    <- which(x_internal$pmf$k == 2)
    x_internal$pmf$xC          <- x_internal$pmf$xC1
    x_internal$pmf$xC[which_2] <- x_internal$pmf$xC1[which_2] +
                                    x_internal$pmf$xC2[which_2]
    x_internal$pmf$xE          <- x_internal$pmf$xE1
    x_internal$pmf$xE[which_2] <- x_internal$pmf$xE1[which_2] +
                                    x_internal$pmf$xE2[which_2]
  }
  x_internal$pmf$k             <- as.numeric(x_internal$pmf$k)
  plots                        <- list()
  counter                      <- 1L
  unique_piC                   <- unique(x_internal$pmf$piC)
  for (i in length(unique_piC)) {
    pmf_i                      <- dplyr::filter(x_internal$pmf,
                                                piC == unique_piC[i])
    unique_piE                 <- unique(pmf_i$piE)
    for (j in 1:length(unique_piE)) {
      pmf_ij                   <- dplyr::filter(pmf_i, piE == unique_piE[j])
      pmf_ij                   <-
        dplyr::summarise(dplyr::group_by(pmf_ij, xC, xE, mC, mE),
                         f = sum(`f(x,m|pi)`),
                         k = sum(k)/dplyr::n())
      pmf_ij$k                 <- factor(pmf_ij$k,
                                         labels = paste("italic(k) ==",
                                                        x_internal$k))
      plots[[counter]]         <-
        ggplot2::ggplot(pmf_ij,
                        ggplot2::aes(x    = xC,
                                     y    = xE,
                                     fill = f)) +
        ggplot2::facet_grid(.~k,
                            labeller = ggplot2::label_parsed) +
        ggplot2::xlab(expression(italic(x[C]))) +
        ggplot2::ylab(expression(italic(x[E]))) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::geom_raster() +
        theme_ph2rand() +
        ggplot2::labs(fill  = expression(paste(italic(f), "(", italic(x), ",",
                                               italic(m), "|", pi, ")",
                                               sep = "")),
                      title =
                        bquote(paste(pi[italic(C)], " = ", .(unique_piC[i]),
                                     ", ", pi[italic(E)], " = ",
                                     .(unique_piE[j]), sep = ""))) +
        ggplot2::theme(legend.position = "right",
                       legend.title    = ggplot2::element_text(),
                       plot.title      = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_fill_viridis_c()
      counter                  <- counter + 1L
    }
  }
  print(plots[[1]])
  
  ##### Outputting #############################################################
  
  if (output) {
    return(list(output  = output,
                plots   = plots,
                summary = summary,
                x       = x))
  }
  
}
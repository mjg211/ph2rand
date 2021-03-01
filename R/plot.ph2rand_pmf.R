#' Plot probability mass functions of a randomised clinical trial design that
#' assumes a Bernoulli distributed primary outcome variable
#'
#' \code{plot.ph2rand_pmf} plots the terminal points of a design returned by
#' \code{\link{pmf}}.
#' 
#' @param x An object of class \code{ph2rand_pmf}, as returned by
#' \code{\link{pmf}}.
#' @param output A \code{\link{logical}} variable indicating whether outputs
#' should be returned by the function.
#' @param ... Not currently used.
#' @return If \code{output = TRUE}, a \code{\link{list}} containing each of the
#' input parameters along with a plot in the slot \code{$plot}, which gives the
#' produced plot of the terminal points.
#' @examples
#' # The default two-stage design
#' des <- des_two_stage()
#' # Its probability mass function under the uninteresting and interesting
#' # scenarios
#' pmf <- pmf(des)
#' # The plot of them
#' plot(pmf)
#' # The same probability mass functions, conditioning on the trial ending in
#' # stage 2
#' pmf <- pmf(des, k = 2)
#' # The plot of them
#' plot(pmf)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{pmf}}, \code{\link{plot.ph2rand_des}}.
#' @method plot ph2rand_pmf
#' @export
plot.ph2rand_pmf <- function(x, output = FALSE, ...) {
  
  ##### Check inputs ###########################################################
  
  check_ph2rand_pmf(x)
  check_logical(output, "output")
  
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
  for (i in 1:length(unique_piC)) {
    pmf_i                      <- dplyr::filter(x_internal$pmf,
                                                .data$piC == unique_piC[i])
    unique_piE                 <- unique(pmf_i$piE)
    for (j in 1:length(unique_piE)) {
      pmf_ij                   <- dplyr::filter(pmf_i,
                                                .data$piE == unique_piE[j])
      if (nrow(pmf_ij) != 0) {
        pmf_ij                 <-
          dplyr::summarise(dplyr::group_by(pmf_ij, .data$xC, .data$xE, .data$mC,
                                           .data$mE),
                           f = sum(.data$`f(x,m|pi)`),
                           k = sum(.data$k)/dplyr::n())
        pmf_ij$k               <- factor(pmf_ij$k,
                                         labels = paste("italic(k) ==",
                                                        x_internal$k))
        if (abs(sum(pmf_ij$f) - 1) > 1e-14) {
          warning("PMF for piC = ", unique_piC[i], " and piE = ", unique_piE[j],
                  " does not sum to 1")
        }
        plots[[counter]]       <-
          ggplot2::ggplot(pmf_ij,
                          ggplot2::aes(x    = .data$xC,
                                       y    = .data$xE,
                                       fill = .data$f)) +
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
        counter                <- counter + 1L
      }
    }
  }
  print(plots[[1]])
  
  ##### Outputting #############################################################
  
  if (output) {
    return(list(output = output,
                plots  = plots,
                x      = x))
  }
  
}
build_des_one_stage_output <- function(alpha, beta, boundaries, delta, feasible,
                                       nC, nCmax, nE, opchar, Pi0, Pi1, ratio,
                                       summary, type) {
  output        <- list(alpha      = alpha,
                        beta       = beta,
                        boundaries = boundaries,
                        delta      = delta,
                        feasible   = feasible,
                        J          = 1L,
                        nC         = nC,
                        nCmax      = nCmax,
                        nE         = nE,
                        opchar     = opchar,
                        Pi0        = Pi0,
                        Pi1        = Pi1,
                        ratio      = ratio,
                        summary    = summary,
                        type       = type)
  class(output) <- c("ph2rand_des", class(output))
  output
}

build_des_two_stage_output <- function(alpha, beta, boundaries, delta, efficacy,
                                       efficacy_param, efficacy_type, equal,
                                       feasible, futility, futility_param,
                                       futility_type, nC, nCmax, nE, opchar,
                                       Pi0, Pi1, piO, ratio, summary, w, type) {
  output        <- list(alpha          = alpha,
                        beta           = beta,
                        boundaries     = boundaries,
                        delta          = delta,
                        efficacy       = efficacy,
                        efficacy_param = efficacy_param,
                        efficacy_type  = efficacy_type,
                        equal          = equal,
                        feasible       = feasible,
                        futility       = futility,
                        futility_param = futility_param,
                        futility_type  = futility_type,
                        J              = 2L,
                        nC             = nC,
                        nCmax          = nCmax,
                        nE             = nE,
                        opchar         = opchar,
                        Pi0            = Pi0,
                        Pi1            = Pi1,
                        piO            = piO,
                        ratio          = ratio,
                        summary        = summary,
                        type           = type,
                        w              = w)
  class(output) <- c("ph2rand_des", class(output))
  output
}

search_parameters          <- function(J, type, nCmax, ratio) {
  poss_nC                                          <- 1:nCmax
  poss_nE                                          <- poss_nE_orig <-
                                                      poss_nC*ratio
  keep                                             <- which(poss_nE%%1 == 0)
  poss_nC                                          <- poss_nC[keep]
  poss_nE                                          <- poss_nE[keep]
  len_poss_nC                                      <- length(poss_nC)
  max_poss_nC                                      <- max(poss_nC)
  if (type == "fisher") {
    maxima                                         <- max(max_poss_nC, poss_nE)
    choose_mat                                     <- matrix(0, maxima,
                                                             maxima + 1)
    for (n in 1:maxima) {
      choose_mat[n, 1:(n + 1)]                     <- choose(n, 0:n)
    }
  } else {
    choose_mat                                     <- NULL
  }
  poss_x <- poss_y <- poss_z <- poss_B <- unique_B <- list()
  all_x                                            <-
    as.matrix(expand.grid(0:poss_nC[len_poss_nC], 0:poss_nE[len_poss_nC]))
  all_y                                            <- all_x[, 2] - all_x[, 1]
  all_z                                            <- all_x[, 1] + all_x[, 2]
  for (n in 1:len_poss_nC) {
    nC                                             <- poss_nC[n]
    nE                                             <- poss_nE[n]
    index                                          <- nC + max_poss_nC*(nE - 1)
    keep                                           <- which(all_x[, 1] <= nC &
                                                              all_x[, 2] <= nE)
    poss_x[[index]]                                <- all_x[keep, ]
    if (type != "barnard") {
      poss_y[[index]]                              <- all_y[keep]
      if (type == "fisher") {
        poss_z[[index]]                            <- all_z[keep]
      }
    } else {
      denom_fact                                   <-
        (poss_x[[index]][, 1] + poss_x[[index]][, 2])/(nC + nE)
      poss_B_index                                 <-
        (poss_x[[index]][, 2]/nE - poss_x[[index]][, 1]/nC)/
        sqrt(denom_fact*(1 - denom_fact)*(1/nC + 1/nE))
      poss_B_index[is.nan(poss_B_index)]           <- 0
      poss_B[[index]]                              <- matrix(0, nC + 1, nE + 1)
      for (i in 1:((nC + 1)*(nE + 1))) {
        poss_B[[index]][poss_x[[index]][i, 1] + 1,
                        poss_x[[index]][i, 2] + 1] <- poss_B_index[i]
      }
      unique_B[[index]]                            <-
        sort(unique(as.vector(poss_B_index)))
      len_unique_B_index                           <- length(unique_B[[index]])
      keep                                         <-
        !logical(len_unique_B_index)
      for (i in 2:len_unique_B_index) {
        if (unique_B[[index]][i] - unique_B[[index]][i - 1] <= 1e-15) {
          keep[i]                                  <- F
        }
      }
      unique_B[[index]]                            <- unique_B[[index]][keep]
      unique_B[[index]]                            <-
        c(unique_B[[index]][1] - 1, unique_B[[index]],
          unique_B[[index]][length(unique_B[[index]])] + 1)
    }
  }
  list(choose_mat   = choose_mat,
       max_poss_nC  = max_poss_nC,
       poss_nC      = poss_nC,
       poss_nE      = poss_nE,
       poss_nE_orig = poss_nE_orig,
       poss_B       = poss_B,
       poss_x       = poss_x,
       poss_y       = poss_y,
       poss_z       = poss_z,
       unique_B     = unique_B)
}

theme_ph2rand              <- function(base_size = 11, base_family = "") {
  ggplot2::theme_grey(base_family = base_family,
                      base_size   = base_size) +
    ggplot2::theme(axis.ticks       = ggplot2::element_line(colour = "grey70",
                                                            size   = 0.25),
                   complete         = T,
                   legend.key       = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   legend.position  = "bottom",
                   legend.title     = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   panel.border     = ggplot2::element_rect(fill   = NA,
                                                            colour = "grey70",
                                                            size   = 0.5),
                   panel.grid.major = ggplot2::element_line(colour = "grey87",
                                                            size   = 0.25),
                   panel.grid.minor = ggplot2::element_line(colour = "grey87",
                                                            size   = 0.125),
                   plot.margin      = ggplot2::unit(c(0.3, 0.5, 0.3, 0.3),
                                                    "cm"),
                   plot.title       = ggplot2::element_text(hjust = 0.5),
                   strip.background = ggplot2::element_rect(fill   = "grey70",
                                                            colour = NA),
                   strip.text       =
                     ggplot2::element_text(colour = "white",
                                           size   = ggplot2::rel(0.8)))
}

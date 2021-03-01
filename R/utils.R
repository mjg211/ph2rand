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
          keep[i]                                  <- FALSE
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

sim_internal               <- function(pi, completed_replicates, des, k,
                                       replicates, summary, total_replicates) {
  summary_i                   <-
    round(seq(1, total_replicates, length.out = 11)[-c(1, 11)])
  J                           <- des$J
  nC                          <- des$nC
  nE                          <- des$nE
  if (des$type %in% c("barnard", "binomial")) {
    e                         <- list(des$boundaries$e1, des$boundaries$e2)
    f                         <- list(des$boundaries$f1, des$boundaries$f2)
  }
  cum_nC                      <- cumsum(nC)
  cum_nE                      <- cumsum(nE)
  seq_J                       <- 1:J
  E                           <- Fu <- numeric(J)
  numeric_2                   <- numeric(2)
  for (i in 1:replicates) {
    x                         <- z <- numeric_2
    for (j in seq_J) {
      x_iterate               <- stats::rbinom(2, c(nC[j], nE[j]), pi)
      z[j]                    <- x_iterate[1] + x_iterate[2]
      x                       <- x + x_iterate
      if (des$type %in% c("binomial", "fisher", "sat")) {
        tD                    <- x[2] - x[1]
        if (des$type == "sat") {
          tS                  <- x[2]
        }
      } else if (des$type == "barnard") {
        if (any(all(x == 0), all(x == c(cum_nC[j], cum_nE[j])))) {
          tB                  <- 0
        } else {
          fact                <- (x[1] + x[2])/(cum_nC[j] + cum_nE[j])
          tB                  <-
            (x[2]/cum_nE[j] - x[1]/cum_nC[j])/
            sqrt(fact*(1 - fact)*(1/cum_nC[j] + 1/cum_nE[j]))
        }
      }
      continue                <- TRUE
      if (des$type == "barnard") {
        if (tB >= e[[j]]) {
          E[j]                <- E[j] + 1
          continue            <- FALSE
        } else if (tB <= f[[j]]) {
          Fu[j]               <- Fu[j] + 1
          continue            <- FALSE
        }
      } else if (des$type == "binomial") {
        if (tD >= e[[j]]) {
          E[j]                <- E[j] + 1
          continue            <- FALSE
        } else if (tD <= f[[j]]) {
          Fu[j]               <- Fu[j] + 1
          continue            <- FALSE
        }
      } else if (des$type == "fisher") {
        if (j == 1) {
          if (tD >= des$boundaries$e1[z[1] + 1]) {
            E[1]              <- E[1] + 1
            continue          <- FALSE
          } else if (tD <= des$boundaries$f1[z[1] + 1]) {
            Fu[1]             <- Fu[1] + 1
            continue          <- FALSE
          }
        } else {
          if (tD >= des$boundaries$e2[z[1] + 1, z[2] + 1]) {
            E[2]              <- E[2] + 1
          } else {
            Fu[2]             <- Fu[2] + 1
          }
        }
      } else if (des$type == "sat") {
        if (j == 1) {
          if (all(tD >= des$boundaries$eT1, tS >= des$boundaries$eS1)) {
            E[1]              <- E[1] + 1
            continue          <- FALSE
          } else if (all(tD <= des$boundaries$fT1, tS <= des$boundaries$fS1)) {
            Fu[1]             <- Fu[1] + 1
            continue          <- FALSE
          }
        } else {
          if (all(tD >= des$boundaries$eT2, tS >= des$boundaries$eS2)) {
            E[2]              <- E[2] + 1
            continue          <- FALSE
          } else {
            Fu[2]             <- Fu[2] + 1
          }
        }
      }
      if (!continue) {
        break
      }
    }
    if (all((completed_replicates + i) %in% summary_i, summary)) {
      message("..approximately ",
              10*which(summary_i == (completed_replicates + i)),
              "% through the required simulations..")
    }
  }
  E          <- E/replicates
  Fu         <- Fu/replicates
  if (J == 1) {
    c(pi, E[1])
  } else {
    S        <- E + Fu
    if (length(k) == 1) {
      E      <- E/S[k]
      Fu     <- Fu/S[k]
      E[-k]  <- 0
      Fu[-k] <- 0
      S      <- E + Fu
    }
    n        <- c(rep(nC[1] + nE[1], replicates*S[1]),
                  rep(cum_nC[2] + cum_nE[2], replicates*S[2]))
    c(pi, sum(E), sum(n)/replicates, stats::sd(n), stats::quantile(n, 0.5), E,
      Fu, S, cum_nC[2] + cum_nE[2])
  }
}

theme_ph2rand              <- function(base_size = 11, base_family = "") {
  ggplot2::theme_grey(base_family = base_family,
                      base_size   = base_size) +
    ggplot2::theme(axis.ticks       = ggplot2::element_line(colour = "grey70",
                                                            size   = 0.25),
                   complete         = TRUE,
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

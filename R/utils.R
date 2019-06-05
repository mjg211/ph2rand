build_des_one_stage_output <- function(alpha, beta, delta, feasible, n0max,
                                       opchar, pi_alt, pi_null, point_alt,
                                       point_null, ratio, summary, type,
                                       type_components) {
  if (type %in% c("bernard", "binomial")) {
    output      <- list(alpha      = alpha,
                        beta       = beta,
                        delta      = delta,
                        e          = type_components$e,
                        feasible   = feasible,
                        n0         = type_components$n0,
                        n0max      = n0max,
                        n1         = type_components$n1,
                        opchar     = opchar,
                        pi_alt     = pi_alt,
                        pi_null    = pi_null,
                        point_alt  = point_alt,
                        point_null = point_null,
                        ratio      = ratio,
                        summary    = summary,
                        type       = type)
  } else if (type == "fisher") {
    output      <- list(alpha      = alpha,
                        beta       = beta,
                        delta      = delta,
                        e          = type_components$e,
                        feasible   = feasible,
                        feasible_e = type_components$feasible_e,
                        n0         = type_components$n0,
                        n0max      = n0max,
                        n1         = type_components$n1,
                        opchar     = opchar,
                        pi_alt     = pi_alt,
                        pi_null    = pi_null,
                        point_alt  = point_alt,
                        point_null = point_null,
                        ratio      = ratio,
                        summary    = summary,
                        type       = type)
  } else if (type == "single_double") {
    output      <- list(alpha      = alpha,
                        beta       = beta,
                        delta      = delta,
                        eS          = type_components$eS,
                        eT          = type_components$eT,
                        feasible   = feasible,
                        n0         = type_components$n0,
                        n0max      = n0max,
                        n1         = type_components$n1,
                        opchar     = opchar,
                        pi_alt     = pi_alt,
                        pi_null    = pi_null,
                        point_alt  = point_alt,
                        point_null = point_null,
                        ratio      = ratio,
                        summary    = summary,
                        type       = type)
  }
  class(output) <- c(class(output), "ph2rand_des_one_stage")
  output
}

build_des_two_stage_output <- function(alpha, beta, delta, equal, feasible,
                                       n0max, opchar, pi_alt, pi_ess, pi_null,
                                       point_alt, point_null, ratio, summary, w,
                                       type, type_components) {
  if (type %in% c("bernard", "binomial")) {
    output      <- list(alpha      = alpha,
                        beta       = beta,
                        delta      = delta,
                        e1         = type_components$e1,
                        e2         = type_components$e2,
                        efficacy   = type_components$efficacy,
                        equal      = equal,
                        f1         = type_components$f1,
                        feasible   = feasible,
                        futility   = type_components$futility,
                        n0         = type_components$n0,
                        n0max      = n0max,
                        n1         = type_components$n1,
                        opchar     = opchar,
                        pi_alt     = pi_alt,
                        pi_ess     = pi_ess,
                        pi_null    = pi_null,
                        point_alt  = point_alt,
                        point_null = point_null,
                        ratio      = ratio,
                        summary    = summary,
                        type       = type,
                        w          = w)
  } else if (type == "fisher") {
    output      <- list(alpha          = alpha,
                        beta           = beta,
                        delta          = delta,
                        e1             = type_components$e1,
                        e2             = type_components$e2,
                        efficacy_param = type_components$efficacy_param,
                        efficacy_type  = type_components$efficacy_type,
                        f1             = type_components$f1,
                        feasible       = feasible,
                        feasible_e1    = type_components$feasible_e1,
                        feasible_e2    = type_components$feasible_e2,
                        feasible_f1    = type_components$feasible_f1,
                        futility_param = type_components$futility_param,
                        futility_type  = type_components$futility_type,
                        n0             = type_components$n0,
                        n0max          = n0max,
                        n1             = type_components$n1,
                        opchar         = opchar,
                        pi_alt         = pi_alt,
                        pi_ess         = pi_ess,
                        pi_null        = pi_null,
                        point_alt      = point_alt,
                        point_null     = point_null,
                        ratio          = ratio,
                        summary        = summary,
                        type           = type,
                        w              = w)
  } else if (type == "single_double") {
    output      <- list(alpha      = alpha,
                        beta       = beta,
                        delta      = delta,
                        efficacy   = type_components$efficacy,
                        equal      = equal,
                        eS1        = type_components$eS1,
                        eS2        = type_components$eS2,
                        eT1        = type_components$eT1,
                        eT2        = type_components$eT2,
                        feasible   = feasible,
                        fS1        = type_components$fS1,
                        fT1        = type_components$fT1,
                        futility   = type_components$futility,
                        n0         = type_components$n0,
                        n0max      = n0max,
                        n1         = type_components$n1,
                        opchar     = opchar,
                        pi_alt     = pi_alt,
                        pi_ess     = pi_ess,
                        pi_null    = pi_null,
                        point_alt  = point_alt,
                        point_null = point_null,
                        ratio      = ratio,
                        summary    = summary,
                        type       = type,
                        w          = w)
  }
  class(output) <- c(class(output), "ph2rand_des_two_stage")
  output
}

search_parameters <- function(J, type, n0max, ratio) {
  poss_n0                                          <- 1:n0max
  poss_n1                                          <- poss_n0*ratio
  keep                                             <- which(poss_n1%%1 == 0)
  poss_n0                                          <- poss_n0[keep]
  poss_n1                                          <- poss_n1[keep]
  if (type == "fisher") {
    maxima                                         <- max(c(poss_n0, poss_n1))
    choose_mat                                     <- matrix(0, maxima,
                                                             maxima + 1)
    for (n in 1:maxima) {
      choose_mat[n, 1:(n + 1)]                     <- choose(n, 0:n)
    }
  } else {
    choose_mat                                     <- NULL
  }
  max_poss_n0                                      <- max(poss_n0)
  poss_x <- poss_y <- poss_z <- poss_B <- unique_B <- list()
  len_poss_n0                                      <- length(poss_n0)
  all_x                                            <-
    as.matrix(expand.grid(0:poss_n0[len_poss_n0], 0:poss_n1[len_poss_n0]))
  all_y                                            <- all_x[, 2] - all_x[, 1]
  all_z                                            <- all_x[, 1] + all_x[, 2]
  for (n in 1:length(poss_n0)) {
    n0                                             <- poss_n0[n]
    n1                                             <- poss_n1[n]
    index                                          <- n0 + max_poss_n0*(n1 - 1)
    keep                                           <-
      which(all_x[, 1] <= n0 & all_x[, 2] <= n1)
    poss_x[[index]]                                <- all_x[keep, ]
      
    if (type != "bernard") {
      poss_y[[index]]                              <- all_y[keep]
      if (type == "fisher") {
        poss_z[[index]]                            <- all_z[keep]
      }
    } else {
      denom_fact                                   <-
        (poss_x[[index]][, 1] + poss_x[[index]][, 2])/(n0 + n1)
      poss_B_index                                 <-
        (poss_x[[index]][, 2]/n1 - poss_x[[index]][, 1]/n0)/
        sqrt(denom_fact*(1 - denom_fact)*(1/n0 + 1/n1))
      poss_B_index[is.nan(poss_B_index)]           <- 0
      poss_B[[index]]                              <- matrix(0, n0 + 1, n1 + 1)
      for (i in 1:((n0 + 1)*(n1 + 1))) {
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
    }
  }
  return(list(choose_mat  = choose_mat,
              max_poss_n0 = max_poss_n0,
              poss_n0     = poss_n0,
              poss_n1     = poss_n1,
              poss_B      = poss_B,
              poss_x      = poss_x,
              poss_y      = poss_y,
              poss_z      = poss_z,
              unique_B    = unique_B))
}

theme_ph2rand     <- function(base_size = 11, base_family = "") {
  ggplot2::theme_grey(base_size   = base_size,
                      base_family = base_family) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   panel.border     = ggplot2::element_rect(fill   = NA,
                                                            colour = "grey70",
                                                            size   = 0.5),
                   panel.grid.major = ggplot2::element_line(colour = "grey87",
                                                            size   = 0.25),
                   panel.grid.minor = ggplot2::element_line(colour = "grey87",
                                                            size   = 0.125),
                   axis.ticks       = ggplot2::element_line(colour = "grey70",
                                                            size   = 0.25),
                   legend.key       = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   strip.background = ggplot2::element_rect(fill   = "grey70",
                                                            colour = NA),
                   strip.text       =
                     ggplot2::element_text(colour = "white",
                                           size   = ggplot2::rel(0.8)),
                   legend.title     = ggplot2::element_blank(),
                   legend.position  = "bottom",
                   plot.margin      = ggplot2::unit(c(0.3, 0.5, 0.3, 0.3),
                                                    "cm"),
                   plot.title       = ggplot2::element_text(hjust = 0.5),
                   complete         = T)
}

row_match         <- function(vec, mat) {
  cvec <- do.call("paste", c(vec[, , drop = FALSE], sep = "\r"))
  cmat <- do.call("paste", c(mat[, , drop = FALSE], sep = "\r"))
  match(cvec, cmat, nomatch = NA_integer_)
}

uc                <- function(char) {
  lookup <- matrix(c("alpha",    "\u03B1",
                     "beta",     "\u03B2",
                     "gamma",    "\u03B3",
                     "delta",    "\u03B4",
                     "epsilon",  "\u03B5",
                     "zeta",     "\u03B6",
                     "eta",      "\u03B7",
                     "theta",    "\u03B8",
                     "iota",     "\u03B9",
                     "kappa",    "\u03BA",
                     "lambda",   "\u03BB",
                     "mu",       "\u03BC",
                     "nu",       "\u03BD",
                     "xi",       "\u03BE",
                     "omicron",  "\u03BF",
                     "pi",       "\u03C0",
                     "rho",      "\u03C1",
                     "sigma",    "\u03C3",
                     "tau",      "\u03C4",
                     "upsilon",  "\u03C5",
                     "phi",      "\u03C6",
                     "chi",      "\u03C7",
                     "psi",      "\u03C8",
                     "omega",    "\u03C9",
                     "Alpha",    "\u0391",
                     "Beta",     "\u0392",
                     "Gamma",    "\u0393",
                     "Delta",    "\u0394",
                     "Epsilon",  "\u0395",
                     "Zeta",     "\u0396",
                     "Eta",      "\u0397",
                     "Theta",    "\u0398",
                     "Iota",     "\u0399",
                     "Kappa",    "\u039A",
                     "Lambda",   "\u039B",
                     "Mu",       "\u039C",
                     "Nu",       "\u039D",
                     "Xi",       "\u039E",
                     "Omicron",  "\u039F",
                     "Pi",       "\u03A0",
                     "Rho",      "\u03A1",
                     "Sigma",    "\u03A3",
                     "Tau",      "\u03A4",
                     "Upsilon",  "\u03A5",
                     "Phi",      "\u03A6",
                     "Chi",      "\u03A7",
                     "Psi",      "\u03A8",
                     "Omega",    "\u03A9",
                     "le",       "\u2264",
                     "third",    "\u2153",
                     "quarter",  "\u00BC",
                     "fifth",    "\u2155",
                     "sixth",    "\u2159",
                     "eigth",    "\u215B",
                     "two_elip", "\u2026\u2026"),
                   ncol = 2, byrow = T)
  lookup[which(lookup[, 1] == char), 2]
}

uc_sub            <- function(n) {
  codes <- c("\u2080", "\u2081", "\u2082", "\u2083", "\u2084", "\u2085",
             "\u2086", "\u2087", "\u2088", "\u2089")
  if (n < 10) {
    codes[n + 1]
  } else {
    paste0(codes[n%/%10 + 1], codes[n%%10 + 1])
  }
}

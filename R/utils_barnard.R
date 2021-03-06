barnard_des_one_stage        <- function(alpha, beta, delta, ratio, Pi0, Pi1,
                                         nCmax, summary) {
  params               <- search_parameters(1, "barnard", nCmax, ratio)
  feasible             <-
    barnard_des_one_stage_cpp(alpha, beta, delta, params$poss_nC,
                              params$poss_nE, params$poss_x, params$poss_B,
                              params$unique_B, (length(Pi0) == 1), Pi0,
                              (length(Pi1) == 1), Pi1, summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message("..feasible designs identified in the range of ",
              "considered sample sizes.\n  Identifying the optimal design..")
    }
    nrow_feasible      <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible         <- matrix(c(1, feasible[, 1],
                                   params$poss_nE_orig[feasible[, 1]],
                                   feasible[, -1]), 1)
    } else {
      feasible         <- cbind(1:nrow_feasible, feasible[, 1],
                                params$poss_nE_orig[feasible[, 1]],
                                feasible[, -1])
    }
    colnames(feasible) <- c("index", "n1C", "n1E", "e1", "argmax alpha",
                            "max alpha", "argmin power", "min power")
    feasible           <- tibble::as_tibble(feasible)
    feasible[, 1:3]    <- dplyr::mutate_if(feasible[, 1:3], is.double,
                                           as.integer)
    feasible           <- dplyr::arrange(feasible, .data$n1C,
                                         dplyr::desc(.data$`min power`))
    e1                 <- feasible$e1[1]
    nC                 <- feasible$n1C[1]
    nE                 <- feasible$n1E[1]
    opchar             <-
      barnard_opchar_one_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                     feasible$`argmin power`[1] + c(0, delta)),
                               nC, nE, e1)
  } else {
    e1                 <- feasible <- nC <- nE <- opchar <- NULL
    if (summary) {
      message("..no feasible designs identified in range of ",
              "considered maximal allowed sample size.\n  Consider increasing ",
              "nCmax..")
    }
  }
  build_des_one_stage_output(alpha, beta, list(e1 = e1, f1 = e1), delta,
                             feasible, nC, nCmax, nE, opchar, Pi0, Pi1, ratio,
                             summary, "barnard")
}

barnard_des_two_stage        <- function(alpha, beta, delta, ratio, Pi0, Pi1,
                                         nCmax, equal, w, piO, efficacy,
                                         futility, summary) {
  params               <- search_parameters(2, "barnard", nCmax, ratio)
  feasible             <-
    barnard_des_two_stage_cpp(alpha, beta, delta, params$poss_nC,
                              params$poss_nE, params$poss_x, params$poss_B,
                              params$unique_B, (length(Pi0) == 1), Pi0,
                              (length(Pi1) == 1), Pi1, equal, efficacy,
                              futility, piO, summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message("..feasible designs identified in the range of ",
              "considered sample sizes.\n  Identifying the optimal design..")
    }
    nrow_feasible      <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible         <- matrix(c(1, feasible[, 1:2],
                                   params$poss_nE_orig[feasible[, 1]],
                                   params$poss_nE_orig[feasible[, 2]],
                                   feasible[, -(1:2)]), 1)
    } else {
      feasible         <- cbind(1:nrow_feasible, feasible[, 1:2],
                                params$poss_nE_orig[feasible[, 1]],
                                params$poss_nE_orig[feasible[, 2]],
                                feasible[, -(1:2)])
    }
    colnames(feasible) <- c("index", "n1C", "n2C", "n1E", "n2E", "e1", "e2",
                            "f1", "argmax alpha", "max alpha", "argmin power",
                            "min power", "ESS(piO,piO)", "ESS(piO,piO+delta)")
    feasible           <- tibble::as_tibble(feasible)
    feasible[, 1:5]    <- dplyr::mutate_if(feasible[, 1:5], is.double,
                                           as.integer)
    feasible           <-
      dplyr::mutate(feasible,
                    `argmax ESS(pi,pi)`       = 0,
                    `max ESS(pi,pi)`          = 0,
                    `argmax_piC ESS(piC,piE)` = 0,
                    `argmax_piE ESS(piC,piE)` = 0,
                    `max ESS(piC,piE)`        = 0,
                    `max N`                   =
                      as.integer(rowSums(feasible[, 2:5])))
    for (i in 1:nrow_feasible) {
      index            <- as.numeric(feasible[i, 2] +
                                       params$max_poss_nC*(feasible[i, 4] - 1))
      max_ESS_1d       <-
        barnard_max_ess_1d_two_stage(as.numeric(feasible[i, 2:3]),
                                     as.numeric(feasible[i, 4:5]),
                                     as.numeric(feasible[i, 6]),
                                     as.numeric(feasible[i, 8]),
                                     params$poss_x[[index]],
                                     params$poss_B[[index]])
      feasible$`argmax ESS(pi,pi)`[i]       <- max_ESS_1d[1]
      feasible$`max ESS(pi,pi)`[i]          <- max_ESS_1d[2]
      max_ESS_2d       <-
        barnard_max_ess_2d_two_stage(as.numeric(feasible[i, 2:3]),
                                     as.numeric(feasible[i, 4:5]),
                                     as.numeric(feasible[i, 6]),
                                     as.numeric(feasible[i, 8]),
                                     params$poss_x[[index]],
                                     params$poss_B[[index]])
      feasible$`argmax_piC ESS(piC,piE)`[i] <- max_ESS_2d[1]
      feasible$`argmax_piE ESS(piC,piE)`[i] <- max_ESS_2d[2]
      feasible$`max ESS(piC,piE)`[i]        <- max_ESS_2d[3]
    }
    feasible$o         <- rowSums(matrix(w, nrow_feasible, 5, TRUE)*
                                    feasible[, c(13:14, 16, 19:20)])
    feasible           <- dplyr::arrange(feasible, .data$o,
                                         dplyr::desc(.data$`min power`))
    if (!efficacy) {
      feasible$e1      <- Inf
    }
    if (!futility) {
      feasible$f1      <- -Inf
    }
    e1                 <- feasible$e1[1]
    e2                 <- feasible$e2[1]
    f1                 <- feasible$f1[1]
    nC                 <- c(feasible$n1C[1], feasible$n2C[1])
    nE                 <- c(feasible$n1E[1], feasible$n2E[1])
    opchar             <-
      barnard_opchar_two_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                     feasible$`argmin power`[1] + c(0, delta)),
                               nC, nE, e1, f1, e2, 1:2)
  } else {
    e1                 <- e2 <- f1 <- feasible <- nC <- nE <- opchar <- NULL
    if (summary) {
      message("..no feasible designs identified in range of ",
              "considered maximal allowed sample size.\n  Consider increasing ",
              "nCmax..")
    }
  }
  build_des_two_stage_output(alpha, beta,
                             list(e1 = e1, e2 = e2, f1 = f1, f2 = e2), delta,
                             efficacy, NULL, NULL, equal, feasible, futility,
                             NULL, NULL, nC, nCmax, nE, opchar, Pi0, Pi1,
                             piO, ratio, summary, w, "barnard")
}

barnard_max_ess_2d_two_stage <- function(nC, nE, e1, f1, poss_x1, poss_B1) {
  max_ESS_2d <- stats::optim(par     = c(0.5, 0.5),
                             fn      = barnard_minus_ess_two_stage,
                             nC      = nC,
                             nE      = nE,
                             e1      = e1,
                             f1      = f1,
                             poss_x1 = poss_x1,
                             poss_B1 = poss_B1)
  c(max_ESS_2d$par, -max_ESS_2d$value)
}

barnard_minus_ess_two_stage  <- function(pi, nC, nE, e1, f1, poss_x1, poss_B1) {
  if (!missing(poss_x1)) {
    -barnard_des_ess_two_stage(pi, nC, nE, e1, f1, poss_x1, poss_B1)
  } else {
    -barnard_ess_two_stage(pi, nC, nE, e1, f1)
  }
}

barnard_opchar_one_stage     <- function(pi, nC, nE, e1, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi <- barnard_pmf_one_stage(pi, nC, nE, e1)
  }
  nrow_pi  <- nrow(pi)
  P        <- numeric(nrow_pi)
  for (i in 1:nrow_pi) {
    P[i]   <- sum(dplyr::filter(pmf_pi, .data$piC == pi[i, 1] &
                                  .data$piE == pi[i, 2] &
                                  .data$decision == "Reject")$`f(x,m|pi)`)
  }
  tibble::tibble(piC     = pi[, 1],
                 piE     = pi[, 2],
                 `P(pi)` = P)
}

barnard_opchar_two_stage     <- function(pi, nC, nE, e1, f1, e2, k, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi         <- barnard_pmf_two_stage(pi, nC, nE, e1, f1, e2, k)
  }
  nrow_pi          <- nrow(pi)
  n                <- c(nC[1] + nE[1], sum(nC) + sum(nE))
  opchar           <- matrix(0, nrow_pi, 13)
  E                <- Fu <- numeric(2)
  for (i in 1:nrow_pi) {
    for (j in k) {
      E[j]         <- sum(dplyr::filter(pmf_pi, .data$piC == pi[i, 1] &
                                          .data$piE == pi[i, 2] &
                                          .data$decision == "Reject" &
                                          .data$k == j)$`f(x,m|pi)`)
      Fu[j]        <- sum(dplyr::filter(pmf_pi, .data$piC == pi[i, 1] &
                                          .data$piE == pi[i, 2] &
                                          .data$decision == "Do not reject" &
                                          .data$k == j)$`f(x,m|pi)`)
    }
    cum_S          <- cumsum(S <- E + Fu)
    MSS            <- ifelse(any(cum_S == 0.5),
                             0.5*(n[which(cum_S == 0.5)] +
                                    n[which(cum_S == 0.5) + 1]),
                             n[which(cum_S > 0.5)[1]])
    ESS            <- sum(n*S)
    if (any(abs(ESS - n) < 1e-6)) {
      SDSS         <- 0
    } else {
      SDSS         <- sqrt(sum(n^2*S) - ESS^2)
    }
    opchar[i, ]    <- c(pi[i, 1], pi[i, 2], sum(E), ESS, SDSS, MSS, E, Fu, S,
                        n[2])
  }
  colnames(opchar) <- c("piC", "piE", "P(pi)", "ESS(pi)", "SDSS(pi)", "MSS(pi)",
                        paste0(rep(c("E", "F", "S"), each = 2), rep(1:2, 3),
                               "(pi)"), "max N")
  opchar           <- tibble::as_tibble(opchar)
  opchar$`max N`   <- as.integer(opchar$`max N`)
  opchar
}

barnard_pmf_one_stage        <- function(pi, nC, nE, e1) {
  x                                        <- expand.grid(0:nC, 0:nE)
  nrow_pmf                                 <- (nC + 1)*(nE + 1)
  nrow_pi                                  <- nrow(pi)
  nrow_total                               <- nrow_pmf*nrow_pi
  f                                        <- numeric(nrow_total)
  seq_0_nC                                 <- 0:nC
  if (nC == nE) {
    seq_0_nE                               <- seq_0_nC
  } else {
    seq_0_nE                               <- 0:nE
  }
  for (i in 1:nrow_pi) {
    dbinom0                                <- stats::dbinom(seq_0_nC, nC,
                                                            pi[i, 1])
    if (all(nC == nE, pi[i, 1] == pi[i, 2])) {
      dbinom1                              <- dbinom0
    } else {
      dbinom1                              <- stats::dbinom(seq_0_nE, nE,
                                                            pi[i, 2])
    }
    f[(1 + (i - 1)*nrow_pmf):(i*nrow_pmf)] <-
      dbinom0[x[, 1] + 1]*dbinom1[x[, 2] + 1]
  }
  fact                                     <- (x[, 1] + x[, 2])/(nC + nE)
  statistic                                <-
    (x[, 2]/nE - x[, 1]/nC)/sqrt(fact*(1 - fact)*(1/nC + 1/nE))
  statistic[is.nan(statistic)]             <- 0
  pmf                                      <-
    tibble::tibble(piC         = rep(pi[, 1], each = nrow_pmf),
                   piE         = rep(pi[, 2], each = nrow_pmf),
                   xC          = rep(as.vector(x[, 1]), nrow_pi),
                   xE          = rep(as.vector(x[, 2]), nrow_pi),
                   mC          = rep(as.integer(nC), nrow_total),
                   mE          = rep(as.integer(nE), nrow_total),
                   statistic   = rep(statistic, nrow_pi),
                   decision    = factor(ifelse(.data$statistic >= e1, "Reject",
                                               "Do not reject")),
                   k           = factor(rep(1, nrow_total), 1),
                   `f(x,m|pi)` = f)
  dplyr::arrange(pmf, .data$piC, .data$piE, .data$xC, .data$xE)
}

barnard_pmf_two_stage        <- function(pi, nC, nE, e1, f1, e2, k) {
  if (e1 == Inf) {
    fact                                         <- nE[1]/(nC[1] + nE[1])
    e1                                           <-
      1/sqrt(fact*(1 - fact)*(1/nC[1] + 1/nE[1])) + 1
  }
  if (f1 == -Inf) {
    fact                                         <- nC[1]/(nC[1] + nE[1])
    f1                                           <-
      -1/sqrt(fact*(1 - fact)*(1/nC[1] + 1/nE[1])) - 1
  }
  pmf                                            <-
    barnard_pmf_two_stage_cpp(pi[1, ], nC, nE, e1, f1, e2, k)
  nrow_pmf                                       <- nrow(pmf)
  nrow_pi                                        <- nrow(pi)
  if (nrow_pi > 1) {
    pmf                                          <-
      rbind(pmf, matrix(0, (nrow_pi - 1)*nrow_pmf, 8))
    for (i in 2:nrow_pi) {
      pmf[(1 + (i - 1)*nrow_pmf):(i*nrow_pmf), ] <-
        barnard_pmf_two_stage_cpp(pi[i, ], nC, nE, e1, f1, e2, k)
    }
  }
  pmf                                            <-
    tibble::tibble(piC         = rep(pi[, 1], each = nrow_pmf),
                   piE         = rep(pi[, 2], each = nrow_pmf),
                   xC          = as.integer(pmf[, 1]),
                   xE          = as.integer(pmf[, 2]),
                   mC          = as.integer(pmf[, 3]),
                   mE          = as.integer(pmf[, 4]),
                   statistic   = pmf[, 5],
                   decision    = factor(c("Do not reject",
                                          "Reject")[pmf[, 6] + 1]),
                   k           = factor(pmf[, 7], k),
                   `f(x,m|pi)` = pmf[, 8])
  dplyr::arrange(pmf, .data$piC, .data$piE, .data$k, .data$xC, .data$xE)
}

barnard_terminal_one_stage   <- function(nC, nE, e1) {
  x                            <- expand.grid(0:nC, 0:nE)
  nrow_pmf                     <- (nC + 1)*(nE + 1)
  fact                         <- (x[, 1] + x[, 2])/(nC + nE)
  statistic                    <-
    (x[, 2]/nE - x[, 1]/nC)/sqrt(fact*(1 - fact)*(1/nC + 1/nE))
  statistic[is.nan(statistic)] <- 0
  return(tibble::tibble(xC        = x[, 1],
                        xE        = x[, 2],
                        mC        = rep(as.integer(nC), nrow_pmf),
                        mE        = rep(as.integer(nE), nrow_pmf),
                        statistic = statistic,
                        decision  = factor(ifelse(statistic >= e1, "Reject",
                                                  "Do not reject")),
                        k         = factor(rep(1, nrow_pmf), 1)))
}

barnard_terminal_two_stage   <- function(nC, nE, e1, f1, e2, k) {
  if (e1 == Inf) {
    fact   <- nE[1]/(nC[1] + nE[1])
    e1     <- 1/sqrt(fact*(1 - fact)*(1/nC[1] + 1/nE[1])) + 1
  }
  if (f1 == -Inf) {
    fact   <- nC[1]/(nC[1] + nE[1])
    f1     <- -1/sqrt(fact*(1 - fact)*(1/nC[1] + 1/nE[1])) - 1
  }
  terminal <- barnard_terminal_two_stage_cpp(nC, nE, e1, f1, e2, k)
  terminal <-
    tibble::tibble(xC        = as.integer(terminal[, 1]),
                   xE        = as.integer(terminal[, 2]),
                   mC        = as.integer(terminal[, 3]),
                   mE        = as.integer(terminal[, 4]),
                   statistic = terminal[, 5],
                   decision  =
                     factor(c("Do not reject", "Reject",
                              "Continue to stage 2")[terminal[, 6]],
                            sort(unique(
                              c("Do not reject", "Reject",
                                "Continue to stage 2")[terminal[, 6]]))),
                   k         = factor(terminal[, 7], k))
  dplyr::arrange(terminal, .data$k, .data$xC, .data$xE)
}
binomial_des_one_stage        <- function(alpha, beta, delta, ratio, Pi0, Pi1,
                                          nCmax, summary) {
  params               <- search_parameters(1, "binomial", nCmax, ratio)
  feasible             <-
    binomial_des_one_stage_cpp(alpha, beta, delta, params$poss_nC,
                               params$poss_nE, params$poss_x, params$poss_y,
                               (length(Pi0) == 1), Pi0, (length(Pi1) == 1),
                               Pi1, summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message(uc("two_elip"), "feasible designs identified in the range of ",
              "considered sample sizes.\n  Identifying the optimal design",
              uc("two_elip"))
    }
    nrow_feasible      <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible         <- matrix(c(1, feasible[, 1],
                                   params$poss_nE[feasible[, 1]],
                                   feasible[, -1]), 1)
    } else {
      feasible         <- cbind(1:nrow(feasible), feasible[, 1],
                                params$poss_nE[feasible[, 1]], feasible[, -1])
    }
    colnames(feasible) <- c("index", "nC1", "nE1", "e1", "argmax alpha",
                            "max alpha", "argmin power", "min power")
    feasible           <- tibble::as_tibble(feasible)
    feasible[, 1:4]    <- dplyr::mutate_if(feasible[, 1:4], is.double,
                                           as.integer)
    feasible           <- dplyr::arrange(feasible, .data$nC1,
                                         dplyr::desc(.data$`min power`))
    e1                 <- feasible$e1[1]
    nC                 <- feasible$nC1[1]
    nE                 <- feasible$nE1[1]
    opchar             <-
      binomial_opchar_one_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                      feasible$`argmin power`[1] + c(0, delta)),
                                nC, nE, e1)
  } else {
    e1                 <- feasible <- nC <- nE <- opchar <- NULL
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size.\n  Consider increasing ",
              "nCmax", uc("two_elip"))
    }
  }
  output               <-
    build_des_one_stage_output(alpha, beta, delta, feasible, nCmax, opchar,
                               Pi0, Pi1, ratio, summary, "binomial",
                               list(e1 = e1, nC = nC, nE = nE))
  output
}

binomial_des_two_stage        <- function(alpha, beta, delta, ratio, Pi0, Pi1,
                                          nCmax, equal, w, piO, efficacy,
                                          futility, summary) {
  params               <- search_parameters(2, "binomial", nCmax, ratio)
  feasible             <-
    binomial_des_two_stage_cpp(alpha, beta, delta, params$poss_nC,
                               params$poss_nE, params$poss_x, params$poss_y,
                               (length(Pi0) == 1), Pi0, (length(Pi1) == 1), Pi1,
                               equal, efficacy, futility, piO, summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message(uc("two_elip"), "feasible designs identified in the range of ",
              "considered sample sizes.\n  Identifying the optimal design",
              uc("two_elip"))
    }
    nrow_feasible      <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible         <- matrix(c(1, feasible[, 1:2],
                                   params$poss_nE[feasible[, 1]],
                                   params$poss_nE[feasible[, 2]],
                                   feasible[, -(1:2)]), 1)
    } else {
      feasible         <- cbind(1:nrow(feasible), feasible[, 1:2],
                                params$poss_nE[feasible[, 1]],
                                params$poss_nE[feasible[, 2]],
                                feasible[, -(1:2)])
    }
    colnames(feasible) <- c("index", "nC1", "nC2", "nE1", "nE2", "e1", "e2",
                            "f1", "argmax alpha", "max alpha", "argmin power",
                            "min power", "ESS(piO,piO)", "ESS(piO,piO+delta)")
    feasible           <- tibble::as_tibble(feasible)
    feasible[, 1:8]    <- dplyr::mutate_if(feasible[, 1:8], is.double,
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
        binomial_max_ess_1d_two_stage(as.numeric(feasible[i, 2:3]),
                                      as.numeric(feasible[i, 4:5]),
                                      as.numeric(feasible[i, 6]),
                                      as.numeric(feasible[i, 8]),
                                      params$poss_x[[index]],
                                      params$poss_y[[index]])
      feasible$`argmax ESS(pi,pi)`[i]       <- max_ESS_1d[1]
      feasible$`max ESS(pi,pi)`[i]          <- max_ESS_1d[2]
      max_ESS_2d       <-
        binomial_max_ess_2d_two_stage(as.numeric(feasible[i, 2:3]),
                                      as.numeric(feasible[i, 4:5]),
                                      as.numeric(feasible[i, 6]),
                                      as.numeric(feasible[i, 8]),
                                      params$poss_x[[index]],
                                      params$poss_y[[index]])
      feasible$`argmax_piC ESS(piC,piE)`[i] <- max_ESS_2d[1]
      feasible$`argmax_piE ESS(piC,piE)`[i] <- max_ESS_2d[2]
      feasible$`max ESS(piC,piE)`[i]        <- max_ESS_2d[3]
    }
    feasible$o         <- rowSums(matrix(w, nrow_feasible, 5, T)*
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
    nC                 <- c(feasible$nC1[1], feasible$nC2[1])
    nE                 <- c(feasible$nE1[1], feasible$nE2[1])
    opchar             <-
      binomial_opchar_two_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                      feasible$`argmin power`[1] + c(0, delta)),
                                nC, nE, e1, f1, e2, 1:2)
    
  } else {
    e1                 <- e2 <- f1 <- feasible <- nC <- nE <- opchar <- NULL
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size.\n  Consider increasing ",
              "nCmax", uc("two_elip"))
    }
  }
  output               <-
    build_des_two_stage_output(alpha, beta, delta, equal, feasible, nCmax,
                               opchar, Pi0, Pi1, piO, ratio, summary, w,
                               "binomial",
                               list(e1 = e1, e2 = e2, efficacy = efficacy,
                                    f1 = f1, futility = futility, nC = nC,
                                    nE = nE))
  output
}

binomial_max_ess_2d_two_stage <- function(nC, nE, e1, f1, poss_x1, poss_y1) {
  max_ESS_2d <- stats::optim(par     = c(0.5, 0.5),
                             fn      = binomial_minus_ess_two_stage,
                             nC      = nC,
                             nE      = nE,
                             e1      = e1,
                             f1      = f1,
                             poss_x1 = poss_x1,
                             poss_y1 = poss_y1)
  c(max_ESS_2d$par, -max_ESS_2d$value)
}

binomial_minus_ess_two_stage  <- function(pi, nC, nE, e1, f1, poss_x1,
                                          poss_y1) {
  if (!missing(poss_x1)) {
    -binomial_des_ess_two_stage(pi, nC, nE, e1, f1, poss_x1, poss_y1)
  } else {
    -binomial_ess_two_stage(pi, nC, nE, e1, f1)
  }
}

binomial_opchar_one_stage     <- function(pi, nC, nE, e, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi <- binomial_pmf_one_stage(pi, nC, nE, e)
  }
  rows_pi  <- nrow(pi)
  P        <- numeric(rows_pi)
  for (i in 1:rows_pi) {
    P[i]   <- sum(dplyr::filter(pmf_pi, .data$piC == pi[i, 1] &
                                  .data$piE == pi[i, 2] &
                                  .data$decision == "Reject")$`f(x,m|pi)`)
  }
  tibble::tibble(piC     = pi[, 1],
                 piE     = pi[, 2],
                 `P(pi)` = P)
}

binomial_opchar_two_stage     <- function(pi, nC, nE, e1, f1, e2, k, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi         <- binomial_pmf_two_stage(pi, nC, nE, e1, f1, e2, k)
  }
  rows_pi          <- nrow(pi)
  n                <- c(nC[1] + nE[1], sum(nC) + sum(nE))
  opchar           <- matrix(0, rows_pi, 13)
  E                <- Fu <- numeric(2)
  for (i in 1:rows_pi) {
    pmf_i          <- dplyr::filter(pmf_pi, .data$piC == pi[i, 1] &
                                      .data$piE == pi[i, 2])
    for (j in k) {
      E[j]         <- sum(dplyr::filter(pmf_i, .data$decision == "Reject" &
                                          .data$k == j)$`f(x,m|pi)`)
      Fu[j]        <- sum(dplyr::filter(pmf_i, .data$decision == "Do not reject" &
                                          .data$k == j)$`f(x,m|pi)`)
    }
    cum_S          <- cumsum(S <- E + Fu)
    MSS            <- ifelse(any(cum_S == 0.5),
                             0.5*(n[which(cum_S == 0.5)] +
                                    n[which(cum_S == 0.5) + 1]),
                             n[which(cum_S > 0.5)[1]])
    ESS            <- sum(n*S)
    if (any(abs(ESS - n) < 1e-13)) {
      SDSS         <- 0
    } else {
      SDSS         <- sqrt(sum(n^2*S) - sum(n*S)^2)
    }
    opchar[i, ]    <- c(pi[i, 1], pi[i, 2], sum(E), ESS, SDSS, MSS, E, Fu, S,
                        n[2])
  }
  opchar           <- tibble::as_tibble(opchar)
  colnames(opchar) <- c("piC", "piE", "P(pi)", "ESS(pi)", "SDSS(pi)", "MSS(pi)",
                        paste0(rep(c("E", "F", "S"), each = 2), rep(1:2, 3),
                               "(pi)"), "max N")
  opchar$`max N`   <- as.integer(opchar$`max N`)
  opchar
}

binomial_pmf_one_stage        <- function(pi, nC, nE, e1) {
  x                                        <- expand.grid(0:nC, 0:nE)
  rows_pmf                                 <- (nC + 1)*(nE + 1)
  rows_pi                                  <- nrow(pi)
  rows_total                               <- rows_pmf*rows_pi
  f                                        <- numeric(rows_total)
  for (i in 1:rows_pi) {
    dbinom0                                <- stats::dbinom(0:nC, nC, pi[i, 1])
    if (all(nC == nE, pi[i, 1] == pi[i, 2])) {
      dbinom1                              <- dbinom0
    } else {
      dbinom1                              <- stats::dbinom(0:nE, nE, pi[i, 2])
    }
    f[(1 + (i - 1)*rows_pmf):(i*rows_pmf)] <-
      dbinom0[x[, 1] + 1]*dbinom1[x[, 2] + 1]
  }
  pmf                                      <-
    tibble::tibble(piC         = rep(pi[, 1], each = rows_pmf),
                   piE         = rep(pi[, 2], each = rows_pmf),
                   xC          = rep(as.vector(x[, 1]), rows_pi),
                   xE          = rep(as.vector(x[, 2]), rows_pi),
                   mC          = rep(as.integer(nC), rows_total),
                   mE          = rep(as.integer(nE), rows_total),
                   statistic   = .data$xE - .data$xC,
                   decision    = factor(ifelse(.data$statistic >= e1, "Reject",
                                               "Do not reject")),
                   k           = factor(rep(1, rows_total), 1),
                   `f(x,m|pi)` = f)
  dplyr::arrange(pmf, .data$piC, .data$piE, .data$xC, .data$xE)
}

binomial_pmf_two_stage        <- function(pi, nC, nE, e1, f1, e2, k) {
  if (e1 == Inf) {
    e1                                           <- nE[1] + 1
  }
  if (f1 == -Inf) {
    f1                                           <- -nC[1] - 1
  }
  pmf                                            <-
    binomial_pmf_two_stage_cpp(pi[1, ], nC, nE, e1, f1, e2, k)
  rows_pmf                                       <- nrow(pmf)
  rows_pi                                        <- nrow(pi)
  if (rows_pi > 1) {
    pmf                                          <-
      rbind(pmf, matrix(0, (rows_pi - 1)*rows_pmf, 8))
    for (i in 2:rows_pi) {
      pmf[(1 + (i - 1)*rows_pmf):(i*rows_pmf), ] <-
        binomial_pmf_two_stage_cpp(pi[i, ], nC, nE, e1, f1, e2, k)
    }
  }
  pmf                                            <-
    tibble::tibble(piC         = rep(pi[, 1], each = rows_pmf),
                   piE         = rep(pi[, 2], each = rows_pmf),
                   xC          = as.integer(pmf[, 1]),
                   xE          = as.integer(pmf[, 2]),
                   mC          = as.integer(pmf[, 3]),
                   mE          = as.integer(pmf[, 4]),
                   statistic   = as.integer(pmf[, 5]),
                   decision    = factor(c("Do not reject",
                                          "Reject")[pmf[, 6] + 1]),
                   k           = factor(pmf[, 7], k),
                   `f(x,m|pi)` = pmf[, 8])
  dplyr::arrange(pmf, .data$piC, .data$piE, .data$k, .data$xC, .data$xE)
}

binomial_terminal_one_stage   <- function(nC, nE, e1) {
  x        <- expand.grid(0:nC, 0:nE)
  rows_pmf <- (nC + 1)*(nE + 1)
  tibble::tibble(xC        = x[, 1],
                 xE        = x[, 2],
                 mC        = rep(as.integer(nC), rows_pmf),
                 mE        = rep(as.integer(nE), rows_pmf),
                 statistic = .data$xE - .data$xC,
                 decision  = factor(ifelse(.data$statistic >= e1, "Reject",
                                           "Do not reject")),
                 k         = factor(rep(1, rows_pmf), 1))
}

binomial_terminal_two_stage   <- function(nC, nE, e1, f1, e2, k) {
  if (e1 == Inf) {
    e1     <- nE[1] + 1
  }
  if (f1 == -Inf) {
    f1     <- -nC[1] - 1
  }
  terminal <- binomial_terminal_two_stage_cpp(nC, nE, e1, f1, e2, k)
  terminal <- tibble::tibble(xC        = as.integer(terminal[, 1]),
                             xE        = as.integer(terminal[, 2]),
                             mC        = as.integer(terminal[, 3]),
                             mE        = as.integer(terminal[, 4]),
                             statistic = as.integer(terminal[, 5]),
                             decision  =
                               factor(c("Do not reject", "Reject",
                                        "Continue to stage 2")[terminal[, 6]],
                                      sort(unique(c("Do not reject", "Reject",
                                                    "Continue to stage 2")[terminal[, 6]]))),
                             k         = factor(terminal[, 7], k))
  dplyr::arrange(terminal, .data$k, .data$xC, .data$xE)
}
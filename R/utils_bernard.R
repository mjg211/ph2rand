bernard_des_one_stage        <- function(alpha, beta, delta, ratio, point_null,
                                         pi_null, point_alt, pi_alt, n0max,
                                         summary) {
  params               <- search_parameters(1, "bernard", n0max, ratio)
  feasible             <-
    bernard_des_one_stage_cpp(alpha, beta, delta, params$poss_n0,
                              params$poss_n1, params$poss_x, params$poss_B,
                              params$unique_B, point_null, pi_null, point_alt,
                              pi_alt, summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message(uc("two_elip"), "feasible designs identified in the range of ",
              "considered sample sizes. Identifying the optimal design",
              uc("two_elip"))
    }
    nrow_feasible      <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible         <- matrix(c(1, feasible[, 1],
                                   params$poss_n1[feasible[, 1]],
                                   feasible[, -1]), 1)
    } else {
      feasible         <- cbind(1:nrow(feasible), feasible[, 1],
                                params$poss_n1[feasible[, 1]], feasible[, -1])
    }
    colnames(feasible) <- c("index", "n0", "n1", "e", "argmax alpha",
                            "max alpha", "argmin power", "min power")
    feasible           <- tibble::as_tibble(feasible)
    feasible[, 1:3]    <- dplyr::mutate_if(feasible[, 1:3], is.double,
                                           as.integer)
    feasible           <- dplyr::arrange(feasible, n0, desc(`min power`))
    e                  <- feasible$e[1]
    n0                 <- feasible$n0[1]
    n1                 <- feasible$n1[1]
    opchar             <-
      bernard_opchar_one_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                     feasible$`argmin power`[1] + c(0, delta)),
                               n0, n1, e)
  } else {
    e                  <- feasible <- n0 <- n1 <- opchar <- NULL
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size. Consider increasing ",
              "n0max", uc("two_elip"))
    }
  }
  output               <-
    build_des_one_stage_output(alpha, beta, delta, feasible, n0max, opchar,
                               pi_alt, pi_null, point_alt, point_null, ratio,
                               summary, "bernard", list(e = e, n0 = n0,
                                                        n1 = n1))
  output
}

bernard_des_two_stage        <- function(alpha, beta, delta, ratio, point_null,
                                         pi_null, point_alt, pi_alt, n0max,
                                         equal, w, pi_ess, efficacy, futility,
                                         summary) {
  params               <- search_parameters(2, "bernard", n0max, ratio)
  feasible             <-
    bernard_des_two_stage_cpp(alpha, beta, delta, params$poss_n0,
                              params$poss_n1, params$poss_x, params$poss_B,
                              params$unique_B, point_null, pi_null, point_alt,
                              pi_alt, equal, efficacy, futility, pi_ess,
                              summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message(uc("two_elip"), "feasible designs identified in the range of ",
              "considered sample sizes. Identifying the optimal design",
              uc("two_elip"))
    }
    nrow_feasible      <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible         <- matrix(c(1, feasible[, 1:2],
                                   params$poss_n1[feasible[, 1]],
                                   params$poss_n1[feasible[, 2]],
                                   feasible[, -(1:2)]), 1)
    } else {
      feasible         <- cbind(1:nrow(feasible), feasible[, 1:2],
                                params$poss_n1[feasible[, 1]],
                                params$poss_n1[feasible[, 2]],
                                feasible[, -(1:2)])
    }
    colnames(feasible) <- c("index", "n01", "n02", "n11", "n12", "e1", "e2",
                            "f1", "argmax alpha", "max alpha", "argmin power",
                            "min power", "ESS0", "ESS1")
    feasible           <- tibble::as_tibble(feasible)
    feasible[, 1:5]    <- dplyr::mutate_if(feasible[, 1:5], is.double,
                                           as.integer)
    feasible           <-
      dplyr::mutate(feasible,
                    `argmax ESS(pi,pi)`       = 0,
                    `max ESS(pi,pi)`          = 0,
                    `argmax_pi0 ESS(pi0,pi1)` = 0,
                    `argmax_pi1 ESS(pi0,pi1)` = 0,
                    `max ESS(pi0,pi1)`        = 0,
                    `max(n)`                  =
                      as.integer(rowSums(feasible[, 2:5])))
    for (i in 1:nrow_feasible) {
      index            <- as.numeric(feasible[i, 2] +
                                       params$max_poss_n0*(feasible[i, 4] - 1))
      max_ESS_1d       <-
        bernard_max_ess_1d_two_stage(as.numeric(feasible[i, 2:3]),
                                     as.numeric(feasible[i, 4:5]),
                                     as.numeric(feasible[i, 6]),
                                     as.numeric(feasible[i, 8]),
                                     params$poss_x[[index]],
                                     params$poss_B[[index]])
      feasible$`argmax ESS(pi,pi)`[i]       <- max_ESS_1d[1]
      feasible$`max ESS(pi,pi)`[i]          <- max_ESS_1d[2]
      max_ESS_2d       <-
        bernard_max_ess_2d_two_stage(as.numeric(feasible[i, 2:3]),
                                     as.numeric(feasible[i, 4:5]),
                                     as.numeric(feasible[i, 6]),
                                     as.numeric(feasible[i, 8]),
                                     params$poss_x[[index]],
                                     params$poss_B[[index]])
      feasible$`argmax_pi0 ESS(pi0,pi1)`[i] <- max_ESS_2d[1]
      feasible$`argmax_pi1 ESS(pi0,pi1)`[i] <- max_ESS_2d[2]
      feasible$`max ESS(pi0,pi1)`[i] <- max_ESS_2d[3]
    }
    feasible$o         <- rowSums(matrix(w, nrow_feasible, 5, byrow = T)*
                                    feasible[, c(13, 14, 16, 19, 20)])
    feasible           <- dplyr::arrange(feasible, o, desc(`min power`))
    if (!efficacy) {
      feasible$e1      <- Inf
    }
    if (!futility) {
      feasible$f1      <- -Inf
    }
    e1                 <- feasible$e1[1]
    e2                 <- feasible$e2[1]
    f1                 <- feasible$f1[1]
    n0                 <- c(feasible$n01[1], feasible$n02[1])
    n1                 <- c(feasible$n11[1], feasible$n12[1])
    opchar             <-
      bernard_opchar_two_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                     feasible$`argmin power`[1] + c(0, delta)),
                               n0, n1, e1, f1, e2, 1:2)
  } else {
    e1                 <- e2 <- f1 <- feasible <- n0 <- n1 <- opchar <- NULL
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size. Consider increasing ",
              "n0max", uc("two_elip"))
    }
  }
  output               <-
    build_des_two_stage_output(alpha, beta, delta, equal, feasible, n0max,
                               opchar, pi_alt, pi_ess, pi_null, point_alt,
                               point_null, ratio, summary, w, "bernard",
                               list(e1 = e1, e2 = e2, efficacy = efficacy,
                                    f1 = f1, futility = futility, n0 = n0,
                                    n1 = n1))
  output
}

bernard_opchar_one_stage     <- function(pi, n0, n1, e, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi <- bernard_pmf_one_stage(pi, n0, n1, e)
  }
  rows_pi  <- nrow(pi)
  P        <- numeric(rows_pi)
  for (i in 1:rows_pi) {
    P[i]   <- sum(dplyr::filter(pmf_pi, pi0 == pi[i, 1] & pi1 == pi[i, 2] &
                                  decision == "Reject")$`f(x,m|pi)`)
  }
  tibble::tibble(pi0           = pi[, 1],
                 pi1           = pi[, 2],
                 `P(pi)`       = P)
}

bernard_opchar_two_stage     <- function(pi, n0, n1, e1, f1, e2, k, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi         <- bernard_pmf_two_stage(pi, n0, n1, e1, f1, e2, k)
  }
  rows_pi          <- nrow(pi)
  n                <- c(n0[1] + n1[1], sum(n0) + sum(n1))
  opchar           <- matrix(0, rows_pi, 15)
  E                <- Fu <- numeric(2)
  for (i in 1:rows_pi) {
    for (j in k) {
      E[j]         <- sum(dplyr::filter(pmf_pi, pi0 == pi[i, 1] &
                                          pi1 == pi[i, 2] &
                                          decision == "Reject" &
                                          k == j)$`f(x,m|pi)`)
      Fu[j]        <- sum(dplyr::filter(pmf_pi, pi0 == pi[i, 1] &
                                          pi1 == pi[i, 2] &
                                          decision == "Do not reject" &
                                          k == j)$`f(x,m|pi)`)
    }
    cum_S          <- cumsum(S <- E + Fu)
    MSS            <- ifelse(any(cum_S == 0.5),
                             0.5*(n[which(cum_S == 0.5)] +
                                    n[which(cum_S == 0.5) + 1]),
                             n[which(cum_S > 0.5)[1]])
    opchar[i, ]    <- c(pi[i, 1], pi[i, 2], sum(E), sum(n*S),
                        sqrt(sum(n^2*S) - sum(n*S)^2), MSS, E, Fu, S, cum_S,
                        n[2])
  }
  opchar           <- tibble::as_tibble(opchar)
  colnames(opchar) <- c("pi0", "pi1", "P(pi)", "ESS(pi)", "SDSS(pi)", "MSS(pi)",
                        paste(rep(c("E", "F", "S"), each = 2), rep(1:2, 3),
                              "(pi)", sep = ""),
                        paste("cum{S", 1:2, "(pi)}", sep = ""), "max(n)")
  opchar$`max(n)`  <- as.integer(opchar$`max(n)`)
  opchar
}

bernard_max_ess_2d_two_stage <- function(n0, n1, e1, f1, poss_x1, poss_B1) {
  max_ESS_2d <- stats::optim(par     = c(0.5, 0.5),
                             fn      = bernard_minus_ess_two_stage,
                             n0      = n0,
                             n1      = n1,
                             e1      = e1,
                             f1      = f1,
                             poss_x1 = poss_x1,
                             poss_B1 = poss_B1)
  c(max_ESS_2d$par, -max_ESS_2d$value)
}

bernard_minus_ess_two_stage  <- function(pi, n0, n1, e1, f1, poss_x1, poss_B1) {
  if (!missing(poss_x1)) {
    -bernard_des_ess_two_stage(pi, n0, n1, e1, f1, poss_x1, poss_B1)
  } else {
    -bernard_ess_two_stage(pi, n0, n1, e1, f1)
  }
}

bernard_pmf_one_stage        <- function(pi, n0, n1, e) {
  x                                        <- expand.grid(0:n0, 0:n1)
  rows_pmf                                 <- (n0 + 1)*(n1 + 1)
  rows_pi                                  <- nrow(pi)
  rows_total                               <- rows_pmf*rows_pi
  f                                        <- numeric(rows_total)
  for (i in 1:rows_pi) {
    dbinom0                                <- dbinom(0:n0, n0, pi[i, 1])
    if (all(n0 == n1, pi[i, 1] == pi[i, 2])) {
      dbinom1                              <- dbinom0
    } else {
      dbinom1                              <- dbinom(0:n1, n1, pi[i, 2])
    }
    f[(1 + (i - 1)*rows_pmf):(i*rows_pmf)] <-
      dbinom0[x[, 1] + 1]*dbinom1[x[, 2] + 1]
  }
  fact                                     <- (x[, 1] + x[, 2])/(n0 + n1)
  statistic                                <-
    (x[, 2]/n1 - x[, 1]/n0)/sqrt(fact*(1 - fact)*(1/n0 + 1/n1))
  statistic[is.nan(statistic)]             <- 0
  pmf                                      <-
    tibble::tibble(pi0         = rep(pi[, 1], each = rows_pmf),
                   pi1         = rep(pi[, 2], each = rows_pmf),
                   x0          = rep(as.vector(x[, 1]), rows_pi),
                   x1          = rep(as.vector(x[, 2]), rows_pi),
                   m0          = rep(as.integer(n0), rows_total),
                   m1          = rep(as.integer(n1), rows_total),
                   statistic   = rep(statistic, rows_pi),
                   decision    = ifelse(statistic >= e, "Reject",
                                        "Do not reject"),
                   k           = factor(rep(1, rows_total), 1),
                   `f(x,m|pi)` = f)
  dplyr::arrange(pmf, pi0, pi1, x0, x1)
}

bernard_pmf_two_stage        <- function(pi, n0, n1, e1, f1, e2, k) {
  if (e1 == Inf) {
    fact                                         <- n1[1]/(n0[1] + n1[1])
    e1                                           <-
      1/sqrt(fact*(1 - fact)*(1/n0[1] + 1/n1[1])) + 1
  }
  if (f1 == -Inf) {
    fact                                         <- n0[1]/(n0[1] + n1[1])
    f1                                           <-
      -1/sqrt(fact*(1 - fact)*(1/n0[1] + 1/n1[1])) - 1
  }
  pmf                                            <-
    bernard_pmf_two_stage_cpp(pi[1, ], n0, n1, e1, f1, e2, k)
  rows_pmf                                       <- nrow(pmf)
  rows_pi                                        <- nrow(pi)
  if (rows_pi > 1) {
    pmf                                          <-
      rbind(pmf, matrix(0, (rows_pi - 1)*rows_pmf, 8))
    for (i in 2:rows_pi) {
      pmf[(1 + (i - 1)*rows_pmf):(i*rows_pmf), ] <-
        bernard_pmf_two_stage_cpp(pi[i, ], n0, n1, e1, f1, e2, k)
    }
  }
  pmf                                            <-
    tibble::tibble(pi0         = rep(pi[, 1], each = rows_pmf),
                   pi1         = rep(pi[, 2], each = rows_pmf),
                   x0          = as.integer(pmf[, 1]),
                   x1          = as.integer(pmf[, 2]),
                   m0          = as.integer(pmf[, 3]),
                   m1          = as.integer(pmf[, 4]),
                   statistic   = pmf[, 5],
                   decision    = ifelse(pmf[, 6] == 1, "Reject",
                                        "Do not reject"),
                   k           = factor(pmf[, 7], k),
                   `f(x,m|pi)` = pmf[, 8])
  dplyr::arrange(pmf, pi0, pi1, k, x0, x1)
}

bernard_terminal_one_stage   <- function(n0, n1, e) {
  x                            <- expand.grid(0:n0, 0:n1)
  rows_pmf                     <- (n0 + 1)*(n1 + 1)
  fact                         <- (x[, 1] + x[, 2])/(n0 + n1)
  statistic                    <-
    (x[, 2]/n1 - x[, 1]/n0)/sqrt(fact*(1 - fact)*(1/n0 + 1/n1))
  statistic[is.nan(statistic)] <- 0
  return(tibble::tibble(x0        = x[, 1],
                        x1        = x[, 2],
                        m0        = rep(as.integer(n0), rows_pmf),
                        m1        = rep(as.integer(n1), rows_pmf),
                        statistic = statistic,
                        decision  = ifelse(statistic >= e, "Reject",
                                           "Do not reject"),
                        k         = factor(rep(1, rows_pmf), 1)))
}

bernard_terminal_two_stage   <- function(n0, n1, e1, f1, e2, k) {
  if (e1 == Inf) {
    fact   <- n1[1]/(n0[1] + n1[1])
    e1     <- 1/sqrt(fact*(1 - fact)*(1/n0[1] + 1/n1[1])) + 1
  }
  if (f1 == -Inf) {
    fact   <- n0[1]/(n0[1] + n1[1])
    f1     <- -1/sqrt(fact*(1 - fact)*(1/n0[1] + 1/n1[1])) - 1
  }
  terminal <- bernard_terminal_two_stage_cpp(n0, n1, e, f, k)
  terminal <- tibble::tibble(x0        = as.integer(terminal[, 1]),
                             x1        = as.integer(terminal[, 2]),
                             m0        = as.integer(terminal[, 3]),
                             m1        = as.integer(terminal[, 4]),
                             statistic = terminal[, 5],
                             decision  = ifelse(terminal[, 6] == 1, "Reject",
                                                "Do not reject"),
                             k  = factor(terminal[, 7], k))
  dplyr::arrange(terminal, k, x0, x1)
}
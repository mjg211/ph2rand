fisher_des_one_stage        <- function(alpha, beta, delta, ratio, point_null,
                                        pi_null, point_alt, pi_alt, n0max,
                                        summary) {
  params                 <- search_parameters(1, "fisher", n0max, ratio)
  feasible               <-
    fisher_des_one_stage_cpp(alpha, beta, delta, params$poss_n0,
                             params$poss_n1, params$poss_x, params$poss_y,
                             params$poss_z, params$choose_mat, point_null,
                             pi_null, point_alt, pi_alt, summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message(uc("two_elip"), "feasible designs identified in the range of ",
              "considered sample sizes. Identifying the optimal design",
              uc("two_elip"))
    }
    ncol_feasible        <- ncol(feasible)
    nrow_feasible        <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible_e         <- matrix(c(1, feasible[, 2:(ncol_feasible - 4)]), 1)
      feasible           <- matrix(c(1, feasible[, 1],
                                     params$poss_n1[feasible[, 1]],
                                     feasible[, (ncol_feasible - 3):
                                                ncol_feasible]), 1)
    } else {
      feasible_e         <- cbind(1:nrow_feasible,
                                  feasible[, 2:(ncol_feasible - 4)])
      feasible           <- cbind(1:nrow(feasible), feasible[, 1],
                                  params$poss_n1[feasible[, 1]],
                                  feasible[, (ncol_feasible - 3):ncol_feasible])
    }
    colnames(feasible)   <- c("index", "n0", "n1", "argmax alpha", "max alpha",
                              "argmin power", "min power")
    feasible             <- tibble::as_tibble(feasible)
    feasible[, 1:3]      <- dplyr::mutate_if(feasible[, 1:3], is.double,
                                             as.integer)
    feasible_e           <- feasible_e[, 1:(max(rowSums(feasible[, 2:3])) + 2)]
    colnames(feasible_e) <- c("index",
                              paste0("e", 0:max(rowSums(feasible[, 2:3]))))
    feasible_e           <- tibble::as_tibble(feasible_e)
    feasible_e           <- dplyr::mutate_if(feasible_e, is.double, as.integer)
    e                    <- as.integer(feasible_e[1, 2:(feasible$n0[1] +
                                                          feasible$n1[1] + 2)])
    n0                   <- feasible$n0[1]
    n1                   <- feasible$n1[1]
    opchar               <-
      fisher_opchar_one_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                    feasible$`argmin power`[1] + c(0, delta)),
                              n0, n1, e)
  } else {
    e                    <- feasible <- feasible_e <- n0 <- n1 <- opchar <- NULL
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size. Consider increasing ",
              "n0max", uc("two_elip"))
    }
  }
  output                 <-
    build_des_one_stage_output(alpha, beta, delta, feasible, n0max, opchar,
                               pi_alt, pi_null, point_alt, point_null, ratio,
                               summary, "fisher",
                               list(e = e, feasible_e = feasible_e, n0 = n0,
                                    n1 = n1))
  output
}

fisher_des_two_stage        <- function(alpha, beta, delta, ratio, point_null,
                                        pi_null, point_alt, pi_alt, n0max,
                                        equal, w, pi_ess, efficacy_type,
                                        efficacy_param, futility_type,
                                        futility_param, summary) {
  params                         <- search_parameters(2, "fisher", n0max, ratio)
  feasible                       <-
    fisher_des_two_stage_cpp(alpha, beta, delta, params$poss_n0, params$poss_n1,
                             params$poss_x, params$poss_y, params$poss_z,
                             params$choose_mat, point_null, pi_null, point_alt,
                             pi_alt, equal, efficacy_type, efficacy_param,
                             futility_type, futility_param, pi_ess, summary)
  if (feasible[[1]][1, 1] > 0) {
    if (summary) {
      message(uc("two_elip"), "feasible designs identified in the range of ",
              "considered sample sizes. Identifying the optimal design",
              uc("two_elip"))
    }
    nrow_feasible                <- nrow(feasible[[1]])
    feasible_e2_vec              <- feasible[[4]]
    if (nrow_feasible == 1) {
      feasible_e1                <- matrix(c(1, feasible[[2]]), 1)
      feasible_f1                <- matrix(c(1, feasible[[3]]), 1)
      feasible                   <- feasible[[1]]
      feasible                   <- matrix(c(1:nrow_feasible, feasible[, 1:2],
                                             params$poss_n1[feasible[, 1]],
                                             params$poss_n1[feasible[, 2]],
                                             feasible[, -(1:2)]), 1)
    } else {
      feasible_e1                <- cbind(1:nrow_feasible, feasible[[2]])
      feasible_f1                <- cbind(1:nrow_feasible, feasible[[3]])
      feasible                   <- feasible[[1]]
      feasible                   <- cbind(1:nrow_feasible, feasible[, 1:2],
                                          params$poss_n1[feasible[, 1]],
                                          params$poss_n1[feasible[, 2]],
                                          feasible[, -(1:2)])
    }
    colnames(feasible)           <- c("index", "n01", "n02", "n11", "n12",
                                      "argmax alpha", "max alpha",
                                      "argmin power", "min power", "ESS0",
                                      "ESS1")
    feasible                     <- tibble::as_tibble(feasible)
    feasible[, 1:5]              <- dplyr::mutate_if(feasible[, 1:5], is.double,
                                                     as.integer)
    feasible_e1                  <-
      feasible_e1[, 1:(max(rowSums(feasible[, c(2, 4)])) + 2)]
    colnames(feasible_e1)        <-
      c("index", paste0("e", 0:max(rowSums(feasible[, c(2, 4)]))))
    feasible_e1                  <- tibble::as_tibble(feasible_e1)
    feasible_f1                  <-
      feasible_f1[, 1:(max(rowSums(feasible[, c(2, 4)])) + 2)]
    colnames(feasible_f1)        <-
      c("index", paste0("f", 0:max(rowSums(feasible[, c(2, 4)]))))
    feasible_f1                  <- tibble::as_tibble(feasible_f1)
    feasible_e2                  <- list()
    for (i in 1:nrow_feasible) {
      n0                         <- as.numeric(feasible[i, 2:3])
      n1                         <- as.numeric(feasible[i, 4:5])
      e2_i                       <- matrix(0L, n0[1] + n1[1] + 1,
                                           n0[2] + n1[2] + 1)
      counter                    <- 1
      for (z1p in 1:(max(params$poss_n0) + max(params$poss_n1) + 1)) {
        for (z2p in 1:(2*(max(params$poss_n0) + max(params$poss_n1)) + 1)) {
          if (all(z1p <= n0[1] + n1[1] + 1, z2p <= n0[2] + n1[2] + 1)) {
            if (feasible_e2_vec[i, counter] != -0.5) {
              e2_i[z1p, z2p]     <- feasible_e2_vec[i, counter]
            } else if (feasible_e2_vec[i, counter] == z1p + z2p - 1) {
              e2_i[z1p, z2p]     <- Inf
            } else {
              e2_i[z1p, z2p]     <- NA_integer_
            }
          }
          counter                <- counter + 1
        }
      }
      feasible_e2[[i]]           <- e2_i
    }
    feasible                     <-
      dplyr::mutate(feasible,
                    `argmax ESS(pi,pi)`       = 0,
                    `max ESS(pi,pi)`          = 0,
                    `argmax_pi0 ESS(pi0,pi1)` = 0,
                    `argmax_pi1 ESS(pi0,pi1)` = 0,
                    `max ESS(pi0,pi1)`        = 0,
                    `max(n)`                  =
                      as.integer(rowSums(feasible[, 2:5])))
    for (i in 1:nrow_feasible) {
      index                      <-
        as.numeric(feasible[i, 2] + params$max_poss_n0*(feasible[i, 4] - 1))
      max_ESS_1d                 <-
        fisher_max_ess_1d_two_stage(as.numeric(feasible[i, 2:3]),
                                    as.numeric(feasible[i, 4:5]),
                                    as.numeric(feasible_e1[i, ]),
                                    as.numeric(feasible_f1[i, ]),
                                    params$poss_x[[index]],
                                    params$poss_y[[index]],
                                    params$poss_z[[index]])
      feasible$`argmax ESS(pi,pi)`[i]       <- max_ESS_1d[1]
      feasible$`max ESS(pi,pi)`[i]          <- max_ESS_1d[2]
      max_ESS_2d                 <-
        fisher_max_ess_2d_two_stage(as.numeric(feasible[i, 2:3]),
                                    as.numeric(feasible[i, 4:5]),
                                    as.numeric(feasible_e1[i, ]),
                                    as.numeric(feasible_f1[i, ]),
                                    params$poss_x[[index]],
                                    params$poss_y[[index]],
                                    params$poss_z[[index]])
      feasible$`argmax_pi0 ESS(pi0,pi1)`[i] <- max_ESS_2d[1]
      feasible$`argmax_pi1 ESS(pi0,pi1)`[i] <- max_ESS_2d[2]
      feasible$`max ESS(pi0,pi1)`[i]        <- max_ESS_2d[3]
    }
    feasible$o                   <-
      rowSums(matrix(w, nrow_feasible, 5, byrow = T)*
                feasible[, c(10, 11, 13, 16, 17)])
    feasible                     <-
      dplyr::arrange(feasible, .data$o, dplyr::desc(.data$`min power`))
    feasible_e1                  <- feasible_e1[as.matrix(feasible[, 1]), ]
    feasible_f1                  <- feasible_f1[as.matrix(feasible[, 1]), ]
    ncol_feasible_e1             <- ncol(feasible_e1)
    for (i in 1:nrow_feasible) {
      n01                        <- feasible$n01[i]
      n11                        <- feasible$n11[i]
      for (z1 in 0:(n01 + n11)) {
        if (feasible_e1[i, z1 + 2] >= z1 + 1) {
          feasible_e1[i, z1 + 2] <- Inf
        }
        if (feasible_f1[i, z1 + 2] <= -z1 - 1) {
          feasible_f1[i, z1 + 2] <- -Inf
        }
      }
      if (n01 + n11 + 2 < ncol_feasible_e1) {
        range                    <- (n01 + n11 + 3):ncol_feasible_e1
        feasible_e1[i, range]    <- feasible_f1[i, range] <- NA_integer_
      }
    }
    e1                           <-
      as.numeric(feasible_e1[feasible$index[1],
                             2:(feasible$n01[1] + feasible$n11[1] + 2)])
    e2                           <- feasible_e2[[feasible$index[1]]]
    f1                           <-
      as.numeric(feasible_f1[feasible$index[1],
                             2:(feasible$n01[1] + feasible$n11[1] + 2)])
    n0                           <- c(feasible$n01[1],feasible$n02[1])
    n1                           <- c(feasible$n11[1],feasible$n12[1])
    opchar                       <-
      fisher_opchar_two_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                    feasible$`argmin power`[1] + c(0, delta)),
                              n0, n1, e1, f1, e2, 1:2)
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size. Consider increasing ",
              "n0max", uc("two_elip"))
    }
  } else {
    
    
    if (summary) {
      message(uc("two_elip"), "no feasible designs identified in range of ",
              "considered maximal allowed sample size. Consider increasing ",
              "n0max", uc("two_elip"))
    }
  }
  output               <-
    build_des_two_stage_output(alpha, beta, delta, equal, feasible, n0max,
                               opchar, pi_alt, pi_ess, pi_null, point_alt,
                               point_null, ratio, summary, w, "fisher",
                               list(e1 = e1, e2 = e2,
                                    efficacy_param = efficacy_param,
                                    efficacy_type = efficacy_type,
                                    f1 = f1, feasible_e1 = feasible_e1,
                                    feasible_e2 = feasible_e2,
                                    feasible_f1 = feasible_f1, n0 = n0,
                                    n1 = n1))
  output
}

fisher_max_ess_2d_two_stage <- function(n0, n1, e1, f1, poss_x1, poss_y1,
                                        poss_z1) {
  max_ESS_2d <- stats::optim(par     = c(0.5, 0.5),
                             fn      = fisher_minus_ess_two_stage,
                             n0      = n0,
                             n1      = n1,
                             e1      = e1,
                             f1      = f1,
                             poss_x1 = poss_x1,
                             poss_y1 = poss_y1,
                             poss_z1 = poss_z1)
  c(max_ESS_2d$par, -max_ESS_2d$value)
}

fisher_minus_ess_two_stage  <- function(pi, n0, n1, e1, f1, poss_x1, poss_y1,
                                        poss_z1) {
  if (!missing(poss_x1)) {
    -fisher_des_ess_two_stage(pi, n0, n1, e1, f1, poss_x1, poss_y1, poss_z1)
  } else {
    -fisher_ess_two_stage(pi, n0, n1, e1, f1)
  }
}

fisher_opchar_one_stage     <- function(pi, n0, n1, e, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi <- fisher_pmf_one_stage(pi, n0, n1, e)
  }
  rows_pi  <- nrow(pi)
  P        <- numeric(rows_pi)
  for (i in 1:rows_pi) {
    P[i]   <- sum(dplyr::filter(pmf_pi, .data$pi0 == pi[i, 1] &
                                  .data$pi1 == pi[i, 2] &
                                  .data$decision == "Reject")$`f(x,m|pi)`)
  }
  tibble::tibble(pi0           = pi[, 1],
                 pi1           = pi[, 2],
                 `P(pi)`       = P)
}

fisher_opchar_two_stage     <- function(pi, n0, n1, e1, f1, e2, k, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi         <- fisher_pmf_two_stage(pi, n0, n1, e1, f1, e2, k)
  }
  rows_pi          <- nrow(pi)
  n                <- c(n0[1] + n1[1], sum(n0) + sum(n1))
  opchar           <- matrix(0, nrow = rows_pi, ncol = 15)
  E                <- Fu <- numeric(2)
  for (i in 1:rows_pi) {
    for (j in k) {
      E[j]         <- sum(dplyr::filter(pmf_pi, .data$pi0 == pi[i, 1] &
                                          .data$pi1 == pi[i, 2] &
                                          .data$decision == "Reject" &
                                          .data$k == j)$`f(x,m|pi)`)
      Fu[j]        <- sum(dplyr::filter(pmf_pi, .data$pi0 == pi[i, 1] &
                                          .data$pi1 == pi[i, 2] &
                                          .data$decision == "Do not reject" &
                                          .data$k == j)$`f(x,m|pi)`)
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

fisher_pmf_one_stage        <- function(pi, n0, n1, e) {
  x                                        <- expand.grid(0:n0, 0:n1)
  rows_pmf                                 <- (n0 + 1)*(n1 + 1)
  rows_pi                                  <- nrow(pi)
  rows_total                               <- rows_pmf*rows_pi
  f                                        <- numeric(rows_total)
  for (i in 1:rows_pi) {
    dbinom0                                <- stats::dbinom(0:n0, n0, pi[i, 1])
    if (all(n0 == n1, pi[i, 1] == pi[i, 2])) {
      dbinom1                              <- dbinom0
    } else {
      dbinom1                              <- stats::dbinom(0:n1, n1, pi[i, 2])
    }
    f[(1 + (i - 1)*rows_pmf):(i*rows_pmf)] <-
      dbinom0[x[, 1] + 1]*dbinom1[x[, 2] + 1]
  }
  pmf                                      <-
    tibble::tibble(pi0         = rep(pi[, 1], each = rows_pmf),
                   pi1         = rep(pi[, 2], each = rows_pmf),
                   x0          = rep(as.vector(x[, 1]), rows_pi),
                   x1          = rep(as.vector(x[, 2]), rows_pi),
                   m0          = rep(as.integer(n0), rows_total),
                   m1          = rep(as.integer(n1), rows_total),
                   z           = .data$x0 + .data$x1,
                   statistic   = .data$x1 - .data$x0,
                   decision    = ifelse(.data$statistic >=
                                          e[.data$x0 + .data$x1 + 1],
                                        "Reject", "Do not reject"),
                   k           = factor(rep(1, rows_total), 1),
                   `f(x,m|pi)` = f)
  dplyr::arrange(pmf, .data$pi0, .data$pi1, .data$x0, .data$x1)
}

fisher_pmf_two_stage        <- function(pi, n0, n1, e1, f1, e2, k) {
  for (z in which(e1 == Inf)) {
    e1[z]                                        <- z
  }
  for (z in which(f1 == -Inf)) {
    f1[z]                                        <- -z
  }
  for (z1 in 1:(n0[1] + n1[1] + 1)) {
    for (z2 in which(e2[z1, ] == Inf)) {
      e2[z1, z2]                                 <- z1 + z2
    }
  }
  pmf                                            <-
    fisher_pmf_two_stage_cpp(pi[1, ], n0, n1, e1, f1, e2, k)
  rows_pmf                                       <- nrow(pmf)
  rows_pi                                        <- nrow(pi)
  if (rows_pi > 1) {
    pmf                                          <-
      rbind(pmf, matrix(0, (rows_pi - 1)*rows_pmf, 12))
    for (i in 2:rows_pi) {
      pmf[(1 + (i - 1)*rows_pmf):(i*rows_pmf), ] <-
        fisher_pmf_two_stage_cpp(pi[i, ], n0, n1, e1, f1, e2, k)
    }
  }
  pmf                                            <-
    tibble::tibble(pi0         = rep(pi[, 1], each = rows_pmf),
                   pi1         = rep(pi[, 2], each = rows_pmf),
                   x01         = as.integer(pmf[, 1]),
                   x11         = as.integer(pmf[, 2]),
                   x02         = as.integer(pmf[, 3]),
                   x12         = as.integer(pmf[, 4]),
                   m0          = as.integer(pmf[, 5]),
                   m1          = as.integer(pmf[, 6]),
                   z1          = as.integer(pmf[, 7]),
                   z2          = as.integer(pmf[, 8]),
                   statistic   = as.integer(pmf[, 9]),
                   decision    = ifelse(pmf[, 10] == 1, "Reject",
                                        "Do not reject"),
                   k           = factor(pmf[, 11], k),
                   `f(x,m|pi)` = pmf[, 12])
  if (2 %in% k) {
    rows                                         <- which(pmf$k == 1)
    pmf$x02[rows]                                <- pmf$x12[rows] <-
                                                    pmf$z2[rows]  <- NA_integer_
  }
  dplyr::arrange(pmf, .data$pi0, .data$pi1, .data$k, .data$x01, .data$x11,
                 .data$x02, .data$x12)
}

fisher_terminal_one_stage   <- function(n0, n1, e) {
  x        <- expand.grid(0:n0, 0:n1)
  rows_pmf <- (n0 + 1)*(n1 + 1)
  tibble::tibble(x0        = x[, 1],
                 x1        = x[, 2],
                 m0        = rep(as.integer(n0), rows_pmf),
                 m1        = rep(as.integer(n1), rows_pmf),
                 z         = .data$x0 + .data$x1,
                 statistic = .data$x1 - .data$x0,
                 decision  = ifelse(.data$statistic >=
                                      e[.data$x0 + .data$x1 + 1],
                                    "Reject", "Do not reject"),
                 k         = factor(rep(1, rows_pmf), 1))
}

fisher_terminal_two_stage   <- function(n0, n1, e1, f1, e2, k) {
  for (z in which(e1 == Inf)) {
    e1[z]              <- z
  }
  for (z in which(f1 == -Inf)) {
    f1[z]              <- -z
  }
  for (z1 in 1:(n0[1] + n1[1] + 1)) {
    for (z2 in which(e2[z1, ] == Inf)) {
      e2[z1, z2]       <- z1 + z2
    }
  }
  terminal             <- fisher_terminal_two_stage_cpp(n0, n1, e1, e2, f1, k)
  terminal             <- tibble::tibble(x01       = as.integer(terminal[, 1]),
                                         x11       = as.integer(terminal[, 2]),
                                         x02       = as.integer(terminal[, 3]),
                                         x12       = as.integer(terminal[, 4]),
                                         m0        = as.integer(terminal[, 5]),
                                         m1        = as.integer(terminal[, 6]),
                                         z1        = as.integer(terminal[, 7]),
                                         z2        = as.integer(terminal[, 8]),
                                         statistic = as.integer(terminal[, 9]),
                                         decision  = ifelse(terminal[, 10] == 1,
                                                            "Reject",
                                                            "Do not reject"),
                                         k         = factor(terminal[, 11], k))
  if (2 %in% k) {
    rows               <- which(terminal$k == 1)
    terminal$x02[rows] <- terminal$x12[rows] <- terminal$z2[rows] <- NA_integer_
  }
  dplyr::arrange(terminal, .data$k, .data$x01, .data$x11, .data$x02, .data$x12)
}
fisher_des_one_stage        <- function(alpha, beta, delta, ratio, Pi0, Pi1,
                                        nCmax, summary) {
  params                  <- search_parameters(1, "fisher", nCmax, ratio)
  feasible                <-
    fisher_des_one_stage_cpp(alpha, beta, delta, params$poss_nC,
                             params$poss_nE, params$poss_x, params$poss_y,
                             params$poss_z, params$choose_mat,
                             (length(Pi0) == 1), Pi0, (length(Pi1) == 1), Pi1,
                             summary)
  if (feasible[1, 1] > 0) {
    if (summary) {
      message("..feasible designs identified in the range of ",
              "considered sample sizes.\n  Identifying the optimal design..")
    }
    ncol_feasible         <- ncol(feasible)
    nrow_feasible         <- nrow(feasible)
    if (nrow_feasible == 1) {
      feasible_e1         <- matrix(c(1, feasible[, 2:(ncol_feasible - 4)]), 1)
      feasible            <- matrix(c(1, feasible[, 1],
                                      params$poss_nE_orig[feasible[, 1]],
                                      feasible[, (ncol_feasible - 3):
                                                 ncol_feasible]), 1)
    } else {
      feasible_e1         <- cbind(1:nrow_feasible,
                                   feasible[, 2:(ncol_feasible - 4)])
      feasible            <-
        cbind(1:nrow_feasible, feasible[, 1],
              params$poss_nE_orig[feasible[, 1]],
              feasible[, (ncol_feasible - 3):ncol_feasible])
    }
    colnames(feasible)    <- c("index", "n1C", "n1E", "argmax alpha",
                               "max alpha", "argmin power", "min power")
    feasible              <- tibble::as_tibble(feasible)
    feasible[, 1:3]       <- dplyr::mutate_if(feasible[, 1:3], is.double,
                                              as.integer)
    feasible_e1           <-
      feasible_e1[, 1:(max(rowSums(feasible[, 2:3])) + 2)]
    colnames(feasible_e1) <- c("index",
                              paste0("e1", 0:max(rowSums(feasible[, 2:3]))))
    feasible_e1           <- tibble::as_tibble(feasible_e1)
    feasible_e1           <- dplyr::mutate_if(feasible_e1, is.double,
                                              as.integer)
    e1                    <-
      as.integer(feasible_e1[1, 2:(feasible$n1C[1] + feasible$n1E[1] + 2)])
    nC                    <- feasible$n1C[1]
    nE                    <- feasible$n1E[1]
    opchar                <-
      fisher_opchar_one_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                    feasible$`argmin power`[1] + c(0, delta)),
                              nC, nE, e1)
  } else {
    e1                    <- feasible <- nC <- nE <- opchar <- NULL
    if (summary) {
      message("..no feasible designs identified in range of ",
              "considered maximal allowed sample size.\n  Consider increasing ",
              "nCmax...")
    }
  }
  build_des_one_stage_output(alpha, beta, list(e1 = e1, f1 = e1), delta,
                             feasible, nC, nCmax, nE, opchar, Pi0, Pi1, ratio,
                             summary, "fisher")
}

fisher_des_two_stage        <- function(alpha, beta, delta, ratio, Pi0, Pi1,
                                        nCmax, equal, w, piO, efficacy_type,
                                        efficacy_param, futility_type,
                                        futility_param, summary) {
  params                         <- search_parameters(2, "fisher", nCmax, ratio)
  feasible                       <-
    fisher_des_two_stage_cpp(alpha, beta, delta, params$poss_nC, params$poss_nE,
                             params$poss_x, params$poss_y, params$poss_z,
                             params$choose_mat, (length(Pi0) == 1), Pi0,
                             (length(Pi1) == 1), Pi1, equal, efficacy_type,
                             efficacy_param, futility_type, futility_param, Pi0,
                             summary)
  if (feasible[[1]][1, 1] > 0) {
    if (summary) {
      message("..feasible designs identified in the range of ",
              "considered sample sizes.\n  Identifying the optimal design..")
    }
    nrow_feasible                <- nrow(feasible[[1]])
    feasible_e2_vec              <- feasible[[4]]
    if (nrow_feasible == 1) {
      feasible_e1                <- matrix(c(1, feasible[[2]]), 1)
      feasible_f1                <- matrix(c(1, feasible[[3]]), 1)
      feasible                   <- feasible[[1]]
      feasible                   <- matrix(c(1:nrow_feasible, feasible[, 1:2],
                                             params$poss_nE_orig[feasible[, 1]],
                                             params$poss_nE_orig[feasible[, 2]],
                                             feasible[, -(1:2)]), 1)
    } else {
      feasible_e1                <- cbind(1:nrow_feasible, feasible[[2]])
      feasible_f1                <- cbind(1:nrow_feasible, feasible[[3]])
      feasible                   <- feasible[[1]]
      feasible                   <- cbind(1:nrow_feasible, feasible[, 1:2],
                                          params$poss_nE_orig[feasible[, 1]],
                                          params$poss_nE_orig[feasible[, 2]],
                                          feasible[, -(1:2)])
    }
    colnames(feasible)           <- c("index", "n1C", "n2C", "n1E", "n2E",
                                      "argmax alpha", "max alpha",
                                      "argmin power", "min power",
                                      "ESS(piO,piO)", "ESS(piO,piO+delta)")
    feasible                     <- tibble::as_tibble(feasible)
    feasible[, 1:5]              <- dplyr::mutate_if(feasible[, 1:5], is.double,
                                                     as.integer)
    feasible_e1                  <-
      feasible_e1[, 1:(max(rowSums(feasible[, c(2, 4)])) + 2)]
    colnames(feasible_e1)        <-
      c("index", paste0("e1", 0:max(rowSums(feasible[, c(2, 4)]))))
    feasible_e1                  <- tibble::as_tibble(feasible_e1)
    feasible_f1                  <-
      feasible_f1[, 1:(max(rowSums(feasible[, c(2, 4)])) + 2)]
    colnames(feasible_f1)        <-
      c("index", paste0("f1", 0:max(rowSums(feasible[, c(2, 4)]))))
    feasible_f1                  <- tibble::as_tibble(feasible_f1)
    feasible_e2                  <- list()
    for (i in 1:nrow_feasible) {
      nC                         <- as.numeric(feasible[i, 2:3])
      nE                         <- as.numeric(feasible[i, 4:5])
      e2_i                       <- matrix(0L, nC[1] + nE[1] + 1,
                                           nC[2] + nE[2] + 1)
      counter                    <- 1L
      for (z1p in 1:(max(params$poss_nC) + max(params$poss_nE) + 1)) {
        for (z2p in 1:(2*(max(params$poss_nC) + max(params$poss_nE)) + 1)) {
          if (all(z1p <= nC[1] + nE[1] + 1, z2p <= nC[2] + nE[2] + 1)) {
            if (feasible_e2_vec[i, counter] != -0.5) {
              e2_i[z1p, z2p]     <- feasible_e2_vec[i, counter]
            } else if (feasible_e2_vec[i, counter] == z1p + z2p - 1) {
              e2_i[z1p, z2p]     <- Inf
            } else {
              e2_i[z1p, z2p]     <- NA_integer_
            }
          }
          counter                <- counter + 1L
        }
      }
      feasible_e2[[i]]           <- e2_i
    }
    feasible                     <-
      dplyr::mutate(feasible,
                    `argmax ESS(pi,pi)`       = 0,
                    `max ESS(pi,pi)`          = 0,
                    `argmax_piC ESS(piC,piE)` = 0,
                    `argmax_piE ESS(piC,piE)` = 0,
                    `max ESS(piC,piE)`        = 0,
                    `max N`                   =
                      as.integer(rowSums(feasible[, 2:5])))
    for (i in 1:nrow_feasible) {
      index                      <-
        as.numeric(feasible[i, 2] + params$max_poss_nC*(feasible[i, 4] - 1))
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
      feasible$`argmax_piC ESS(piC,piE)`[i] <- max_ESS_2d[1]
      feasible$`argmax_piE ESS(piC,piE)`[i] <- max_ESS_2d[2]
      feasible$`max ESS(piC,piE)`[i]        <- max_ESS_2d[3]
    }
    feasible$o                   <-
      rowSums(matrix(w, nrow_feasible, 5, TRUE)*feasible[, c(10:11, 13, 16:17)])
    feasible                     <-
      dplyr::arrange(feasible, .data$o, dplyr::desc(.data$`min power`))
    feasible_e1                  <- feasible_e1[as.matrix(feasible[, 1]), ]
    feasible_f1                  <- feasible_f1[as.matrix(feasible[, 1]), ]
    ncol_feasible_e1             <- ncol(feasible_e1)
    for (i in 1:nrow_feasible) {
      n1C                        <- feasible$n1C[i]
      n1E                        <- feasible$n1E[i]
      for (z1 in 0:(n1C + n1E)) {
        if (feasible_e1[i, z1 + 2] >= z1 + 1) {
          feasible_e1[i, z1 + 2] <- Inf
        }
        if (feasible_f1[i, z1 + 2] <= -z1 - 1) {
          feasible_f1[i, z1 + 2] <- -Inf
        }
      }
      if (n1C + n1E + 2 < ncol_feasible_e1) {
        range                    <- (n1C + n1E + 3):ncol_feasible_e1
        feasible_e1[i, range]    <- feasible_f1[i, range] <- NA_integer_
      }
    }
    e1                           <-
      as.numeric(feasible_e1[feasible$index[1],
                             2:(feasible$n1C[1] + feasible$n1E[1] + 2)])
    e2                           <- feasible_e2[[feasible$index[1]]]
    f1                           <-
      as.numeric(feasible_f1[feasible$index[1],
                             2:(feasible$n1C[1] + feasible$n1E[1] + 2)])
    nC                           <- c(feasible$n1C[1], feasible$n2C[1])
    nE                           <- c(feasible$n1E[1], feasible$n2E[1])
    opchar                       <-
      fisher_opchar_two_stage(rbind(rep(feasible$`argmax alpha`[1], 2),
                                    feasible$`argmin power`[1] + c(0, delta)),
                              nC, nE, e1, f1, e2, 1:2)
  } else {
    e1                           <-          e2 <-          f1 <-    feasible <-
                                    feasible_e1 <- feasible_e2 <- feasible_f1 <-
                                             nC <-          nE <-      opchar <-
                                           NULL
    if (summary) {
      message("..no feasible designs identified in range of ",
              "considered maximal allowed sample size.\n  Consider increasing ",
              "nCmax..")
    }
  }
  build_des_two_stage_output(alpha, beta,
                             list(e1 = e1, e2 = e2, f1 = f1, f2 = e2), delta,
                             NULL, efficacy_param, efficacy_type, equal,
                             feasible, NULL, futility_param, futility_type, nC,
                             nCmax, nE, opchar, Pi0, Pi1, piO, ratio, summary,
                             w, "fisher")
}

fisher_max_ess_2d_two_stage <- function(nC, nE, e1, f1, poss_x1, poss_y1,
                                        poss_z1) {
  max_ESS_2d <- stats::optim(par     = c(0.5, 0.5),
                             fn      = fisher_minus_ess_two_stage,
                             nC      = nC,
                             nE      = nE,
                             e1      = e1,
                             f1      = f1,
                             poss_x1 = poss_x1,
                             poss_y1 = poss_y1,
                             poss_z1 = poss_z1)
  c(max_ESS_2d$par, -max_ESS_2d$value)
}

fisher_minus_ess_two_stage  <- function(pi, nC, nE, e1, f1, poss_x1, poss_y1,
                                        poss_z1) {
  if (!missing(poss_x1)) {
    -fisher_des_ess_two_stage(pi, nC, nE, e1, f1, poss_x1, poss_y1, poss_z1)
  } else {
    -fisher_ess_two_stage(pi, nC, nE, e1, f1)
  }
}

fisher_opchar_one_stage     <- function(pi, nC, nE, e, pmf_pi) {
  if (missing(pmf_pi)) {
    pmf_pi <- fisher_pmf_one_stage(pi, nC, nE, e)
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

fisher_opchar_two_stage     <- function(pi, nC, nE, e1, f1, e2, k, pmf_pi) {
  nrow_pi          <- nrow(pi)
  n                <- c(nC[1] + nE[1], sum(nC) + sum(nE))
  opchar           <- matrix(0, nrow_pi, 13)
  E                <- Fu <- numeric(2)
  if (missing(pmf_pi)) {
    for (i in 1:nrow_pi) {
      pmf_pi         <- fisher_pmf_two_stage(pi[i, , drop = FALSE], nC, nE, e1,
                                             f1, e2, k)
      for (j in k) {
        E[j]         <- sum(dplyr::filter(pmf_pi, .data$decision == "Reject" &
                                            .data$k == j)$`f(x,m|pi)`)
        Fu[j]        <- sum(dplyr::filter(pmf_pi, .data$decision == "Do not reject" &
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
  } else {
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
  }
  colnames(opchar) <- c("piC", "piE", "P(pi)", "ESS(pi)", "SDSS(pi)", "MSS(pi)",
                        paste0(rep(c("E", "F", "S"), each = 2), rep(1:2, 3),
                              "(pi)"), "max N")
  opchar           <- tibble::as_tibble(opchar)
  opchar$`max N`   <- as.integer(opchar$`max N`)
  opchar
}

fisher_pmf_one_stage        <- function(pi, nC, nE, e1) {
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
  pmf                                      <-
    tibble::tibble(piC         = rep(pi[, 1], each = nrow_pmf),
                   piE         = rep(pi[, 2], each = nrow_pmf),
                   xC          = rep(as.vector(x[, 1]), nrow_pi),
                   xE          = rep(as.vector(x[, 2]), nrow_pi),
                   mC          = rep(as.integer(nC), nrow_total),
                   mE          = rep(as.integer(nE), nrow_total),
                   z           = .data$xC + .data$xE,
                   statistic   = .data$xE - .data$xC,
                   decision    = factor(ifelse(.data$statistic >=
                                                 e1[.data$xC + .data$xE + 1],
                                               "Reject", "Do not reject")),
                   k           = factor(rep(1, nrow_total), 1),
                   `f(x,m|pi)` = f)
  dplyr::arrange(pmf, .data$piC, .data$piE, .data$xC, .data$xE)
}

fisher_pmf_two_stage        <- function(pi, nC, nE, e1, f1, e2, k) {
  for (z in which(e1 == Inf)) {
    e1[z]                                        <- z
  }
  for (z in which(f1 == -Inf)) {
    f1[z]                                        <- -z
  }
  for (z1 in 1:(nC[1] + nE[1] + 1)) {
    for (z2 in which(e2[z1, ] == Inf)) {
      e2[z1, z2]                                 <- z1 + z2
    }
  }
  pmf                                            <-
    fisher_pmf_two_stage_cpp(pi[1, ], nC, nE, e1, f1, e2, k)
  nrow_pmf                                       <- nrow(pmf)
  nrow_pi                                        <- nrow(pi)
  if (nrow_pi > 1) {
    pmf                                          <-
      rbind(pmf, matrix(0, (nrow_pi - 1)*nrow_pmf, 12))
    for (i in 2:nrow_pi) {
      pmf[(1 + (i - 1)*nrow_pmf):(i*nrow_pmf), ] <-
        fisher_pmf_two_stage_cpp(pi[i, ], nC, nE, e1, f1, e2, k)
    }
  }
  pmf                                            <-
    tibble::tibble(piC         = rep(pi[, 1], each = nrow_pmf),
                   piE         = rep(pi[, 2], each = nrow_pmf),
                   xC1         = as.integer(pmf[, 1]),
                   xE1         = as.integer(pmf[, 2]),
                   xC2         = as.integer(pmf[, 3]),
                   xE2         = as.integer(pmf[, 4]),
                   mC          = as.integer(pmf[, 5]),
                   mE          = as.integer(pmf[, 6]),
                   z1          = as.integer(pmf[, 7]),
                   z2          = as.integer(pmf[, 8]),
                   statistic   = as.integer(pmf[, 9]),
                   decision    = factor(c("Do not reject",
                                          "Reject")[pmf[, 10] + 1]),
                   k           = factor(pmf[, 11], k),
                   `f(x,m|pi)` = pmf[, 12])
  if (2 %in% k) {
    rows                                         <- which(pmf$k == 1)
    pmf$xC2[rows]                                <- pmf$xE2[rows] <-
                                                    pmf$z2[rows]  <- NA_integer_
  }
  dplyr::arrange(pmf, .data$piC, .data$piE, .data$k, .data$xC1, .data$xE1,
                 .data$xC2, .data$xE2)
}

fisher_terminal_one_stage   <- function(nC, nE, e1) {
  x        <- expand.grid(0:nC, 0:nE)
  nrow_pmf <- (nC + 1)*(nE + 1)
  tibble::tibble(xC        = x[, 1],
                 xE        = x[, 2],
                 mC        = rep(as.integer(nC), nrow_pmf),
                 mE        = rep(as.integer(nE), nrow_pmf),
                 z         = .data$xC + .data$xE,
                 statistic = .data$xE - .data$xC,
                 decision  = factor(ifelse(.data$statistic >=
                                             e1[.data$xC + .data$xE + 1],
                                           "Reject", "Do not reject")),
                 k         = factor(rep(1, nrow_pmf), 1))
}

fisher_terminal_two_stage   <- function(nC, nE, e1, f1, e2, k) {
  for (z in which(e1 == Inf)) {
    e1[z]              <- z
  }
  for (z in which(f1 == -Inf)) {
    f1[z]              <- -z
  }
  for (z1 in 1:(nC[1] + nE[1] + 1)) {
    for (z2 in which(e2[z1, ] == Inf)) {
      e2[z1, z2]       <- z1 + z2
    }
  }
  terminal             <- fisher_terminal_two_stage_cpp(nC, nE, e1, e2, f1, k)
  terminal             <-
    tibble::tibble(xC1       = as.integer(terminal[, 1]),
                   xE1       = as.integer(terminal[, 2]),
                   xC2       = as.integer(terminal[, 3]),
                   xE2       = as.integer(terminal[, 4]),
                   mC        = as.integer(terminal[, 5]),
                   mE        = as.integer(terminal[, 6]),
                   z1        = as.integer(terminal[, 7]),
                   z2        = as.integer(terminal[, 8]),
                   statistic = as.integer(terminal[, 9]),
                   decision  =
                     factor(c("Do not reject", "Reject",
                              "Continue to stage 2")[terminal[, 10]],
                            sort(unique(c("Do not reject", "Reject",
                                          "Continue to stage 2")[terminal[, 10]]))),
                   k         = factor(terminal[, 11], k))
  if (2 %in% k) {
    rows               <- which(terminal$k == 1)
    terminal$xC2[rows] <- terminal$xE2[rows] <- terminal$z2[rows] <- NA_integer_
  }
  dplyr::arrange(terminal, .data$k, .data$xC1, .data$xE1, .data$xC2, .data$xE2)
}
summary_des         <- function(J, type, alpha, beta, delta, ratio, Pi0, Pi1,
                                nCmax, equal, w, piO, efficacy, futility,
                                efficacy_type, efficacy_param, futility_type,
                                futility_param) {
  if (J == 1) {
    stage  <- "one-stage"
  } else {
    stage  <- "two-stage"
  }
  if (type == "barnard") {
    dashes <- 57
    design <- "Barnard's exact test"
  } else if (type == "binomial") {
    dashes <- 59
    design <- "an exact binomial test"
  } else if (type == "fisher") {
    dashes <- 56
    design <- "Fisher's exact test"
  } else {
    dashes <- 74
    design <- "one-arm and two-arm testing decisions"
  }
  message("  ", rep("-", dashes))
  message("  Design of a ", stage, " trial based on ", design)
  message("  ", rep("-", dashes))
  message("\n  ---------------")
  message("  Hypothesis test")
  message("  ---------------")
  message("  You have chosen to test the following hypothesis")
  message("    H0 : piE <= piC")
  message("  with the following type-I error constraint")
  if (length(Pi0) == 1) {
    message("    P(", Pi0, ",", Pi0, ") <= alpha = ", alpha)
  } else {
    message("    max_{pi in Pi0} P(pi,pi) <= alpha = ", alpha, ", Pi0 = [",
            Pi0[1], ",", Pi0[2], "]")
  }
  message("  and the following type-II error constraint")
  if (length(Pi1) == 1) {
    message("    P(", Pi1, ",", Pi1 + delta, ") >= 1 - beta = ", 1 - beta)
  } else {
    message("    max_{pi in Pi1} P(pi,pi+delta) >= 1 - beta = ", 1 - beta,
            ", Pi1 = [", Pi1[1], ",", Pi1[2], "], delta = ", delta)
  }
  message("\n  ------------")
  message("  Restrictions")
  message("  ------------")
  message("  - You have chosen to limit the allowed maximal sample ",
          "size in the control arm, nC,\n    to: ", nCmax)
  if (J == 1) {
    message("  - The sample size in the experimental arm, nE, will be set",
            " to: r x nC, with r = ", ratio)
  } else {
    if (equal) {
      message("  - You have chosen to restrict the sample sizes in the",
              " control arm in each stage, n1C\n    and n2C, such that: n1C = ",
              "n2C")
    } else {
      message("  - You have chosen to allow the sample sizes in the",
              " control arm in each stage, n1C\n    and n2C, to take different",
              " values")
    }
    message("  - The sample sizes in the experimental arm in each ",
            "stage, n1E and n2E, will be se\n    to: r x n1C and r x n2C ",
            "respectively, with r = ", ratio)
    if (type == "fisher") {
      if (efficacy_type == 0) {
        message("  - You have chosen to prevent early stopping for ",
                "efficacy. Thus e1z1 = Inf, for all z1, in all",
                "\n    considered designs")
      } else if (efficacy_type == 1) {
        if (efficacy_param == -0.5) {
          message("  - You have chosen to include early stopping for ",
                  "efficacy, with e1z1 = [0.5*(n1C + n1E)*delta]_* + 1, for ",
                  "all z1, in all considered designs")
        } else {
          message("  - You have chosen to include early stopping for ",
                  "efficacy, with e1z1 = ", efficacy_param, ", for all z1, in ",
                  "all considered designs.")
        }
      } else {
        message("  - You have chosen to include early stopping for ",
                "efficacy, with e1z1 chosen for each z1, in each considered ",
                "design, to control the probability of committing a type-I ",
                "error at the end of stage one to ", efficacy_param)
      }
      if (futility_type == 0) {
        message("  - You have chosen to prevent early stopping for ",
                "futility. Thus f1z1 = -Inf, for all z1, in all considered ",
                "designs")
      } else if (futility_type == 1) {
        message("  - You have chosen to include early stopping for ",
                "futility, with f1z1 = ", futility_param, ", for all z1, in ",
                "all considered designs.")
      } else {
        message("  - You have chosen to include early stopping for ",
                "futility, with f1z1 chosen for each z1, in each considered ",
                "design, to control the probability of committing a type-II ",
                "error at the end of stage one to ", futility_param)
      }
    } else {
      if (futility) {
        message("  - You have chosen to allow early stopping for ",
                "futility")
      } else {
        if (type != "sat") {
          message("  - You have chosen to prevent early stopping for ",
                  "futility. Thus f1 = -Inf in all considered designs")
        } else {
          message("  - You have chosen to prevent early stopping for ",
                  "futility. Thus fS1 = fT1 = -inf in all considered ",
                  "designs")
        }
      }
      if (efficacy) {
        message("  - You have chosen to allow early stopping for ",
                "efficacy")
      } else {
        if (type != "sat") {
          message("  - You have chosen to prevent early stopping for ",
                  "efficacy. Thus e1 = Inf in all considered\n    designs")
        } else {
          message("  - You have chosen to prevent early stopping for ",
                  "efficacy. Thus eS1 = eT1 = Inf in all\n    considered ",
                  "designs")
        }
      }
    }
    message("\n  The design will be optimised for:")
    message("    w1ESS(piO,piO) + w2ESS(piO,piO + delta) + w3max_pi ESS(pi,",
            "pi)\n    + w4max_{piC,piE} ESS(piC,piE) + w5max N")
    message("  with:")
    message("    w1 = ", w[1], ", w2 = ", w[2], ", w3 = ", w[3], ", w4 = ",
            w[4], ", w5 = ", w[5])
    message("  and piO = ", piO)
  }
}

summary_opchar      <- function(des, pi, k) {
  if (des$J == 1) {
    stage  <- "one-stage"
  } else {
    stage  <- "two-stage"
  }
  if (des$type == "barnard") {
    design <- "barnard's exact test"
    dashes <- 77
  } else if (des$type == "binomial") {
    design <- "an exact binomial test"
    dashes <- 79
  } else if (des$type == "fisher") {
    design <- "Fisher's exact test"
    dashes <- 76
  } else if (des$type == "sat") {
    design <- "one-arm and two-arm testing\n  decisions"
    dashes <- 84
  }
  message("  ", rep("-", dashes))
  message("  Operating characteristics of a ", stage, " design based on ",
          design)
  message("  ", rep("-", dashes))
  message("\n  You have chosen to analytically determine the operating ",
          "characteristics of a design with")
  if (des$J == 1) {
    message("    - n1C = ", des$nC)
    message("    - n1E = ", des$nE)
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(des$boundaries$e1, 3))
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
    }
  } else if (des$J == 2) {
    message("    - n1C = ", des$nC[1])
    message("    - n2C = ", des$nC[2])
    message("    - n1E = ", des$nE[1])
    message("    - n2E = ", des$nE[2])
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(des$boundaries$e1, 3))
      message("    - f1 = ", round(des$boundaries$f1, 3))
      message("    - e2 = ", round(des$boundaries$e2, 3))
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
      message("    - f10 = ", des$boundaries$f1[1], ", ..., f11 = ",
              des$boundaries$f1[2], ", ..., f1", des$nC + des$nE, " = ",
              des$boundaries$f1[des$nC + des$nE + 1])
      message("    - e200 = ", des$boundaries$e2[1, 1], ", ..., ",
              "e2", des$nC[1] + des$nE[1], des$nC[2] + des$nE[2], " = ",
              des$boundaries$e2[des$nC[1] + des$nE[1], des$nC[2] + des$nE[2]])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
      message("    - fS1 = ", des$boundaries$fS1)
      message("    - fT1 = ", des$boundaries$fT1)
      message("    - eS2 = ", des$boundaries$eS2)
      message("    - eT2 = ", des$boundaries$eT2)
    }
  }
  if (nrow(pi) == 1) {
    message("  when pi = (", pi[1, 1], ", ", pi[1, 2], ")'.")
  } else if (nrow(pi) == 2) {
    message("  when pi in {(", pi[1, 1], ", ", pi[1, 2], ")', (", pi[2, 1],
            ", ", pi[2, 2], ")'}.")
  } else {
    message("  when pi in {(", pi[1, 1], ", ", pi[1, 2], ")', ..., (",
            pi[nrow(pi), 1], ", ", pi[nrow(pi), 2], ")'}.")
  }
}

summary_ph2rand_des <- function(x) {
  J          <- x$J
  type       <- x$type
  alpha      <- x$alpha
  beta       <- x$beta
  delta      <- x$delta
  Pi0        <- x$Pi0
  Pi1        <- x$Pi1
  if (J == 1) {
    stage    <- "single-stage"
  } else {
    stage    <- "two-stage"
  }
  if (type == "barnard") {
    dashes   <- 47
    design   <- "Barnard's exact test"
  } else if (type == "binomial") {
    dashes   <- 49
    design   <- "an exact binomial test"
  } else if (type == "fisher") {
    dashes   <- 46
    design   <- "Fisher's exact test"
  } else {
    dashes   <- 57
    design <- "one-arm and two-arm testing decisions"
  }
  message("  ", rep("-", dashes))
  message("  A ", stage, " trial based on ", design)
  message("  ", rep("-", dashes))
  message("\n  ---------------")
  message("  Hypothesis test")
  message("  ---------------")
  message("  You have chosen to test the following hypothesis")
  message("    H0 : piE <= piC")
  message("  with the following type-I error constraint")
  if (length(Pi0) == 1) {
    message("    P(", Pi0, ",", Pi0, ") <= alpha = ", alpha)
  } else {
    message("    max_{pi in Pi0} P(pi,pi) <= alpha = ", alpha, ", Pi0 = [",
            Pi0[1], ",", Pi0[2], "]")
  }
  message("  and the following type-II error constraint")
  if (length(Pi1) == 1) {
    message("    P(", Pi1, ",", Pi1 + delta, ") >= 1 - beta = ", 1 - beta)
  } else {
    message("    max_{pi in Pi1} P(pi,pi+delta) >= 1 - beta = ", 1 - beta,
            ", Pi1 = [", Pi1[1], ",", Pi1[2], "], delta = ", delta)
  }
  message("\n  -----------------")
  message("  Design parameters")
  message("  -----------------")
  message("  The design has:")
  if (x$J == 1) {
    message("    - n1C = ", x$nC)
    message("    - n1E = ", x$nE)
    if (x$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(x$boundaries$e1, 3))
    } else if (x$type == "fisher") {
      message("    - e10 = ", x$boundaries$e1[1], ", ..., e11 = ",
              x$boundaries$e1[2], ", ..., e1", x$nC + x$nE, " = ",
              x$boundaries$e1[x$nC + x$nE + 1])
    } else if (x$type == "sat") {
      message("    - eS1 = ", x$boundaries$eS1)
      message("    - eT1 = ", x$boundaries$eT1)
    }
  } else if (x$J == 2) {
    message("    - n1C = ", x$nC[1])
    message("    - n2C = ", x$nC[2])
    message("    - n1E = ", x$nE[1])
    message("    - n2E = ", x$nE[2])
    if (x$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(x$boundaries$e1, 3))
      message("    - f1 = ", round(x$boundaries$f1, 3))
      message("    - e2 = ", round(x$boundaries$e2, 3))
    } else if (x$type == "fisher") {
      message("    - e10 = ", x$boundaries$e1[1], ", ..., e11 = ",
              x$boundaries$e1[2], ", ..., e1", x$nC + x$nE, " = ",
              x$boundaries$e1[x$nC + x$nE + 1])
      message("    - f10 = ", x$boundaries$f1[1], ", ..., f11 = ",
              x$boundaries$f1[2], ", ..., f1", x$nC + x$nE, " = ",
              x$boundaries$f1[x$nC + x$nE + 1])
      message("    - e200 = ", x$boundaries$e2[1, 1], ", ..., ",
              "e2", x$nC[1] + x$nE[1], x$nC[2] + x$nE[2], " = ",
              x$boundaries$e2[x$nC[1] + x$nE[1], x$nC[2] + x$nE[2]])
    } else if (x$type == "sat") {
      message("    - eS1 = ", x$boundaries$eS1)
      message("    - eT1 = ", x$boundaries$eT1)
      message("    - fS1 = ", x$boundaries$fS1)
      message("    - fT1 = ", x$boundaries$fT1)
      message("    - eS2 = ", x$boundaries$eS2)
      message("    - eT2 = ", x$boundaries$eT2)
    }
  }
  message("\n  -------------------------")
  message("  Operating Characteristics")
  message("  -------------------------")
  message("  Key operating characteristics include")
  if (x$J == 1) {
    print(x$opchar)
  } else {
    print(x$opchar[, -(11:13)])
  }
}

summary_pmf         <- function(des, pi, k) {
  if (des$J == 1) {
    stage  <- "one-stage"
  } else {
    stage  <- "two-stage"
  }
  if (des$type == "barnard") {
    design <- "barnard's exact test"
    dashes <- 55
  } else if (des$type == "binomial") {
    design <- "an exact binomial test"
    dashes <- 57
  } else if (des$type == "fisher") {
    design <- "Fisher's exact test"
    dashes <- 54
  } else if (des$type == "sat") {
    design <- "single-arm and two-arm testing\n  decisions"
    dashes <- 65
  }
  message("  ", rep("-", dashes))
  message("  PMF of a ", stage, " design based on ", design)
  message("  ", rep("-", dashes))
  message("\n  You have chosen to find the PMF of a design with")
  if (des$J == 1) {
    message("    - n1C = ", des$nC)
    message("    - n1E = ", des$nE)
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(des$boundaries$e1, 3))
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
    }
  } else if (des$J == 2) {
    message("    - n1C = ", des$nC[1])
    message("    - n2C = ", des$nC[2])
    message("    - n1E = ", des$nE[1])
    message("    - n2E = ", des$nE[2])
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(des$boundaries$e1, 3))
      message("    - f1 = ", round(des$boundaries$f1, 3))
      message("    - e2 = ", round(des$boundaries$e2, 3))
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
      message("    - f10 = ", des$boundaries$f1[1], ", ..., f11 = ",
              des$boundaries$f1[2], ", ..., f1", des$nC + des$nE, " = ",
              des$boundaries$f1[des$nC + des$nE + 1])
      message("    - e200 = ", des$boundaries$e2[1, 1], ", ..., ",
              "e2", des$nC[1] + des$nE[1], des$nC[2] + des$nE[2], " = ",
              des$boundaries$e2[des$nC[1] + des$nE[1], des$nC[2] + des$nE[2]])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
      message("    - fS1 = ", des$boundaries$fS1)
      message("    - fT1 = ", des$boundaries$fT1)
      message("    - eS2 = ", des$boundaries$eS2)
      message("    - eT2 = ", des$boundaries$eT2)
    }
  }
  if (nrow(pi) == 1) {
    message("  when pi = (", pi[1, 1], ", ", pi[1, 2], ")'.")
  } else if (nrow(pi) == 2) {
    message("  when pi in {(", pi[1, 1], ", ", pi[1, 2], ")', (", pi[2, 1],
            ", ", pi[2, 2], ")'}.")
  } else {
    message("  when pi in {(", pi[1, 1], ", ", pi[1, 2], ")', ..., (",
            pi[nrow(pi), 1], ", ", pi[nrow(pi), 2], ")'}.")
  }
}

summary_sim         <- function(des, pi, k, replicates) {
  if (des$J == 1) {
    stage  <- "one-stage"
  } else {
    stage  <- "two-stage"
  }
  if (des$type == "barnard") {
    design <- "barnard's exact test"
    dashes <- 77
  } else if (des$type == "binomial") {
    design <- "an exact binomial test"
    dashes <- 79
  } else if (des$type == "fisher") {
    design <- "Fisher's exact test"
    dashes <- 76
  } else if (des$type == "sat") {
    design <- "one-arm and two-arm testing\n  decisions"
    dashes <- 84
  }
  message("  ", rep("-", dashes))
  message("  Operating characteristics of a ", stage, " design based on ",
          design)
  message("  ", rep("-", dashes))
  message("\n  You have chosen to estimate via simulation the operating ",
          "characteristics of a design\n  with")
  if (des$J == 1) {
    message("    - n1C = ", des$nC)
    message("    - n1E = ", des$nE)
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(des$boundaries$e1, 3))
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
    }
  } else if (des$J == 2) {
    message("    - n1C = ", des$nC[1])
    message("    - n2C = ", des$nC[2])
    message("    - n1E = ", des$nE[1])
    message("    - n2E = ", des$nE[2])
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", round(des$boundaries$e1, 3))
      message("    - f1 = ", round(des$boundaries$f1, 3))
      message("    - e2 = ", round(des$boundaries$e2, 3))
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
      message("    - f10 = ", des$boundaries$f1[1], ", ..., f11 = ",
              des$boundaries$f1[2], ", ..., f1", des$nC + des$nE, " = ",
              des$boundaries$f1[des$nC + des$nE + 1])
      message("    - e200 = ", des$boundaries$e2[1, 1], ", ..., ",
              "e2", des$nC[1] + des$nE[1], des$nC[2] + des$nE[2], " = ",
              des$boundaries$e2[des$nC[1] + des$nE[1], des$nC[2] + des$nE[2]])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
      message("    - fS1 = ", des$boundaries$fS1)
      message("    - fT1 = ", des$boundaries$fT1)
      message("    - eS2 = ", des$boundaries$eS2)
      message("    - eT2 = ", des$boundaries$eT2)
    }
  }
  if (nrow(pi) == 1) {
    message("  when pi = (", pi[1, 1], ", ", pi[1, 2], ")'.")
  } else if (nrow(pi) == 2) {
    message("  when pi in {(", pi[1, 1], ", ", pi[1, 2], ")', (", pi[2, 1],
            ", ", pi[2, 2], ")'}.")
  } else {
    message("  when pi in {(", pi[1, 1], ", ", pi[1, 2], ")', ..., (",
            pi[nrow(pi), 1], ", ", pi[nrow(pi), 2], ")'}.")
  }
  message("\n  ", replicates, " simulations will be used for each value of pi.")
}

summary_terminal    <- function(des, k) {
  if (des$J == 1) {
    stage  <- "one-stage"
  } else {
    stage  <- "two-stage"
  }
  if (des$type == "barnard") {
    design <- "barnard's exact test"
    dashes <- 67
  } else if (des$type == "binomial") {
    design <- "an exact binomial test"
    dashes <- 69
  } else if (des$type == "fisher") {
    design <- "Fisher's exact test"
    dashes <- 66
  } else if (des$type == "sat") {
    design <- "single-arm and two-arm testing\n  decisions"
    dashes <- 77
  }
  message("  ", rep("-", dashes))
  message("  Terminal points of a ", stage, " based on ", design)
  message("  ", rep("-", dashes))
  message("\n  You have chosen to find the terminal points of a design with")
  if (des$J == 1) {
    message("    - n1C = ", des$nC)
    message("    - n1E = ", des$nE)
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", des$boundaries$e1)
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
    }
  } else if (des$J == 2) {
    message("    - n1C = ", des$nC[1])
    message("    - n2C = ", des$nC[2])
    message("    - n1E = ", des$nE[1])
    message("    - n2E = ", des$nE[2])
    if (des$type %in% c("barnard", "binomial")) {
      message("    - e1 = ", des$boundaries$e1)
      message("    - f1 = ", des$boundaries$f1)
      message("    - e2 = ", des$boundaries$e2)
    } else if (des$type == "fisher") {
      message("    - e10 = ", des$boundaries$e1[1], ", ..., e11 = ",
              des$boundaries$e1[2], ", ..., e1", des$nC + des$nE, " = ",
              des$boundaries$e1[des$nC + des$nE + 1])
      message("    - f10 = ", des$boundaries$f1[1], ", ..., f11 = ",
              des$boundaries$f1[2], ", ..., f1", des$nC + des$nE, " = ",
              des$boundaries$f1[des$nC + des$nE + 1])
      message("    - e200 = ", des$boundaries$e2[1, 1], ", ..., ",
              "e2", des$nC[1] + des$nE[1], des$nC[2] + des$nE[2], " = ",
              des$boundaries$e2[des$nC[1] + des$nE[1], des$nC[2] + des$nE[2]])
    } else if (des$type == "sat") {
      message("    - eS1 = ", des$boundaries$eS1)
      message("    - eT1 = ", des$boundaries$eT1)
      message("    - fS1 = ", des$boundaries$fS1)
      message("    - fT1 = ", des$boundaries$fT1)
      message("    - eS2 = ", des$boundaries$eS2)
      message("    - eT2 = ", des$boundaries$eT2)
    }
  }
}
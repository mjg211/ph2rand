summary_des      <- function(J, type, alpha, beta, delta, ratio, Pi0, Pi1,
                             nCmax, equal, w, piO, efficacy, futility,
                             efficacy_type, efficacy_param, futility_type,
                             futility_param) {
  if (J == 1) {
    stage    <- "single-stage"
  } else {
    stage    <- "two-stage"
  }
  if (type == "barnard") {
    dashes   <- 60 - 3*(J == 2)
    design   <- "Barnard's exact test"
  } else if (type == "binomial") {
    dashes   <- 62 - 3*(J == 2)
    design   <- "an exact binomial test"
  } else if (type == "fisher") {
    dashes   <- 59 - 3*(J == 2)
    design   <- "Fisher's exact test"
  } else {
    dashes   <- 70 + 7*(J == 2)
    if (J == 1) {
      design <- "single-arm and two-arm testing\n  decisions"
    } else {
      design <- "single-arm and two-arm testing decisions"
    }
  }
  message("  ", rep("-", dashes))
  message("  Design of a ", stage, " trial based on ", design)
  message("  ", rep("-", dashes))
  message("\n  ---------------")
  message("  Hypothesis test")
  message("  ---------------")
  message("  You have chosen to test the following hypothesis:")
  if (length(Pi0) == 1) {
    message("    H", uc_sub(0), " : ", uc("pi"), "C = ", uc("pi"), "E \u2208 ",
            uc("Pi"), uc_sub(0), " = ", Pi0)
  } else {
    message("    H", uc_sub(0), " : ", uc("pi"), "C = ", uc("pi"), "E \u2208 ",
            uc("Pi"), uc_sub(0), " = [", Pi0[1], ",", Pi0[2], "]")
  }
  message("  with the following type-I error constraint:")
  if (length(Pi0) == 1) {
    message("    P(", Pi0, ",", Pi0, ") ", uc("le"), " ", uc("alpha"), " = ",
            alpha)
  } else {
    message("    max_{", uc("pi"), " \u2208 ", uc("Pi"), uc_sub(0), "} P(",
            uc("pi"), ",", uc("pi"), ") ", uc("le"), " ", uc("alpha"), " = ",
            alpha, ", ", uc("Pi"), uc_sub(0), " = [", Pi0[1], ",", Pi0[2], "]")
  }
  message("  and the following type-II error constraint:")
  if (length(Pi1) == 1) {
    message("    P(", Pi1, ",", Pi1 + delta, ") ", uc("ge"), " 1 - ",
            uc("beta"), " = ", 1 - beta)
  } else {
    message("    max_{", uc("pi"), " \u2208 ", uc("Pi"), uc_sub(1), "} P(",
            uc("pi"), ",", uc("pi"), " + ", uc("delta"), ") ", uc("ge"),
            " 1 - ", uc("beta"), " = ", 1 - beta, ", ", uc("Pi"), uc_sub(1),
            " = [", Pi1[1], ",", Pi1[2], "], ", uc("delta"), " = ", delta)
  }
  message("\n  ------------")
  message("  Restrictions")
  message("  ------------")
  message("  \u2022 You have chosen to limit the allowed maximal sample ",
          "size in the control arm, nC,\n    to: ", nCmax)
  if (J == 1) {
    message("  \u2022 The sample size in the experimental arm, nE, will be set",
            " to: r \u00D7 nC, with r = ", ratio)
  } else {
    if (equal) {
      message("  \u2022 You have chosen to restrict the sample sizes in the",
              " control arm in each stage, nC1\n    and nC2, such that: nC1 = ",
              "nC2")
    } else {
      message("  \u2022 You have chosen to allow the sample sizes in the",
              " control arm in each stage, nC1\n    and nC2, to take different",
              " values")
    }
    message("  \u2022 The sample sizes in the experimental arm in each ",
            "stage, nE1 and nE2, will be se\n    to: r \u00D7 nC1 and r \u00D7",
            " nC2 respectively, with r = ", ratio)
    if (type == "fisher") {
      if (efficacy_type == 0) {
        message("  \u2022 You have chosen to prevent early stopping for ",
                "efficacy. Thus e1z1 = \u221E, for all z1, in all",
                "\n    considered designs")
      } else if (efficacy_type == 1) {
        if (efficacy_param == -0.5) {
          message("  \u2022 You have chosen to include early stopping for ",
                  "efficacy, with e1z1 = [0.5(n1C + n1E)", uc("delta"),
                  "]_* + 1, for all z1, in all considered designs")
        } else {
          message("  \u2022 You have chosen to include early stopping for ",
                  "efficacy, with e1z1 = ", efficacy_param, ", for all z1, in ",
                  "all considered designs.")
        }
      } else {
        message("  \u2022 You have chosen to include early stopping for ",
                "efficacy, with e1z1 chosen for each z1, in each considered ",
                "design, to control the probability of committing a type-I ",
                "error at the end of stage one to ", efficacy_param)
      }
      if (futility_type == 0) {
        message("  \u2022 You have chosen to prevent early stopping for ",
                "futility. Thus f1z1 = -\u221E, for all z1, in all considered ",
                "designs")
      } else if (futility_type == 1) {
        message("  \u2022 You have chosen to include early stopping for ",
                "futility, with f1z1 = ", futility_param, ", for all z1, in ",
                "all considered designs.")
      } else {
        message("  \u2022 You have chosen to include early stopping for ",
                "futility, with f1z1 chosen for each z1, in each considered ",
                "design, to control the probability of committing a type-II ",
                "error at the end of stage one to ", futility_param)
      }
    } else {
      if (futility) {
        message("  \u2022 You have chosen to allow early stopping for ",
                "futility")
      } else {
        if (type != "single_double") {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "futility. Thus f1 = -\u221E in all considered designs")
        } else {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "futility. Thus fS1 = fT1 = -\u221E in all considered ",
                  "designs")
        }
      }
      if (efficacy) {
        message("  \u2022 You have chosen to allow early stopping for ",
                "efficacy")
      } else {
        if (type != "single_double") {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "efficacy. Thus e1 = \u221E in all considered designs")
        } else {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "efficacy. Thus eS1 = eT1 = \u221E in all considered designs")
        }
      }
    }
    message("\n  The design will be optimised for:")
    message("    w", uc_sub(1), "ESS(", uc("pi"), "O,", uc("pi"), "O) + w",
            uc_sub(2), "ESS(", uc("pi"), "O,", uc("pi"), "O + ", uc("delta"),
            ") + w", uc_sub(3), "max_", uc("pi"), " ESS(", uc("pi"), ",",
            uc("pi"), ") +\n      w", uc_sub(4), "max_{", uc("pi"), "C,",
            uc("pi"), "E} ESS(", uc("pi"), "C,", uc("pi"), "E) + w", uc_sub(5),
            "max N")
    message("  with:")
    message("    w", uc_sub(1)," = ", w[1], ", w", uc_sub(2), " = ", w[2],
            ", w", uc_sub(3), " = ", w[3], ", w", uc_sub(4), " = ", w[4], ", w",
            uc_sub(5), " = ", w[5])
    message("  and ", uc("pi"), "O = ", piO)
  }
}

summary_terminal <- function(des, k) {
  if (des$J == 1) {
    stage  <- "single-stage"
  } else {
    stage  <- "two-stage"
  }
  if (des$type == "barnard") {
    design <- "barnard's exact test"
    dashes <- switch(stage, "single-stage" = 70, "two-stage" = 67)
  } else if (des$type == "binomial") {
    design <- "an exact binomial test"
    dashes <- switch(stage, "single-stage" = 72, "two-stage" = 69)
  } else if (des$type == "fisher") {
    design <- "Fisher's exact test"
    dashes <- switch(stage, "single-stage" = 69, "two-stage" = 66)
  } else if (des$type == "single_double") {
    design <- "single-arm and two-arm testing\n  decisions"
    dashes <- switch(stage, "single-stage" = 80, "two-stage" = 77)
  }
  message("  ", rep("-", dashes))
  message("  Terminal points of a ", stage, " based on ", design)
  message("\n  You have chosen to find the terminal points of a design with:\n")
  if (any(class(des) %in% "ph2rand_des_one_stage")) {
    message("    \u2022 nC1 = ", des$nC)
    message("    \u2022 nE1 = ", des$nE)
    if (des$type %in% c("barnard", "binomial")) {
      message("    \u2022 e1 = ", des$e1)
    } else if (des$type == "fisher") {
      message("    \u2022 e10 = ", des$e1[1], ", ..., e11 = ", des$e1[2],
              ", ..., e1", des$nC + des$nE, " = ", des$e1[des$nC + des$nE + 1])
    } else if (des$type == "single_double") {
      message("    \u2022 eS1 = ", des$eS1)
      message("    \u2022 eT1 = ", des$eT1)
    }
  } else if (any(class(des) %in% "ph2rand_des_two_stage")) {
    message("    \u2022 nC1 = ", des$nC[1])
    message("    \u2022 nC2 = ", des$nC[2])
    message("    \u2022 nE1 = ", des$nE[1])
    message("    \u2022 nE2 = ", des$nE[2])
    if (des$type %in% c("barnard", "binomial")) {
      message("    \u2022 e1 = ", des$e1)
      message("    \u2022 f1 = ", des$f1)
      message("    \u2022 e2 = ", des$e2)
    } else if (des$type == "fisher") {
      message("    \u2022 e10 = ", des$e1[1], ", ..., e11 = ", des$e1[2],
              ", ..., e1", des$nC + des$nE, " = ", des$e1[des$nC + des$nE + 1])
      message("    \u2022 f10 = ", des$f1[1], ", ..., f11 = ", des$f1[2],
              ", ..., f1", des$nC + des$nE, " = ", des$f1[des$nC + des$nE + 1])
      message("    \u2022 e200 = ", des$e2[1, 1], ", ..., ",
              "e2", des$nC[1] + des$nE[1], des$nC[2] + des$nE[2], " = ",
              des$e2[des$nC[1] + des$nE[1], des$nC[2] + des$nE[2]])
    } else if (des$type == "single_double") {
      message("    \u2022 eS1 = ", des$eS1)
      message("    \u2022 eT1 = ", des$eT1)
      message("    \u2022 fS1 = ", des$fS1)
      message("    \u2022 fT1 = ", des$fT1)
      message("    \u2022 eS2 = ", des$eS2)
      message("    \u2022 eT2 = ", des$eT2)
    }
  }
}

summary_pmf      <- function(des, pi, k) {
  if (des$J == 1) {
    stage  <- "single-stage"
  } else {
    stage  <- "two-stage"
  }
  if (des$type == "barnard") {
    design <- "barnard's exact test"
    dashes <- switch(stage, "single-stage" = 58, "two-stage" = 55)
  } else if (des$type == "binomial") {
    design <- "an exact binomial test"
    dashes <- switch(stage, "single-stage" = 60, "two-stage" = 57)
  } else if (des$type == "fisher") {
    design <- "Fisher's exact test"
    dashes <- switch(stage, "single-stage" = 57, "two-stage" = 54)
  } else if (des$type == "single_double") {
    design <- "single-arm and two-arm testing\n  decisions"
    dashes <- switch(stage, "single-stage" = 68, "two-stage" = 65)
  }
  message("  ", rep("-", dashes))
  message("  PMF of a ", stage, " design based on ", design)
  message("  ", rep("-", dashes))
  message("\n  You have chosen to find the PMF of a design with:")
  if (any(class(des) %in% "ph2rand_des_one_stage")) {
    message("      \u2022 nC1 = ", des$nC)
    message("      \u2022 nE1 = ", des$nE)
    if (des$type %in% c("barnard", "binomial")) {
      message("      \u2022 e1 = ", des$e1)
    } else if (des$type == "fisher") {
      message("      \u2022 e10 = ", des$e1[1], ", ..., e11 = ", des$e1[2],
              ", ..., e1", des$nC + des$nE, " = ", des$e1[des$nC + des$nE + 1])
    } else if (des$type == "single_double") {
      message("      \u2022 eS1 = ", des$eS1)
      message("      \u2022 eT1 = ", des$eT1)
    }
  } else if (any(class(des) %in% "ph2rand_des_two_stage")) {
    message("      \u2022 nC1 = ", des$nC[1])
    message("      \u2022 nC2 = ", des$nC[2])
    message("      \u2022 nE1 = ", des$nE[1])
    message("      \u2022 nE2 = ", des$nE[2])
    if (des$type %in% c("barnard", "binomial")) {
      message("      \u2022 e1 = ", des$e1)
      message("      \u2022 f1 = ", des$f1)
      message("      \u2022 e2 = ", des$e2)
    } else if (des$type == "fisher") {
      message("      \u2022 e10 = ", des$e1[1], ", ..., e11 = ", des$e1[2],
              ", ..., e1", des$nC + des$nE, " = ", des$e1[des$nC + des$nE + 1])
      message("      \u2022 f10 = ", des$f1[1], ", ..., f11 = ", des$f1[2],
              ", ..., f1", des$nC + des$nE, " = ", des$f1[des$nC + des$nE + 1])
      message("      \u2022 e200 = ", des$e2[1, 1], ", ..., ",
              "e2", des$nC[1] + des$nE[1], des$nC[2] + des$nE[2], " = ",
              des$e2[des$nC[1] + des$nE[1], des$nC[2] + des$nE[2]])
    } else if (des$type == "single_double") {
      message("      \u2022 eS1 = ", des$eS1)
      message("      \u2022 eT1 = ", des$eT1)
      message("      \u2022 fS1 = ", des$fS1)
      message("      \u2022 fT1 = ", des$fT1)
      message("      \u2022 eS2 = ", des$eS2)
      message("      \u2022 eT2 = ", des$eT2)
    }
  }
  if (nrow(pi) == 1) {
    message("  when ", uc("pi"), " = (", pi[1, 1], ", ", pi[1, 2], ")'.")
  } else if (nrow(pi) == 2) {
    message("  when ", uc("pi"), " \u2208 {(", pi[1, 1], ", ", pi[1, 2],
            ")', (", pi[2, 1], ", ", pi[2, 2], ")'}.")
  } else {
    message("  when ", uc("pi"), " \u2208 {(", pi[1, 1], ", ", pi[1, 2],
            ")', ..., (", pi[nrow(pi), 1], ", ", pi[nrow(pi), 2], ")'}.")
  }
}
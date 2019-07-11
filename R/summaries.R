summary_des      <- function(J, type, alpha, beta, delta, ratio, pi0_null,
                             pi0_alt, n0max, equal, w, pi0_ess, efficacy,
                             futility, efficacy_type, efficacy_param,
                             futility_type, futility_param) {
  if (J == 1) {
    stage    <- "single-stage"
  } else {
    stage    <- "two-stage"
  }
  if (type == "bernard") {
    dashes   <- 60 - 3*(J == 2)
    design   <- "Bernard's exact test"
  } else if (type == "binomial") {
    dashes   <- 62 - 3*(J == 2)
    design   <- "an exact binomial test"
  } else if (type == "fisher") {
    dashes   <- 59 - 3*(J == 2)
    design   <- "Fisher's exact test"
  } else {
    dashed   <- 70 + 7*(J == 2)
    if (J == 1) {
      design <- "single-arm and two-arm testing\n  decisions"
    } else {
      design <- "single-arm and two-arm testing decisions"
    }
  }
  message("  ", rep("-", dashes))
  message("  Design of a ", stage, " trial based on ", design)
  message("  ", rep("-", dashes))
  message("\n  Hypothesis test")
  message("  ---------------")
  message("  You have chosen to test the following hypothesis:")
  if (length(pi0_null) == 1) {
    message("      H", uc_sub(0), ": ", uc("pi"), uc_sub(0), " = ", uc("pi"),
            uc_sub(1), " \u2208 ", uc("Pi"), uc_sub(0), " = ", pi0_null)
  } else {
    message("      H", uc_sub(0), ": ", uc("pi"), uc_sub(0), " = ", uc("pi"),
            uc_sub(1), "\u2208 ", uc("Pi"), uc_sub(0), " = [", pi0_null[1], ",",
            pi0_null[2], "]")
  }
  message("  with the following type-I error constraint:")
  if (length(pi0_null) == 1) {
    message("      P(", pi0_null, ",", pi0_null, ") ", uc("le"), " ",
            uc("alpha"), " = ", alpha)
  } else {
    message("      max_{", uc("pi"), " \u2208 ", uc("Pi"), uc_sub(0), "} P(",
            uc("pi"), ",", uc("pi"), ") ", uc("le"), " ", uc("alpha"), " = ",
            alpha, ", ", uc("Pi"), uc_sub(0), " = [", pi0_null[1], ",",
            pi0_null[2], "]")
  }
  message("  and the following type-II error constraint:")
  if (length(pi0_alt) == 1) {
    message("      P(", pi_alt, ",", pi_alt + delta, ") ", uc("ge"), " 1 - ",
            uc("beta"), " = ", 1 - beta)
  } else {
    message("      max_{", uc("pi"), " \u2208 ", uc("Pi"), uc_sub(0), "} P(",
            uc("pi"), ",", uc("pi"), " + ", uc("delta"), ") ", uc("ge"),
            " 1 - ", uc("beta"), " = ", 1 - beta, ", ", uc("Pi"), uc_sub(1),
            " = [", pi_alt[1], ",", pi_alt[2], "], ", uc("delta"), " = ",
            delta)
  }
  message("\n  Restrictions")
  message("  ------------")
  message("  \u2022 You have chosen to limit the allowed maximal sample ",
          "size in the control arm, n", uc_sub(0), ", to ", n0max)
  if (J == 1) {
    message("  \u2022 The sample size in the experimental arm, n", uc_sub(1),
            ", will be set to rn", uc_sub(0), ", with r = ", ratio)
  } else {
    if (equal) {
      message("  \u2022 You have chosen to restrict the sample sizes in the",
              " control arm in each stage, n", uc_sub(0), uc_sub(1),
              " and n", uc_sub(0), uc_sub(2), ", such that n", uc_sub(0),
              uc_sub(1), " = ", "n", uc_sub(0), uc_sub(2))
    } else {
      message("  \u2022 You have chosen to allow the sample sizes in the",
              " control arm in each stage, n", uc_sub(0), uc_sub(1),
              " and n", uc_sub(0), uc_sub(2), ", to take different values")
    }
    message("  \u2022 The sample sizes in the experimental arm in each ",
            "stage, n", uc_sub(1), uc_sub(1), " and n", uc_sub(1), uc_sub(2),
            "will be set to rn", uc_sub(0), uc_sub(1), " and rn", uc_sub(0),
            uc_sub(2), " respectively, with r = ", ratio)
    if (type == "fisher") {
      if (efficacy_type == 0) {
        message("  \u2022 You have chosen to prevent early stopping for ",
                "efficacy. Thus e_z", uc_sub(1), " = \u221E, for all z",
                uc_sub(1), " in all considered designs.")
      } else if (efficacy_type == 1) {
        if (efficacy_param == -0.5) {
          message("  \u2022 You have chosen to include early stopping for ",
                  "efficacy, with e_z", uc_sub(1), " = [0.5(n", uc_sub(0),
                  uc_sub(1), " + n", uc_sub(1), uc_sub(1), ")", uc("delta"),
                  "]_* + 1, for all z", uc_sub(1), ", in all considered ",
                  "designs")
        } else {
          message("  \u2022 You have chosen to include early stopping for ",
                  "efficacy, with e_z", uc_sub(1), " chosen for each z",
                  uc_sub(1), ", in each considered design, to control the ",
                  "probability of committing a type-I error at the end of ",
                  "stage one to ", efficacy_param)
        }
      }
      if (futility_type == 0) {
        message("  \u2022 You have chosen to prevent early stopping for ",
                "futility. Thus f_z", uc_sub(1), " = -\u221E, for all z",
                uc_sub(1), ", in all considered designs")
      } else if (futility_type == 1) {
        message("  \u2022 You have chosen to include early stopping for ",
                "futility, with f_z", uc_sub(1), " = ",
                futility_param, ", for all z", uc_sub(1),
                ", in all considered designs.")
      } else {
        message("  \u2022 You have chosen to include early stopping for ",
                "futility, with f_z", uc_sub(1), " chosen for each z",
                uc_sub(1), ", in each considered design, to control the ",
                "probability of committing a type-II error at the end of stage",
                " one to ", futility_param)
      }
    } else {
      if (futility) {
        message("  \u2022 You have chosen to allow early stopping for ",
                "futility")
      } else {
        if (type != "single_double") {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "futility. Thus f", uc_sub(1), " = -\u221E in all considered",
                  " designs.")
        } else {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "futility. Thus f_S", uc_sub(1), " = ", "f_T", uc_sub(1),
                  "-\u221E in all considered designs.")
        }
      }
      if (efficacy) {
        message("  \u2022 You have chosen to allow early stopping for ",
                "efficacy")
      } else {
        if (type != "single_double") {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "efficacy. Thus e", uc_sub(1), " = \u221E in all considered",
                  " designs.")
        } else {
          message("  \u2022 You have chosen to prevent early stopping for ",
                  "efficacy. Thus e_S", uc_sub(1), " = ", "e_T", uc_sub(1),
                  "\u221E in all considered designs.")
        }
      }
    }
    message("  The design will be optimised for:")
    message("      w", uc_sub(1), "ESS(", uc("pi"), "_ESS", uc("pi"), "_ESS) +",
            " w", uc_sub(2), "ESS(", uc("pi"), "_ESS", uc("pi"), "_ESS + ",
            uc("delta"), ") + w", uc_sub(3), "max_", uc("pi"), " ESS(",
            uc("pi"), ",", uc("pi"), ") + w", uc_sub(4), "max_{", uc("pi"),
            uc_sub(0), ",", uc("pi"), uc_sub(1), "} ESS(", uc("pi"), uc_sub(0),
            ",", uc("pi"), uc_sub(1), ") + w", uc_sub(5), "max N")
    message("  with:")
    message("      w", uc_sub(1), " = ", w[1], ", w", uc_sub(2), " = ", w[2],
            ", w", uc_sub(3), " = ", w[3], ", w", uc_sub(4), " = ", w[4], ", w",
            uc_sub(5), " = ", w[5])
  }
}

summary_terminal <- function(des, k) {
  if (any(class(des) %in% "ph2rand_des_one_stage")) {
    stage  <- "single-stage"
  } else if (any(class(des) %in% "ph2rand_des_two_stage")) {
    stage  <- "two-stage"
  }
  if (des$type == "bernard") {
    design <- "Bernard's exact test"
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
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    if (des$type %in% c("bernard", "binomial")) {
      message("      \u2022 e = ", des$e)
    } else if (des$type == "fisher") {
      message("      \u2022 e\u2080 = ", des$e[1], ", ..., e\u2081 = ",
              des$e[2], ", ..., e", uc_sub(des$n0 + des$n1), " = ",
              des$e[des$n0 + des$n1 + 1])
    } else if (des$type == "single_double") {
      message("      \u2022 e_S = ", des$eS)
      message("      \u2022 e_T = ", des$eT)
    }
  } else if (any(class(des) %in% "ph2rand_des_two_stage")) {
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    if (des$type %in% c("bernard", "binomial")) {
      message("      \u2022 e\u2081 = ", des$e1)
      message("      \u2022 f\u2081 = ", des$f1)
      message("      \u2022 e\u2082 = ", des$e2)
    } else if (des$type == "fisher") {
      message("      \u2022 e\u2081\u2080 = ", des$e1[1], ", ..., e\u2081",
              "\u2081 = ", des$e1[2], ", ..., e\u2081", uc_sub(des$n0 + des$n1),
              " = ", des$e1[des$n0 + des$n1 + 1])
      message("      \u2022 f\u2081\u2080 = ", des$f1[1], ", ..., f\u2081",
              "\u2081 = ", des$f1[2], ", ..., f\u2081", uc_sub(des$n0 + des$n1),
              " = ", des$f1[des$n0 + des$n1 + 1])
      message("      \u2022 e\u2082\u2080\u2080 = ", des$e2[1, 1], ", ..., ",
              "e\u2082", uc_sub(des$n0[1] + des$n1[1]),
              uc_sub(des$n0[2] + des$n1[2]), " = ",
              des$e2[des$n0[1] + des$n1[1], des$n0[2] + des$n1[2]])
    } else if (des$type == "single_double") {
      message("      \u2022 e_S\u2081 = ", des$eS1)
      message("      \u2022 e_T\u2081 = ", des$eT1)
      message("      \u2022 f_S\u2081 = ", des$fS1)
      message("      \u2022 f_T\u2081 = ", des$fT1)
      message("      \u2022 e_S\u2082 = ", des$eS2)
      message("      \u2022 e_T\u2082 = ", des$eT2)
    }
  }
}

summary_pmf      <- function(des, pi, k) {
  if (any(class(des) %in% "ph2rand_des_one_stage")) {
    stage  <- "single-stage"
  } else if (any(class(des) %in% "ph2rand_des_two_stage")) {
    stage  <- "two-stage"
  }
  if (des$type == "bernard") {
    design <- "Bernard's exact test"
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
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    if (des$type %in% c("bernard", "binomial")) {
      message("      \u2022 e = ", des$e)
    } else if (des$type == "fisher") {
      message("      \u2022 e\u2080 = ", des$e[1], ", ..., e\u2081 = ",
              des$e[2], ", ..., e", uc_sub(des$n0 + des$n1), " = ",
              des$e[des$n0 + des$n1 + 1])
    } else if (des$type == "single_double") {
      message("      \u2022 e_S = ", des$eS)
      message("      \u2022 e_T = ", des$eT)
    }
  } else if (any(class(des) %in% "ph2rand_des_two_stage")) {
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    if (des$type %in% c("bernard", "binomial")) {
      message("      \u2022 e\u2081 = ", des$e1)
      message("      \u2022 f\u2081 = ", des$f1)
      message("      \u2022 e\u2082 = ", des$e2)
    } else if (des$type == "fisher") {
      message("      \u2022 e\u2081\u2080 = ", des$e1[1], ", ..., e\u2081",
              "\u2081 = ", des$e1[2], ", ..., e\u2081", uc_sub(des$n0 + des$n1),
              " = ", des$e1[des$n0 + des$n1 + 1])
      message("      \u2022 f\u2081\u2080 = ", des$f1[1], ", ..., f\u2081",
              "\u2081 = ", des$f1[2], ", ..., f\u2081", uc_sub(des$n0 + des$n1),
              " = ", des$f1[des$n0 + des$n1 + 1])
      message("      \u2022 e\u2082\u2080\u2080 = ", des$e2[1, 1], ", ..., ",
              "e\u2082", uc_sub(des$n0[1] + des$n1[1]),
              uc_sub(des$n0[2] + des$n1[2]), " = ",
              des$e2[des$n0[1] + des$n1[1], des$n0[2] + des$n1[2]])
    } else if (des$type == "single_double") {
      message("      \u2022 e_S\u2081 = ", des$eS1)
      message("      \u2022 e_T\u2081 = ", des$eT1)
      message("      \u2022 f_S\u2081 = ", des$fS1)
      message("      \u2022 f_T\u2081 = ", des$fT1)
      message("      \u2022 e_S\u2082 = ", des$eS2)
      message("      \u2022 e_T\u2082 = ", des$eT2)
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
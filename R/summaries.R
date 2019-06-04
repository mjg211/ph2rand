summary_des_one_stage <- function(type, alpha, beta, delta, ratio, point_null,
                                  pi_null, point_alt, pi_alt, n0max) {
  if (type %in% c("bernard", "binomial")) {
    repeats <- 57
  } else if (type == "fisher") {
    repeats <- 56
  } else {
    repeats <- 77
  }
  message("   ", rep("-", repeats))
  if (all(J == 1, type == "bernard")) {
    message("   Design of a one-stage trial based on Bernard's exact test")
  } else if (all(J == 1, type == "binomial")) {
    message("   Design of a one-stage trial based on exact binomial tests")
  } else if (all(J == 1, type == "fisher")) {
    message("   Design of a one-stage trial based on Fisher's exact test")
  } else if (all(J == 1, type == "single_double")) {
    message("   Design of a one-stage trial based on single-arm and ",
            "two-arm testing decisions")
  } else if (all(J == 2, type == "bernard")) {
    message("   Design of a two-stage trial based on Bernard's exact test")
  } else if (all(J == 2, type == "binomial")) {
    message("   Design of a two-stage trial based on exact binomial tests")
  } else if (all(J == 2, type == "fisher")) {
    message("   Design of a two-stage trial based on Fisher's exact test")
  } else if (all(J == 2, type == "single_double")) {
    message("   Design of a two-stage trial based on single-arm and ",
            "two-arm testing decisions")
  }
  message("   ", rep("-", repeats))
  Sys.sleep(1)
  message("\n   Hypothesis test")
  message("   ---------------")
  message("   You have chosen to test the following hypothesis")
  Sys.sleep(1)
  if (point_null) {
    message("      H\u2080: \u03c0\u2080 = \u03c0\u2081 \u2208 \u03A0\u2080 ",
            "= ", pi_null, ",")
  } else {
    message("      H\u2080: \u03c0\u2080 = \u03c0\u2081 \u2208 \u03A0\u2080 ",
            "= [", pi_null[1], ",", pi_null[2], "],")
  }
  Sys.sleep(1)
  message("   with the following type-I error constraint")
  Sys.sleep(1)
  if (point_null) {
    message("      P(", pi_null, ",", pi_null, ") \u2264 \u03b1 = ", alpha,
            ",")
  } else {
    message("      max_{\u03c0 \u2208 \u03A0\u2080} P(\u03c0,\u03c0) \u2264 ",
            "\u03b1 = ", alpha, ",")
  }
  Sys.sleep(1)
  message("   and the following type-II error constraint")
  Sys.sleep(1)
  if (point_alt) {
    message("      P(", pi_alt, ",", pi_alt + delta, ") \u2265 1 - \u03b2 = ",
            1 - beta, ".")
  } else {
    message("      max_{\u03c0 \u2208 \u03A0\u2081} P(\u03c0,\u03c0 + \u03b4",
            ") \u2265 1 - \u03b2 = ", 1 - beta, ", \u03A0\u2081 = [",
            pi_alt[1], ",", pi_alt[2], "], \u03b4 = ", delta, ".")
  }
  Sys.sleep(1)
  message("\n   Restrictions")
  message("   ------------")
  message("   \u2022 You have chosen to limit the allowed maximal sample ",
          "size in the control arm, n\u2080, to ", n0max, ".")
  Sys.sleep(1)
  if (J == 1) {
    if (ratio == 1) {
      message("   \u2022 The sample size in the experimental arm, n\u2081, ",
              "will be set to n\u2080.")
    } else {
      message("   \u2022 The sample size in the experimental arm, n\u2081, ",
              "will be set to ", ratio, "n\u2080.")
    }
  } else {
    if (two_stage$equal) {
      message("   \u2022 You have chosen to restrict the sample sizes in the",
              " control arm in each stage, n\u2080\u2081 and n\u2080\u2082, ",
              "such that n\u2080\u2081 = n\u2080\u2082.")
    } else {
      message("   \u2022 You have chosen to allow the sample sizes in the ",
              "control arm in each stage, n\u2080\u2081 and n\u2080\u2082, ",
              "to take different values.")
    }
    Sys.sleep(1)
    if (ratio == 1) {
      message("   \u2022 The sample sizes in the experimental arm in each ",
              "stage, n\u2081\u2081 and n\u2081\u2082, will be set to ",
              "n\u2081\u2081 = n\u2080\u2081", " and n\u2081\u2082 = ",
              "n\u2080\u2082.")
    } else {
      message("   \u2022 The sample sizes in the experimental arm in each ",
              "stage, n\u2081\u2081 and n\u2081\u2082, will be set to ",
              "n\u2081\u2081 = ", ratio, "n\u2080\u2081", " and ",
              "n\u2081\u2082 = ", ratio, "n\u2080\u2082.")
    }
    Sys.sleep(1)
    if (type == "fisher") {
      if (two_stage$efficacy_type == 0) {
        message("   \u2022 You have chosen to prevent early stopping for ",
                "efficacy. Thus e_z\u2081 = \u221E, for all z\u2081, in all ",
                "considered designs.")
      } else if (two_stage$efficacy_type == 1) {
        if (two_stage$efficacy_param == -0.5) {
          message("   \u2022 You have chosen to include early stopping for ",
                  "efficacy, with e_z\u2081 = [0.5(n\u2080\u2081 + ",
                  "n\u2081\u2081)\u03B4]_* + 1, for all z\u2081, in all ",
                  "considered designs.")
        } else {
          message("   \u2022 You have chosen to include early stopping for ",
                  "efficacy, with e_z\u2081 chosen for each z\u2081, in each ",
                  "considered design, to control the probability of ",
                  "committing a type-I error at the end of stage one to ",
                  two_stage$efficacy_param)
        }
      } else {
        message("   \u2022 You have chosen to prevent early stopping for ",
                "efficacy. Thus e_z\u2081 = \u221E, for all z\u2081, in all ",
                "considered designs.")
      }
      Sys.sleep(1)
      if (two_stage$futility_type == 0) {
        message("   \u2022 You have chosen to prevent early stopping for ",
                "futility. Thus f_z\u2081 = = -\u221E, for all z\u2081, in ",
                "all considered designs.")
      } else if (two_stage$futility_type == 1) {
        message("   \u2022 You have chosen to include early stopping for ",
                "futility, with f_z\u2081 = ", two_stage$futility_param,
                " for all z\u2081, in all considered designs.")
      } else {
        message("   \u2022 You have chosen to include early stopping for ",
                "futility, with f_z\u2081 chosen for each z\u2081, in each ",
                "considered design, to control the probability of committing",
                " a type-II error at the end of stage one to ",
                two_stage$futility_param)
      }
    } else {
      
    }
  }
  Sys.sleep(1)
  message("\n   Beginning the required calculations...")
}

summary_terminal <- function(des, k) {
  if (all(any(class(des) %in% "ph2rand_bernard"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a single-stage design based on an exact ",
            "Bernard test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e = ", des$e)
  } else if (all(any(class(des) %in% "ph2rand_binomial"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a single-stage design based on an exact ",
            "binomial test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e = ", des$e)
  } else if (all(any(class(des) %in% "ph2rand_fisher"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a single-stage design based on Fisher's ",
            "exact test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e\u2080 = ", des$e[1], ", ..., e\u2081 = ",
            des$e[2], ", ..., e", subnum(des$n0 + des$n1), " = ",
            des$e[des$n0 + des$n1 + 1])
  } else if (all(any(class(des) %in% "ph2rand_single_double"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a single-stage design based on single-arm",
            " and two-arm testing decisions")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e_S = ", des$eS)
    message("      \u2022 e_T = ", des$eT)
  } else if (all(any(class(des) %in% "ph2rand_bernard"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a two-stage design based on exact ",
            "Bernard tests")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e\u2081 = ", des$e1)
    message("      \u2022 f\u2081 = ", des$f1)
    message("      \u2022 e\u2082 = ", des$e2)
  } else if (all(any(class(des) %in% "ph2rand_binomial"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a two-stage design based on exact ",
            "binomial tests")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e\u2081 = ", des$e1)
    message("      \u2022 f\u2081 = ", des$f1)
    message("      \u2022 e\u2082 = ", des$e2)
  } else if (all(any(class(des) %in% "ph2rand_fisher"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a two-stage design based on Fisher's ",
            "exact test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e\u2081\u2080 = ", des$e1[1], ", ..., e\u2081",
            "\u2081 = ", des$e1[2], ", ..., e\u2081", subnum(des$n0 + des$n1),
            " = ", des$e1[des$n0 + des$n1 + 1])
    message("      \u2022 f\u2081\u2080 = ", des$f1[1], ", ..., f\u2081",
            "\u2081 = ", des$f1[2], ", ..., f\u2081", subnum(des$n0 + des$n1),
            " = ", des$f1[des$n0 + des$n1 + 1])
    message("      \u2022 e\u2082\u2080\u2080 = ", des$e2[1, 1], ", ..., ",
            "e\u2082", subnum(des$n0[1] + des$n1[1]),
            subnum(des$n0[2] + des$n1[2]), " = ",
            des$e2[des$n0[1] + des$n1[1], des$n0[2] + des$n1[2]])
  } else if (all(any(class(des) %in% "ph2rand_single_double"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Terminal points of a two-stage design based on single-arm",
            " and two-arm testing decisions")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the terminal points of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e_S\u2081 = ", des$eS1)
    message("      \u2022 e_T\u2081 = ", des$eT1)
    message("      \u2022 f_S\u2081 = ", des$fS1)
    message("      \u2022 f_T\u2081 = ", des$fT1)
    message("      \u2022 e_S\u2082 = ", des$eS2)
    message("      \u2022 e_T\u2082 = ", des$eT2)
  }
  Sys.sleep(1)
  message("   Beginning the required calculations...")
}

summary_pmf <- function(des, pi, k) {
  if (all(any(class(des) %in% "ph2rand_bernard"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a single-stage design based on an exact ",
            "Bernard test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e = ", des$e)
  } else if (all(any(class(des) %in% "ph2rand_binomial"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a single-stage design based on an exact ",
            "binomial test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e = ", des$e)
  } else if (all(any(class(des) %in% "ph2rand_fisher"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a single-stage design based on Fisher's ",
            "exact test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e\u2080 = ", des$e[1], ", ..., e\u2081 = ",
            des$e[2], ", ..., e", subnum(des$n0 + des$n1), " = ",
            des$e[des$n0 + des$n1 + 1])
  } else if (all(any(class(des) %in% "ph2rand_single_double"), des$J == 1)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a single-stage design based on single-arm",
            " and two-arm testing decisions")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080 = ", des$n0)
    message("      \u2022 n\u2081 = ", des$n1)
    message("      \u2022 e_S = ", des$eS)
    message("      \u2022 e_T = ", des$eT)
  } else if (all(any(class(des) %in% "ph2rand_bernard"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a two-stage design based on exact ",
            "Bernard tests")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e\u2081 = ", des$e1)
    message("      \u2022 f\u2081 = ", des$f1)
    message("      \u2022 e\u2082 = ", des$e2)
  } else if (all(any(class(des) %in% "ph2rand_binomial"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a two-stage design based on exact ",
            "binomial tests")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e\u2081 = ", des$e1)
    message("      \u2022 f\u2081 = ", des$f1)
    message("      \u2022 e\u2082 = ", des$e2)
  } else if (all(any(class(des) %in% "ph2rand_fisher"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a two-stage design based on Fisher's ",
            "exact test")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e\u2081\u2080 = ", des$e1[1], ", ..., e\u2081",
            "\u2081 = ", des$e1[2], ", ..., e\u2081", subnum(des$n0 + des$n1),
            " = ", des$e1[des$n0 + des$n1 + 1])
    message("      \u2022 f\u2081\u2080 = ", des$f1[1], ", ..., f\u2081",
            "\u2081 = ", des$f1[2], ", ..., f\u2081", subnum(des$n0 + des$n1),
            " = ", des$f1[des$n0 + des$n1 + 1])
    message("      \u2022 e\u2082\u2080\u2080 = ", des$e2[1, 1], ", ..., ",
            "e\u2082", subnum(des$n0[1] + des$n1[1]),
            subnum(des$n0[2] + des$n1[2]), " = ",
            des$e2[des$n0[1] + des$n1[1], des$n0[2] + des$n1[2]])
  } else if (all(any(class(des) %in% "ph2rand_single_double"), des$J == 2)) {
    message("   ", rep("-", 10))
    message("   Probability mass function of a two-stage design based on single-arm",
            " and two-arm testing decisions")
    message("   ", rep("-", 10))
    Sys.sleep(1)
    message("\n   You have chosen to find the probability mass function of a design ",
            "with")
    Sys.sleep(1)
    message("      \u2022 n\u2080\u2081 = ", des$n0[1])
    message("      \u2022 n\u2080\u2082 = ", des$n0[2])
    message("      \u2022 n\u2081\u2081 = ", des$n1[1])
    message("      \u2022 n\u2081\u2082 = ", des$n1[2])
    message("      \u2022 e_S\u2081 = ", des$eS1)
    message("      \u2022 e_T\u2081 = ", des$eT1)
    message("      \u2022 f_S\u2081 = ", des$fS1)
    message("      \u2022 f_T\u2081 = ", des$fT1)
    message("      \u2022 e_S\u2082 = ", des$eS2)
    message("      \u2022 e_T\u2082 = ", des$eT2)
  }
  Sys.sleep(1)
  if (rows(pi) == 1) {
    message("   when pi = (", pi[1, 1], ", ", pi[1, 2], ")'.")
  } else if (rows(pi) == 2) {
    message("   when pi \u2208 {(", pi[1, 1], ", ", pi[1, 2], ")', (",
            pi[2, 1], ", ", pi[2, 2], ")'}.")
  } else {
    message("   when pi \u2208 {(", pi[1, 1], ", ", pi[1, 2], ")', ..., (",
            pi[rows(pi), 1], ", ", pi[rows(pi), 2], ")'}.")
  }
  Sys.sleep(1)
}
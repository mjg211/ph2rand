#' @export
des_two_stage <- function(type = "binomial", alpha = 0.1, beta = 0.2,
                          delta = 0.2, ratio = 1, point_null = T, pi_null = 0.1,
                          point_alt = T, pi_alt = 0.1, n0max = 100L, equal = T,
                          w = c(1, 0, 0, 0, 0), pi_ess = pi_null[1],
                          efficacy = F, futility = T, efficacy_type = 0L,
                          efficacy_param = NULL, futility_type = 1L,
                          futility_param = 0L, summary = F) {

  ##### Check inputs ###########################################################

  check_belong(type, "type", c("bernard", "binomial", "fisher",
                               "single_double"), 1)
  check_real_range_strict(alpha, "alpha", c(0, 1),   1)
  check_real_range_strict(beta,  "beta",  c(0, 1),   1)
  check_real_range_strict(delta, "delta", c(0, 1),   1)
  check_real_range_strict(ratio, "ratio", c(0, Inf), 1)
  check_logical(point_null, "point_null")
  check_pi_null(pi_null, point_null)
  check_logical(point_alt, "point_alt")
  check_pi_alt(pi_alt, point_alt, delta)
  n0max <- check_integer_range(n0max, "n0max", c(1, Inf), 1)
  check_logical(equal, "equal")
  check_w(w)
  check_real_range(pi_ess, "pi_ess", c(0, 1 - delta), 1)
  if (type == "fisher") {
    check_default(efficacy, "efficacy", F)
    check_default(futility, "futility", T)
    check_integer_range(efficacy_type, "efficacy_type", c(-1, 3), 1)
    check_integer_range(futility_type, "futility_type", c(-1, 3), 1)
    efficacy_param <- check_fisher_params(efficacy_type, efficacy_param,
                                          futility_type, futility_param)
  } else {
    check_logical(efficacy, "efficacy")
    check_logical(futility, "futility")
    check_default(efficacy_type, "efficacy_type", 0)
    check_default(efficacy_param, "efficacy_param", NULL)
    check_default(futility_type, "futility_type", 1)
    check_default(futility_param, "futility_param", 0)
  }
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    #summary_des_two_stage(type, alpha, beta, delta, ratio, point_null, pi_null,
    #                      point_alt, pi_alt, n0max, equal, w, pi_ess, efficacy,
    #                      futility, efficacy_type, efficacy_param,
    #                      futility_type, futility_param)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("  Identifying feasible designs", uc("two_elip"))
  }
  output <- switch(type,
                   bernard       =
                     bernard_des_two_stage(alpha, beta, delta, ratio,
                                           point_null, pi_null, point_alt,
                                           pi_alt, n0max, equal, w, pi_ess,
                                           efficacy, futility, summary),
                   binomial      =
                     binomial_des_two_stage(alpha, beta, delta, ratio,
                                            point_null, pi_null, point_alt,
                                            pi_alt, n0max, equal, w, pi_ess,
                                            efficacy, futility, summary),
                   fisher        =
                     fisher_des_two_stage(alpha, beta, delta, ratio,
                                          point_null, pi_null, point_alt,
                                          pi_alt, n0max, equal, w, pi_ess,
                                          efficacy_type, efficacy_param,
                                          futility_type, futility_param,
                                          summary),
                   single_double =
                     single_double_des_two_stage(alpha, beta, delta, ratio,
                                                 point_null, pi_null,
                                                 point_alt, pi_alt, n0max,
                                                 equal, w, pi_ess, efficacy,
                                                 futility, summary))
  if (summary) {
    message(uc("two_elip"), "outputting.")
  }

  ##### Output #################################################################

  return(output)

}
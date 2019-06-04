#' @export
des_one_stage <- function(type = "binomial", alpha = 0.1, beta = 0.2,
                          delta = 0.2, ratio = 1, point_null = T, pi_null = 0.1,
                          point_alt = T, pi_alt = 0.1, n0max = 100L,
                          summary = F) {

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
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    #summary_des_one_stage(type, alpha, beta, delta, ratio, point_null, pi_null,
    #                      point_alt, pi_alt, n0max)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("  Identifying feasible designs", uc("two_elip"))
  }
  output <- switch(type,
                   bernard       =
                     bernard_des_one_stage(alpha, beta, delta, ratio,
                                           point_null, pi_null, point_alt,
                                           pi_alt, n0max, summary),
                   binomial      =
                     binomial_des_one_stage(alpha, beta, delta, ratio,
                                            point_null, pi_null, point_alt,
                                            pi_alt, n0max, summary),
                   fisher        =
                     fisher_des_one_stage(alpha, beta, delta, ratio,
                                          point_null, pi_null, point_alt,
                                          pi_alt, n0max, summary),
                   single_double =
                     single_double_des_one_stage(alpha, beta, delta, ratio,
                                                 point_null, pi_null,
                                                 point_alt, pi_alt, n0max,
                                                 summary))
  if (summary) {
    message(uc("two_elip"), "outputting.")
  }
  
  ##### Output #################################################################

  return(output)

}
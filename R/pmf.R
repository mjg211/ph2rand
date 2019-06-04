#' @export
pmf <- function(des, pi, k, summary = F) {
  
  ##### Check inputs ###########################################################
  
  check_des(des, "any", 1:2)
  pi <- check_pi(pi)
  k  <- check_k(k, des)
  check_logical(summary, "summary")
  
  ##### Print summary ##########################################################
  
  if (summary) {
    #summary_pmf(des, pi, k)
  }
  
  ##### Perform main computations ##############################################
  
  if (summary) {
    message("  Identifying PMFs", uc("two_elip"))
  }
  if ("ph2rand_des_one_stage" %in% class(des)) {
    pmf <- switch(type,
                  bernard       =
                    bernard_pmf_one_stage(pi, des$n0, des$n1, des$e),
                  binomial      =
                    binomial_pmf_one_stage(pi, des$n0, des$n1, des$e),
                  fisher        =
                    fisher_pmf_one_stage(pi, des$n0, des$n1, des$e),
                  single_double =
                    single_double_pmf_one_stage(pi, des$n0, des$n1, des$eS,
                                                   des$eT))
  } else {
    pmf <- switch(type,
                  bernard       =
                    bernard_pmf_two_stage(pi, des$n0, des$n1, des$e1, des$f1,
                                          des$e2, k),
                  binomial      =
                    binomial_pmf_two_stage(pi, des$n0, des$n1, des$e1, des$f1,
                                           des$e2, k),
                  fisher        =
                    fisher_pmf_two_stage(pi, des$n0, des$n1, des$e1, des$f1,
                                         des$e2, k),
                  single_double =
                    single_double_pmf_two_stage(pi, des$n0, des$n1, des$eS1,
                                                des$eT1, des$fS1, des$fT1,
                                                des$eS2, des$eT2, k))
  }
  if (summary) {
    message(uc("two_elip"), "outputting.")
  }
  
  ##### Output results #########################################################
  
  output        <- list(des     = des,
                        k       = k,
                        pmf     = pmf,
                        summary = summary)
  class(output) <- c(class(output), "ph2rand_pmf")
  output
}
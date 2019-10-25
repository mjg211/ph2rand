#' @export
pmf <- function(des, pi, k, summary = F) {

  ##### Check inputs ###########################################################

  check_des(des, "any")
  pi <- check_pi(pi, des)
  k  <- check_k(k, des)
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_pmf(des, pi, k)
  }

  ##### Perform main computations ##############################################

  if (summary) {
    message("  Identifying PMFs", uc("two_elip"))
  }
  if (des$J == 1) {
    pmf <- switch(des$type,
                  barnard       =
                    barnard_pmf_one_stage(pi, des$nC, des$nE, des$e1),
                  binomial      =
                    binomial_pmf_one_stage(pi, des$nC, des$nE, des$e1),
                  fisher        =
                    fisher_pmf_one_stage(pi, des$nC, des$nE, des$e1),
                  single_double =
                    single_double_pmf_one_stage(pi, des$nC, des$nE, des$eS1,
                                                   des$eT1))
  } else {
    pmf <- switch(des$type,
                  barnard       =
                    barnard_pmf_two_stage(pi, des$nC, des$nE, des$e1, des$f1,
                                          des$e2, k),
                  binomial      =
                    binomial_pmf_two_stage(pi, des$nC, des$nE, des$e1, des$f1,
                                           des$e2, k),
                  fisher        =
                    fisher_pmf_two_stage(pi, des$nC, des$nE, des$e1, des$f1,
                                         des$e2, k),
                  single_double =
                    single_double_pmf_two_stage(pi, des$nC, des$nE, des$eS1,
                                                des$eT1, des$fS1, des$fT1,
                                                des$eS2, des$eT2, k))
  }
  if (summary) {
    message(uc("two_elip"), "outputting")
  }

  ##### Output results #########################################################

  output        <- list(des     = des,
                        k       = k,
                        pmf     = pmf,
                        summary = summary)
  class(output) <- c("ph2rand_pmf", class(output))
  output

}
#' Probability mass functions of a two-arm randomised clinical trial design for
#' a binary primary outcome variable
#'
#' \code{pmf} determines probability mass functions of a design returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}, under given
#' response rate scenarios (see \code{pi}).
#' 
#' @param des An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' @param pi A \code{\link{matrix}} with two columns, giving the response rate
#' scenarios to consider. The first column should correspond to the control arm
#' and the second column to the experimental arm. Defaults internally to the
#' values in \code{des$pi}.
#' @param k A \code{\link{numeric}} \code{\link{vector}} indicating which stages
#' to consider in determining the probability mass functions. That is, it will
#' condition the calculations on the trial ending in the stages given in
#' \code{k}. Defaults internally to all stages of the given design.
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{F}.
#' @return An object of class \code{"ph2rand_pmf"}, containing each of the
#' input parameters along with a tibble in the slot \code{$pmf}, which gives the
#' determined probability mass functions.
#' @examples
#' # The default two-stage design
#' des <- des_two_stage()
#' # Its probability mass function under the uninteresting and interesting
#' # scenarios
#' pmf <- pmf(des)
#' # The same probability mass functions, conditioning on the trial ending in
#' # stage 2
#' pmf <- pmf(des, k = 2)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{plot.ph2rand_terminal}}.
#' @export
pmf <- function(des, pi, k, summary = F) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(des, "any")
  pi <- check_pi(pi, des)
  k  <- check_k(k, des)
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_pmf(des, pi, k)
  }

  ##### Perform main computations ##############################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying PMFs...")
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
    message("..outputting")
  }

  ##### Output results #########################################################

  output        <- list(des     = des,
                        k       = k,
                        pmf     = pmf,
                        summary = summary)
  class(output) <- "ph2rand_pmf"
  output

}
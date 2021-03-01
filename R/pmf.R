#' Probability mass functions of a randomised clinical trial design that assumes
#' a Bernoulli distributed primary outcome variable
#'
#' \code{pmf} determines probability mass functions of a design returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}, under given
#' response rate scenarios (see \code{pi}).
#' 
#' @param des An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}. Defaults to
#' \code{ph2rand::des_one_stage()}.
#' @param pi A \code{\link{numeric}} \code{\link{vector}} with two elements, or
#' a \code{\link{numeric}} \code{\link{matrix}} or \code{\link{data.frame}} with
#' two columns, giving the response rate scenarios to consider. The first
#' element/column should correspond to the control arm and the second
#' element/column to the experimental arm. Defaults to \code{des$opchar[, 1:2]}.
#' @param k A \code{\link{numeric}} \code{\link{vector}} indicating which stages
#' to consider in determining the probability mass functions. That is, it will
#' condition the calculations on the trial ending in the stages given in
#' \code{k}. Defaults to \code{1:des$J} (i.e., to all stages of the given
#' design).
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{FALSE}.
#' @return A \code{\link{list}} with additional class \code{"ph2rand_pmf"},
#' containing each of the input parameters along with a tibble in the slot
#' \code{$pmf}, which gives the determined probability mass functions.
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
pmf <- function(des = ph2rand::des_one_stage(), pi = des$opchar[, 1:2],
                k = 1:des$J, summary = FALSE) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(des, "any")
  pi <- check_pi(pi, des)
  check_k(k, des)
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_pmf(des, pi, k)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying probability mass functions..")
  }
  if (des$J == 1) {
    pmf <- switch(des$type,
                  barnard  = barnard_pmf_one_stage(pi, des$nC, des$nE,
                                                   des$boundaries$e1),
                  binomial = binomial_pmf_one_stage(pi, des$nC, des$nE,
                                                    des$boundaries$e1),
                  fisher   = fisher_pmf_one_stage(pi, des$nC, des$nE,
                                                  des$boundaries$e1),
                  sat      = sat_pmf_one_stage(pi, des$nC, des$nE,
                                               des$boundaries$eS1,
                                               des$boundaries$eT1))
  } else {
    pmf <- switch(des$type,
                  barnard  = barnard_pmf_two_stage(pi, des$nC, des$nE,
                                                   des$boundaries$e1,
                                                   des$boundaries$f1,
                                                   des$boundaries$e2, k),
                  binomial = binomial_pmf_two_stage(pi, des$nC, des$nE,
                                                    des$boundaries$e1,
                                                    des$boundaries$f1,
                                                    des$boundaries$e2, k),
                  fisher   = fisher_pmf_two_stage(pi, des$nC, des$nE,
                                                  des$boundaries$e1,
                                                  des$boundaries$f1,
                                                  des$boundaries$e2, k),
                  sat      = sat_pmf_two_stage(pi, des$nC, des$nE,
                                               des$boundaries$eS1,
                                               des$boundaries$eT1,
                                               des$boundaries$fS1,
                                               des$boundaries$fT1,
                                               des$boundaries$eS2,
                                               des$boundaries$eT2, k))
  }
  if (summary) {
    message("..outputting")
  }

  ##### Output results #########################################################

  output        <- list(des     = des,
                        k       = k,
                        pi      = pi,
                        pmf     = pmf,
                        summary = summary)
  class(output) <- c("ph2rand_pmf", class(output))
  output

}
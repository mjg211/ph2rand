#' Determine operating characteristics of a randomised clinical trial design
#' that assumes a Bernoulli distributed primary outcome variable
#'
#' \code{opchar} determines the operating characteristics (analytically) of a
#' design returned by \code{\link{des_one_stage}} or
#' \code{\link{des_two_stage}}, under given response rate scenarios (see
#' \code{pi}).
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
#' to consider in determining the operating characteristics. That is, it will
#' condition the calculations on the trial ending in the stages given in
#' \code{k}. Defaults to \code{1:des$J} (i.e., to all stages of the given
#' design).
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{FALSE}.
#' @return A \code{\link{list}} with additional class \code{"ph2rand_opchar"},
#' containing each of the input parameters along with a tibble in the slot
#' \code{$opchar}, which gives the determined operating characteristics.
#' @examples
#' # The default two-stage design
#' des    <- des_two_stage()
#' # Its operating characteristics under the uninteresting and interesting
#' # scenarios
#' opchar <- opchar(des)
#' # The same operating characteristics, conditioning on the trial ending in
#' # stage 2
#' opchar <- opchar(des, k = 2)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}}.
#' @export
opchar <- function(des = ph2rand::des_one_stage(), pi = des$opchar[, 1:2],
                   k = 1:des$J, summary = FALSE) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(des, "any")
  pi <- check_pi(pi, des)
  check_k(k, des)
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_opchar(des, pi, k)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying operating characteristics..")
  }
  if (des$J == 1) {
    opchar <- switch(des$type,
                     barnard  = barnard_opchar_one_stage(pi, des$nC, des$nE,
                                                         des$boundaries$e1),
                     binomial = binomial_opchar_one_stage(pi, des$nC, des$nE,
                                                          des$boundaries$e1),
                     fisher   = fisher_opchar_one_stage(pi, des$nC, des$nE,
                                                        des$boundaries$e1),
                     sat      = sat_opchar_one_stage(pi, des$nC, des$nE,
                                                     des$boundaries$eS1,
                                                     des$boundaries$eT1))
  } else {
    opchar <- switch(des$type,
                     barnard  = barnard_opchar_two_stage(pi, des$nC, des$nE,
                                                         des$boundaries$e1,
                                                         des$boundaries$f1,
                                                         des$boundaries$e2, k),
                     binomial = binomial_opchar_two_stage(pi, des$nC, des$nE,
                                                          des$boundaries$e1,
                                                          des$boundaries$f1,
                                                          des$boundaries$e2, k),
                     fisher   = fisher_opchar_two_stage(pi, des$nC, des$nE,
                                                        des$boundaries$e1,
                                                        des$boundaries$f1,
                                                        des$boundaries$e2, k),
                     sat      = sat_opchar_two_stage(pi, des$nC, des$nE,
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
                        opchar  = opchar,
                        pi      = pi,
                        summary = summary)
  class(output) <- c("ph2rand_opchar", class(output))
  output

}
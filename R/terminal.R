#' Terminal points of a randomised clinical trial design that assumes a
#' Bernoulli distributed primary outcome variable
#'
#' \code{terminal} determines the 'terminal' points of a design returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' 
#' @param des An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}. Defaults to
#' \code{ph2rand::des_one_stage()}.
#' @param k A \code{\link{numeric}} \code{\link{vector}} indicating which stages
#' to consider when determining the terminal points. Defaults to \code{1:des$J}
#' (i.e., to all stages of the given design).
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{FALSE}.
#' @return A \code{\link{list}} with additional class \code{"ph2rand_terminal"},
#' containing each of the input parameters along with a tibble in the slot
#' \code{$terminal}, which gives the determined terminal points.
#' @examples
#' # The default two-stage design
#' des     <- des_two_stage()
#' # Its terminal points across stages 1 and 2
#' term_12 <- terminal(des)
#' # Its terminal points from stage 2 only
#' term_2  <- terminal(des, 2)
#' # A plot of these points
#' plot(term_2)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{plot.ph2rand_terminal}}.
#' @export
terminal <- function(des = ph2rand::des_one_stage(), k = 1:des$J,
                     summary = FALSE) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(des, "any")
  check_k(k, des)
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_terminal(des, k)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying terminal points..")
  }
  if (des$J == 1) {
    terminal <- switch(des$type,
                       barnard  = barnard_terminal_one_stage(des$nC, des$nE,
                                                             des$boundaries$e1),
                       binomial =
                         binomial_terminal_one_stage(des$nC, des$nE,
                                                     des$boundaries$e1),
                       fisher   = fisher_terminal_one_stage(des$nC, des$nE,
                                                            des$boundaries$e1),
                       sat      = sat_terminal_one_stage(des$nC, des$nE,
                                                         des$boundaries$e1))
  } else {
    terminal <- switch(des$type,
                       barnard  = barnard_terminal_two_stage(des$nC, des$nE,
                                                             des$boundaries$e1,
                                                             des$boundaries$f1,
                                                             des$boundaries$e2,
                                                             k),
                       binomial = binomial_terminal_two_stage(des$nC, des$nE,
                                                              des$boundaries$e1,
                                                              des$boundaries$f1,
                                                              des$boundaries$e2,
                                                              k),
                       fisher   = fisher_terminal_two_stage(des$nC, des$nE,
                                                            des$boundaries$e1,
                                                            des$boundaries$f1,
                                                            des$boundaries$e2,
                                                            k),
                       sat       = sat_terminal_two_stage(des$nC, des$nE,
                                                          des$boundaries$eS1,
                                                          des$boundaries$eT1,
                                                          des$boundaries$fS1,
                                                          des$boundaries$fT1,
                                                          des$boundaries$eS2,
                                                          des$boundaries$eT2,
                                                          k))
  }
  if (summary) {
    message("..outputting")
  }

  ##### Output results #########################################################

  output        <- list(des      = des,
                        k        = k,
                        summary  = summary,
                        terminal = terminal)
  class(output) <- c("ph2rand_terminal", class(output))
  output

}
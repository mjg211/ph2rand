#' Terminal points of a two-arm randomised clinical trial design for a binary
#' primary outcome variable
#'
#' \code{terminal} determines terminal points of a design returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' 
#' @param des An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' @param k A \code{\link{numeric}} \code{\link{vector}} indicating which stages
#' to consider in determining the terminal points. Defaults internally to all
#' stages of the given design.
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{F}.
#' @return An object of class \code{"ph2rand_terminal"}, containing each of the
#' input parameters along with a tibble in the slot \code{$terminal}, which
#' gives the determined terminal points.
#' @examples
#' # The default two-stage design
#' des  <- des_two_stage()
#' # Its terminal points across stages 1 and 2
#' term <- terminal(des)
#' # Its terminal points from stage 2 only
#' term <- terminal(des, 2)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{plot.ph2rand_terminal}}.
#' @export
terminal <- function(des, k, summary = F) {
  
  ##### Check inputs ###########################################################
  
  check_ph2rand_des(des, "any")
  k <- check_k(k, des)
  check_logical(summary, "summary")
  
  ##### Print summary ##########################################################
  
  if (summary) {
    summary_terminal(des, k)
  }

  ##### Perform main computations ##############################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying terminal points...")
  }
  if (des$J == 1) {
    terminal <- switch(des$type,
                       barnard       =
                         barnard_terminal_one_stage(des$nC, des$nE, des$e1),
                       binomial      =
                         binomial_terminal_one_stage(des$nC, des$nE, des$e1),
                       fisher        =
                         fisher_terminal_one_stage(des$nC, des$nE, des$e1),
                       single_double =
                         single_double_terminal_one_stage(des$nC, des$nE,
                                                          des$e1))
  } else {
    terminal <- switch(des$type,
                       barnard       =
                         barnard_terminal_two_stage(des$nC, des$nE, des$e1,
                                                    des$f1, des$e2, k),
                       binomial      =
                         binomial_terminal_two_stage(des$nC, des$nE, des$e1,
                                                     des$f1, des$e2, k),
                       fisher        =
                         fisher_terminal_two_stage(des$nC, des$nE, des$e1,
                                                   des$f1, des$e2, k),
                       single_double =
                         single_double_terminal_two_stage(des$nC, des$nE,
                                                          des$eS1, des$eT1,
                                                          des$fS1, des$fT1,
                                                          des$eS2, des$eT2, k))
  }
  if (summary) {
    message("..outputting")
  }

  ##### Output results #########################################################

  output        <- list(des      = des,
                        k        = k,
                        summary  = summary,
                        terminal = terminal)
  class(output) <- "ph2rand_terminal"
  output

}
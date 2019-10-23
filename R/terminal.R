#' @export
terminal <- function(des, k, summary = F) {
  
  ##### Check inputs ###########################################################
  
  check_des(des, "any")
  k <- check_k(k, des)
  check_logical(summary, "summary")
  
  ##### Print summary ##########################################################
  
  if (summary) {
    summary_terminal(des, k)
  }

  ##### Perform main computations ##############################################

  if (summary) {
    message("  Identifying terminal points", uc("two_elip"))
  }
  if (des$J == 1) {
    terminal <- switch(des$type,
                       bernard       =
                         bernard_terminal_one_stage(des$nC, des$nE, des$e1),
                       binomial      =
                         binomial_terminal_one_stage(des$nC, des$nE, des$e1),
                       fisher        =
                         fisher_terminal_one_stage(des$nC, des$nE, des$e1),
                       single_double =
                         single_double_terminal_one_stage(des$nC, des$nE,
                                                          des$e1))
  } else {
    terminal <- switch(des$type,
                       bernard       =
                         bernard_terminal_two_stage(des$nC, des$nE, des$e1,
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
    message(uc("two_elip"), "outputting")
  }

  ##### Output results #########################################################

  output        <- list(des      = des,
                        k        = k,
                        summary  = summary,
                        terminal = terminal)
  class(output) <- c("ph2rand_terminal", class(output))
  output

}
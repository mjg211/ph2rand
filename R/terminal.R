#' @export
terminal <- function(des, k, summary = F) {
  
  ##### Check inputs ###########################################################
  
  check_des(des, "any", 1:2)
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
  if ("ph2rand_des_one_stage" %in% class(des)) {
    terminal <- switch(des$type,
                       bernard       =
                         bernard_terminal_one_stage(des$n0, des$n1, des$e),
                       binomial      =
                         binomial_terminal_one_stage(des$n0, des$n1, des$e),
                       fisher        =
                         fisher_terminal_one_stage(des$n0, des$n1, des$e),
                       single_double =
                         single_double_terminal_one_stage(des$n0, des$n1,
                                                          des$e))
  } else {
    terminal <- switch(des$type,
                       bernard       =
                         bernard_terminal_two_stage(des$n0, des$n1, des$e1,
                                                    des$f1, des$e2, k),
                       binomial      =
                         binomial_terminal_two_stage(des$n0, des$n1, des$e1,
                                                     des$f1, des$e2, k),
                       fisher        =
                         fisher_terminal_two_stage(des$n0, des$n1, des$e1,
                                                   des$f1, des$e2, k),
                       single_double =
                         single_double_terminal_two_stage(des$n0, des$n1,
                                                          des$eS1, des$eT1,
                                                          des$fS1, des$fT1,
                                                          des$eS2, des$eT2, k))
  }
  if (summary) {
    message(uc("two_elip"), "outputting.")
  }

  ##### Output results #########################################################

  output        <- list(des      = des,
                        k        = k,
                        summary  = summary,
                        terminal = terminal)
  class(output) <- c(class(output), "ph2rand_terminal")
  output

}
#' Summarise a two-arm randomised clinical trial design for a binary primary
#' outcome variable
#'
#' \code{summary.ph2rand_des} prints a summary of a design returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' 
#' @param x An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' @examples
#' # The default two-stage design
#' des   <- des_two_stage()
#' # Print a summary
#' summary(des)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{plot.ph2rand_des}}.
#' @export
summary.ph2rand_des <- function(x, ...) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(x, "any", "x")

  ##### Print summary ##########################################################

  summary_ph2rand_des(x)

}
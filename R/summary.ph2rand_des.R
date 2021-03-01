#' Summarise a randomised clinical trial design that assumes a Bernoulli
#' distributed primary outcome variable
#'
#' \code{summary.ph2rand_des} prints a summary of a design returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' 
#' @param object An object of class \code{ph2rand_des}, as returned by
#' \code{\link{des_one_stage}} or \code{\link{des_two_stage}}.
#' @param ... Not currently used.
#' @return Currently not used.
#' @examples
#' # The default two-stage design
#' des <- des_two_stage()
#' # Print a summary
#' summary(des)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}},
#' \code{\link{plot.ph2rand_des}}.
#' @method summary ph2rand_des
#' @export
summary.ph2rand_des <- function(object, ...) {

  ##### Check inputs ###########################################################

  check_ph2rand_des(object, "any", "object")

  ##### Print summary ##########################################################

  summary_ph2rand_des(object)

}
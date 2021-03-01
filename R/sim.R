#' Estimate operating characteristics of a randomised clinical trial design that
#' assumes a Bernoulli distributed primary outcome variable
#'
#' \code{sim} estimated the operating characteristics (via simulation) of a
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
#' @param replicates A \code{\link{numeric}} indicating the number of replicate
#' simulations to use for each value of
#' \ifelse{html}{\out{<b><i>&pi;</i></b>}}{\eqn{\bold{\pi}}}. Defaults to
#' \code{1e4}.
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{FALSE}.
#' @return A \code{\link{list}} with additional class \code{"ph2rand_sim"},
#' containing each of the input parameters along with a tibble in the slot
#' \code{$sim}, which gives the estimated operating characteristics.
#' @examples
#' # The default two-stage design
#' des <- des_two_stage()
#' # Its operating characteristics under the uninteresting and interesting
#' # scenarios
#' sim <- sim(des)
#' # The same operating characteristics, conditioning on the trial ending in
#' # stage 2
#' sim <- sim(des, k = 2)
#' @seealso \code{\link{des_one_stage}}, \code{\link{des_two_stage}}.
#' @export
sim <- function(des = ph2rand::des_one_stage(), pi = des$opchar[, 1:2],
                k = 1:des$J, replicates = 1e4, summary = FALSE) {
  
  ##### Check inputs ###########################################################
  
  check_ph2rand_des(des, "any")
  pi         <- check_pi(pi, des)
  check_k(k, des)
  replicates <- check_integer_range(replicates, "replicates", c(1, Inf), 1)
  check_logical(summary, "summary")
  
  ##### Print summary ##########################################################
  
  if (summary) {
    summary_sim(des, pi, k, replicates)
  }
  
  ##### Main computations ######################################################
  
  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Estimating operating characteristics..")
  }
  nrow_pi          <- nrow(pi)
  total_replicates <- nrow_pi*replicates
  sim              <- matrix(0, nrow_pi, ifelse(des$J == 1, 3, 13))
  for (i in 1:nrow_pi) {
    sim[i, ]       <- sim_internal(pi[i, ], (i - 1)*replicates, des, k,
                                   replicates, summary, total_replicates)
  }
  if (des$J == 1) {
    sim            <- tibble::tibble(piC     = sim[, 1],
                                     piE     = sim[, 2],
                                     `P(pi)` = sim[, 3])
  } else {
    colnames(sim)  <- c("piC", "piE", "P(pi)", "ESS(pi)", "SDSS(pi)", "MSS(pi)",
                        paste0(rep(c("E", "F", "S"), each = 2), rep(1:2, 3),
                               "(pi)"), "max N")
    sim            <- tibble::as_tibble(sim)
    sim$`max N`    <- as.integer(sim$`max N`)
  }
  if (summary) {
    message("..completed the required simulations..")
  }
  
  ##### Output results #########################################################
  
  if (summary) {
    message("..outputting")
  }
  output        <- list(des        = des,
                        k          = k,
                        pi         = pi,
                        replicates = replicates,
                        sim        = sim,
                        summary    = summary)
  class(output) <- c("ph2rand_sim", class(output))
  output
  
}
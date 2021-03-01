#' Design a one-stage two-arm randomised clinical trial assuming a Bernoulli
#' distributed primary outcome variable
#'
#' \code{des_one_stage} determines one-stage two-arm randomised clinical trial
#' designs, assuming the primary outcome variable is Bernoulli distributed. It
#' supports a flexible framework for specifying which scenarios to control the
#' type-I and type-II error-rates for, and allows for design determination
#' assuming a variety of test statistics. In all instances, \code{des_one_stage}
#' computes the optimal required sample size in each arm, the associated optimal
#' stopping boundaries, and returns information on key operating
#' characteristics.
#'
#' @param type A \code{\link{character}} string indicating the chosen design
#' framework/test statistic to assume. Must be one of \code{"barnard"},
#' \code{"binomial"}, \code{"fisher"}, or \code{"sat"}. Defaults to
#' \code{"binomial"}.
#' @param alpha A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>&alpha;</i>}}{\eqn{\alpha}}, the significance level
#' (i.e., the type-I error-rate). Defaults to \code{0.1}.
#' @param beta A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>&beta;</i>}}{\eqn{\beta}}, used in the definition of
#' the desired power (i.e., the type-II error-rate). Defaults to \code{0.2}.
#' @param delta A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>&delta;</i>}}{\eqn{\delta}}, the treatment effect
#' assumed in the power calculation. Defaults to \code{0.2}.
#' @param ratio A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>r</i>}}{\eqn{r}}, the allocation ratio to
#' the experimental arm, relative to the control arm. Defaults to \code{1}.
#' @param Pi0 A \code{\link{numeric}} \code{\link{vector}} indicating the
#' chosen value for \ifelse{html}{\out{<i>&Pi;</i><sub>0</sub>}}{\eqn{\Pi_0}},
#' the control arm response rates to control the type-I error-rate to level
#' \ifelse{html}{\out{<i>&alpha;</i>}}{\eqn{\alpha}} for. Must either be of
#' \code{\link{length}} one, indicating a single point, or of
#' \code{\link{length}} two. In this case, the elements indicate the
#' range of possible response rates to allow for. Defaults to \code{0.1}.
#' @param Pi1 A \code{\link{numeric}} \code{\link{vector}} indicating the
#' chosen value for \ifelse{html}{\out{<i>&Pi;</i><sub>1</sub>}}{\eqn{\Pi_1}},
#' the control arm response rates to allow for in the power calculations. Must
#' either be of \code{\link{length}} one, indicating a single point, or of
#' \code{\link{length}} two. In this case, the elements indicate the range of
#' possible response rates to allow for. Defaults to \code{Pi0[1]}.
#' @param nCmax A \code{\link{numeric}} indicating the maximum value of the
#' sample size in the control arm to consider in the search procedure. Defaults
#' to \code{50L}.
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{FALSE}.
#' @return A \code{\link{list}} with additional class \code{"ph2rand_des"},
#' containing each of the input parameters along with several additional
#' variables, including
#' \itemize{
#' \item A \code{\link{list}} in the slot \code{$boundaries} giving the
#' rejection boundary/boundaries of the optimal design. The names of these
#' elements depends on the value of \code{type}.
#' \item A \code{\link{tibble}} in the slot \code{$feasible} summarising the
#' operating characteristics of the feasible designs.
#' \item A \code{\link{numeric}} in the slot \code{$nC} giving the sample size
#' in the control arm for the optimal design.
#' \item A \code{\link{numeric}} in the slot \code{$nE} giving the sample size
#' in the experimental arm for the optimal design.
#' \item A \code{\link{tibble}} in the slot \code{$opchar} summarising the
#' operating characteristics of the optimal design.
#' }
#' @examples
#' # The design for the default parameters
#' des       <- des_one_stage()
#' # Controlling the type-I/II error-rates over a range of possible response
#' # rates
#' des_range <- des_one_stage(Pi0 = c(0, 1),
#'                            Pi1 = c(0, 0.8))
#' @seealso \code{\link{des_two_stage}}, \code{\link{opchar}},
#' \code{\link{pmf}}, \code{\link{sim}}, \code{\link{terminal}},
#' \code{\link{plot.ph2rand_des}}, \code{\link{summary.ph2rand_des}}.
#' @export
des_one_stage <- function(type = "binomial", alpha = 0.1, beta = 0.2,
                          delta = 0.2, ratio = 1, Pi0 = 0.1, Pi1 = Pi0[1],
                          nCmax = 50L, summary = FALSE) {

  ##### Check inputs ###########################################################

  check_belong(type, "type", c("barnard", "binomial", "fisher", "sat"), 1)
  check_real_range_strict(alpha, "alpha", c(0, 1),   1)
  check_real_range_strict(beta,  "beta",  c(0, 1),   1)
  check_real_range_strict(delta, "delta", c(0, 1),   1)
  check_real_range_strict(ratio, "ratio", c(0, Inf), 1)
  check_Pi0(Pi0)
  check_Pi1(Pi1, delta)
  nCmax <- check_integer_range(nCmax, "nCmax", c(1, Inf), 1)
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_des(1, type, alpha, beta, delta, ratio, Pi0, Pi1, nCmax)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying feasible designs..")
  }
  output <- switch(type,
                   barnard  = barnard_des_one_stage(alpha, beta, delta, ratio,
                                                    Pi0, Pi1, nCmax, summary),
                   binomial = binomial_des_one_stage(alpha, beta, delta, ratio,
                                                     Pi0, Pi1, nCmax, summary),
                   fisher   = fisher_des_one_stage(alpha, beta, delta, ratio,
                                                   Pi0, Pi1, nCmax, summary),
                   sat      = sat_des_one_stage(alpha, beta, delta, ratio, Pi0,
                                                Pi1, nCmax, summary))
  if (summary) {
    message("..outputting")
  }

  ##### Output #################################################################

  output

}
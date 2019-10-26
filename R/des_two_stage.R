#' Design a two-stage two-arm randomised clinical trial for a binary primary
#' outcome variable
#'
#' \code{des_two_stage} determines one-stage two-arm randomised clinical trial
#' designs, assuming the primary outcome variable is binary. More specifically,
#' under the assumption that it is Bernoulli distributed. It supports a flexible
#' framework for specifying which scenarios to control the type-I and type-II
#' error-rates for, and allows for design determination assuming a variety of
#' test statistics. In all instances, \code{des_two_stage} computes the relevant
#' required sample size in each arm, and returns information on key operating
#' characteristics.
#'
#' @param type A \code{\link{character}} string indicating the chosen design
#' framework/test statistic to assume. Must be one of \code{"barnard"},
#' \code{"binomial"}, \code{"fisher"}, or \code{"single_double"}. Defaults to
#' \code{"binomial"}.
#' @param alpha A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>&alpha;</i>}}{\eqn{\alpha}}, the significance level.
#' Defaults to \code{0.1}.
#' @param beta A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>&beta;</i>}}{\eqn{\beta}}, used in the definition of
#' the desired power. Defaults to \code{0.2}.
#' @param delta A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<i>&delta;</i>}}{\eqn{\delta}}, the desired treatment
#' effect. Defaults to \code{0.2}.
#' @param ratio A \code{\link{numeric}} indicating the chosen value for
#' \ifelse{html}{\out{<b><i>r</i></b>}}{\eqn{\bold{r}}}, the allocation ratio to
#' the experimental arm, relative to the control arm. Defaults to \code{1}.
#' @param Pi0 A \code{\link{numeric}} \code{\link{vector}} indicating the
#' chosen values of the control arm response rate to allow for in the null
#' hypothesis. Must either be of \code{\link{length}} one, indicating a point
#' null, or of \code{\link{length}} two. In this case, the elements indicate the
#' range of possible response rates to allow for. Defaults to \code{0.1}.
#' @param Pi1 A \code{\link{numeric}} \code{\link{vector}} indicating the
#' chosen values of the control arm response rate to allow for in the
#' alternative hypothesis. Must either be of \code{\link{length}} one,
#' indicating a point null, or of \code{\link{length}} two. Then, the elements
#' indicate the range of possible response rates to allow for. Defaults to
#' \code{Pi0[1]}.
#' @param nCmax A \code{\link{numeric}} indicating the maximum value of the
#' sample size in the control arm (across both stages) to consider. Defaults to
#' \code{50L}.
#' @param equal A \code{\link{logical}} variable indicating whether the sample
#' size for a given arm should be equal in the two stages. Defaults to \code{T}.
#' @param w A \code{\link{numeric}} \code{\link{vector}} indicating the weights
#' to use in the optimality criteria. Must be of \code{\link{length}} five, with
#' all elements greater than or equal to zero, and at least one of the first
#' four elements strictly positive. Defaults to \code{c(1, 0, 0, 0, 0)}.
#' @param piO A \code{\link{numeric}} indicating the value of the control
#' arm response rate to assume in the optimality criteria. Defaults to
#' \code{Pi0[1]}.
#' @param efficacy Only used if \code{type} is one of \code{"barnard"},
#' \code{"binomial"}, or \code{"single_double"}. Then, it is a
#' \code{\link{logical}} variable indicating whether to include early stopping
#' for efficacy in the design. Defaults to \code{F}.
#' @param futility Only used if \code{type} is one of \code{"barnard"},
#' \code{"binomial"}, or \code{"single_double"}. Then, it is a
#' \code{\link{logical}} variable indicating whether to include early stopping
#' for futility in the design. Defaults to \code{F}.
#' @param efficacy_type Only used if \code{type} is \code{"fisher"}. Then, it is
#' a \code{\link{numeric}} indicating whether, and which type of, early stopping
#' for efficacy to include in the design. See the vignette for details. Defaults
#' to \code{0L}.
#' @param efficacy_param Only used if \code{type} is \code{"fisher"} and
#' \code{efficacy_type} is not equal to \code{0L}. Then, it is a
#' \code{\link{numeric}} that influences the precise way in which an efficacy
#' boundary is specified. See the vignette for details. Defaults to \code{NULL}.
#' @param futility_type Only used if \code{type} is \code{"fisher"}. Then, it is
#' a \code{\link{numeric}} indicating whether, and which type of, early stopping
#' for futility to include in the design. See the vignette for details. Defaults
#' to \code{1L}.
#' @param futility_param Only used if \code{type} is \code{"fisher"} and
#' \code{futility_type} is not equal to \code{0L}. Then, it is a
#' \code{\link{numeric}} that influences the precise way in which a futility
#' boundary is specified. See the vignette for details. Defaults to \code{1L}.
#' @param summary A \code{\link{logical}} variable indicating whether a summary
#' of the function's progress should be printed to the console. Defaults to
#' \code{F}.
#' @return An object of class \code{"ph2rand_des"}, containing each of the input
#' parameters along with several additional variables, including:
#' \itemize{
#' \item A selection of elements that prescribe the rejection boundaries of the
#' optimal design. The names of these elements depends on the value of
#' \code{type}.
#' \item A \code{\link{tibble}} in the slot \code{$feasible} summarising the
#' operating characteristics of the feasible designs.
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot \code{$nC}
#' giving the sample sizes in the control arm in each stage for the optimal
#' design.
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot \code{$nE}
#' giving the sample sizes in the experimental arm in each stage for the optimal
#' design.
#' \item A \code{\link{tibble}} in the slot \code{$opchar} summarising the
#' operating characteristics of the optimal design.
#' }
#' @examples
#' # The design for the default parameters
#' des       <- des_two_stage()
#' # Controlling the type-I/II error-rates over a range of possible response
#' # rates
#' des_range <- des_two_stage(Pi0 = c(0, 1),
#'                            Pi1 = c(0, 0.8))
#' @seealso \code{\link{des_one_stage}}, \code{\link{opchar}},
#' \code{\link{pmf}}, \code{\link{terminal}}, \code{\link{plot.ph2rand_des}},
#' \code{\link{summary.ph2rand_des}}.
#' @export
des_two_stage <- function(type = "binomial", alpha = 0.1, beta = 0.2,
                          delta = 0.2, ratio = 1, Pi0 = 0.1, Pi1 = Pi0[1],
                          nCmax = 50L, equal = T, w = c(1, 0, 0, 0, 0),
                          piO = Pi0[1], efficacy = F, futility = T,
                          efficacy_type = 0L, efficacy_param = NULL,
                          futility_type = 1L, futility_param = 0L,
                          summary = F) {

  ##### Check inputs ###########################################################

  check_belong(type, "type", c("barnard", "binomial", "fisher",
                               "single_double"), 1)
  check_real_range_strict(alpha, "alpha", c(0, 1),   1)
  check_real_range_strict(beta,  "beta",  c(0, 1),   1)
  check_real_range_strict(delta, "delta", c(0, 1),   1)
  check_real_range_strict(ratio, "ratio", c(0, Inf), 1)
  check_Pi0(Pi0)
  check_Pi1(Pi1, delta)
  nCmax <- check_integer_range(nCmax, "nCmax", c(1, Inf), 1)
  check_logical(equal, "equal")
  check_w(w)
  check_real_range(piO, "piO", c(0, 1 - delta), 1)
  if (type == "fisher") {
    check_default(efficacy, "efficacy", F)
    check_default(futility, "futility", T)
    check_integer_range(efficacy_type, "efficacy_type", c(-1, 3), 1)
    check_integer_range(futility_type, "futility_type", c(-1, 3), 1)
    efficacy_param <- check_fisher_params(efficacy_type, efficacy_param,
                                          futility_type, futility_param)
  } else {
    check_logical(efficacy, "efficacy")
    check_logical(futility, "futility")
    check_default(efficacy_type,  "efficacy_type",  0)
    check_default(efficacy_param, "efficacy_param", NULL)
    check_default(futility_type,  "futility_type",  1)
    check_default(futility_param, "futility_param", 0)
  }
  check_logical(summary, "summary")

  ##### Print summary ##########################################################

  if (summary) {
    summary_des(2, type, alpha, beta, delta, ratio, Pi0, Pi1, nCmax, equal, w,
                piO, efficacy, futility, efficacy_type, efficacy_param,
                futility_type, futility_param)
  }

  ##### Main computations ######################################################

  if (summary) {
    message("\n  ------------")
    message("  Computations")
    message("  ------------")
    message("  Identifying feasible designs...")
  }
  output <- switch(type,
                   barnard       =
                     barnard_des_two_stage(alpha, beta, delta, ratio, Pi0, Pi1,
                                           nCmax, equal, w, piO, efficacy,
                                           futility, summary),
                   binomial      =
                     binomial_des_two_stage(alpha, beta, delta, ratio, Pi0, Pi1,
                                            nCmax, equal, w, piO, efficacy,
                                            futility, summary),
                   fisher        =
                     fisher_des_two_stage(alpha, beta, delta, ratio, Pi0, Pi1,
                                          nCmax, equal, w, piO, efficacy_type,
                                          efficacy_param, futility_type,
                                          futility_param, summary),
                   single_double =
                     single_double_des_two_stage(alpha, beta, delta, ratio,
                                                 Pi0, Pi1, nCmax, equal, w, piO,
                                                 efficacy, futility, summary))
  if (summary) {
    message("..outputting")
  }

  ##### Output #################################################################

  output

}
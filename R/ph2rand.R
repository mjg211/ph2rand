#' ph2rand: Design of Randomized Comparative Phase II Oncology Trials with a
#' Bernoulli Primary Outcome
#'
#' \strong{ph2rand} provides a suite of functions to assist with the design of
#' randomized comparative phase II oncology trials with a Bernoulli primary
#' outcome variable. Specifically, support is provided to (a) perform a sample
#' size calculation when using one of several published designs, (b) evaluate
#' the operating characteristics of a given design, and (c) produce informative
#' plots.
#'
#' @section Getting started:
#'
#' You can install the latest development version of \strong{ph2rand} from
#' \href{https://github.com/}{Github} with:
#'
#' \code{devtools::install_github("mjg211/ph2rand")}
#'
#' An introductory example of how to make use of the package's core
#' functionality can be found \href{https://github.com/mjg211/multiarm}{here}.
#' For further help, please email \email{michael.grayling@@newcastle.ac.uk}.
#'
#' @section Details:
#' Currently, the following functions are available (exported)
#'
#' \itemize{
#' \item \code{\link{des_one_stage}}: Determine a one-stage randomized
#' comparative phase II oncology trial design, assuming a Bernoulli primary
#' outcome variable.
#' \item \code{\link{des_two_stage}}: Determine a two-stage randomized
#' comparative phase II oncology trial design, assuming a Bernoulli primary
#' outcome variable.
#' \item \code{\link{opchar}}: Evaluate the operating characteristics of a
#' randomized comparative phase II oncology trial design that assumes a
#' Bernoulli primary outcome variable.
#' \item \code{\link{plot.ph2rand_des}}: Plot the operating characteristics of a
#' randomized comparative phase II oncology trial design that assumes a
#' Bernoulli primary outcome variable.
#' \item \code{\link{plot.ph2rand_pmf}}: Plot the probability mass function of a
#' randomized comparative phase II oncology trial design that assumes a
#' Bernoulli primary outcome variable.
#' \item \code{\link{plot.ph2rand_terminal}}: Plot the terminal points of a
#' randomized comparative phase II oncology trial design that assumes a
#' Bernoulli primary outcome variable.
#' \item \code{\link{pmf}}: Find the probability mass function of a randomized
#' comparative phase II oncology trial design that assumes a Bernoulli primary
#' outcome variable.
#' \item \code{\link{summary.ph2rand_des}}: Display a summary of a randomized
#' comparative phase II oncology trial design that assumes a Bernoulli primary
#' outcome variable.
#' \item \code{\link{terminal}}: Find the terminal points of a randomized
#' comparative phase II oncology trial design that assumes a Bernoulli primary
#' outcome variable.
#' }
#'
#' @docType package
#' @name ph2rand
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(".data")
}

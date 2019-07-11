#' ph2rand: Design of randomized comparative phase II oncology trials
#'
#' \strong{multiarm} provides a suite of functions to assist with the design of
#' randomized comparative phase II oncology trials. Specifically, support is
#' provided to perform a sample size calculation for each of the most popular
#' randomized comparative phase II oncology trial designs, for both point and
#' composite null hypotheses.
#'
#' @section Getting started:
#'
#' You can install the latest development version of ph2rand from
#' \href{https://github.com/}{Github} with:
#'
#' \code{devtools::install_github("mjg211/ph2rand")}
#'
#' An introductory example of how to make use of the package's core
#' functionality can be found \href{https://github.com/mjg211/multiarm}{here}.
#' For further help, please contact Michael Grayling at
#' \email{michael.grayling@@newcastle.ac.uk}.
#'
#' @section Details:
#' Currently, the following functions are provided:
#'
#' \itemize{
#' \item \code{\link{des_one_stage}}: Design a single-stage randomized
#' comparative phase II oncology trial.
#' \item \code{\link{des_two_stage}}: Design a two-stage randomized comparative
#' phase II oncology trial.
#' \item \code{\link{pmf}}: Find the probability mass function of a randomized
#' comparative phase II oncology trial.
#' \item \code{\link{terminal}}: Find the terminal points of a randomized
#' comparative phase II oncology trial.
#' }
#'
#' @docType package
#' @name ph2rand
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(".data")
}

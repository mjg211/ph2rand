#' ph2rand: Randomized Phase II Oncology Trials with Bernoulli Outcomes
#'
#' \strong{ph2rand} provides functions to assist with the design of randomized
#' comparative phase II oncology trials that assume their primary outcome
#' variable is Bernoulli distributed. Specifically, support is provided to (a)
#' perform a sample size calculation when using one of several published
#' designs, (b) evaluate the operating characteristics of a given design (both
#' analytically and via simulation), and (c) produce informative plots.
#'
#' @section Getting started:
#'
#' You can install the latest development version of \strong{ph2rand} from
#' \href{https://github.com/}{Github} with:
#'
#' \code{devtools::install_github("mjg211/ph2rand")}
#'
#' An introductory example of how to make use of the package's core
#' functionality can be found \href{https://github.com/mjg211/ph2rand}{here}.
#' For further help, please see the package vignettes or email
#' \email{michael.grayling@@newcastle.ac.uk}.
#'
#' @section Details:
#' 
#' Currently, the following functions are available (exported)
#'
#' \itemize{
#' \item \code{\link{des_one_stage}}: Determine a one-stage randomized
#' comparative phase II oncology trial design, assuming a Bernoulli primary
#' outcome variable.
#' \item \code{\link{des_two_stage}}: Determine a two-stage randomized
#' comparative phase II oncology trial design, assuming a Bernoulli primary
#' outcome variable.
#' \item \code{\link{opchar}}: Evaluate the operating characteristics
#' (analytically) of a randomized comparative phase II oncology trial design
#' that assumes a Bernoulli primary outcome variable.
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
#' \item \code{\link{sim}}: Evaluate the operating characteristics (via
#' simulation) of a randomized comparative phase II oncology trial design
#' that assumes a Bernoulli primary outcome variable.
#' \item \code{\link{summary.ph2rand_des}}: Display a summary of a randomized
#' comparative phase II oncology trial design that assumes a Bernoulli primary
#' outcome variable.
#' \item \code{\link{terminal}}: Find the terminal points of a randomized
#' comparative phase II oncology trial design that assumes a Bernoulli primary
#' outcome variable.
#' }
#'
#' @section References:
#' Jung SH (2008) Randomized phase II trials with a prospective control.
#' \emph{Stat Med} \strong{27}(4)\strong{:}568--83.
#' DOI: \doi{10.1002/sim.2961}.
#' PMID: \href{https://pubmed.ncbi.nlm.nih.gov/17573688}{17573688}.
#' 
#' Jung SH, Sargent DJ (2014) Randomized phase II clinical trials.
#' \emph{J Biopharm Stat} \strong{24}(4)\strong{:}802--16.
#' DOI: \doi{10.1080/10543406.2014.901343}.
#' PMID: \href{https://pubmed.ncbi.nlm.nih.gov/24697589}{24697589}.
#' 
#' Kepner JL (2010) On group sequential designs comparing two binomial
#' proportions.
#' \emph{J Biopharm Stat} \strong{20}(1)\strong{:}145--59.
#' DOI: \doi{10.1080/10543400903280621}.
#' PMID: \href{https://pubmed.ncbi.nlm.nih.gov/20077254}{20077254}.
#' 
#' Litwin S, Basickes S, Ross EA (2017) Two-sample binary phase 2 trials with
#' low type I error and low sample size.
#' \emph{Stat Med} \strong{36}(9)\strong{:}1383--94.
#' DOI: \doi{10.1002/sim.7226}.
#' PMID: \href{https://pubmed.ncbi.nlm.nih.gov/28118686}{28118686}.
#' 
#' Shan G, Ma C, Hutson AD, Wilding GE (2013) Randomized two-stage phase II
#' clinical trial designs based on Barnard's exact test.
#' \emph{J Biopharm Stat} \strong{23}(5)\strong{:}1081--90.
#' DOI: \doi{10.1080/10543406.2013.813525}.
#' PMID: \href{https://pubmed.ncbi.nlm.nih.gov/23957517}{23957517}.
#' 
#' @docType package
#' @name ph2rand
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(".data")
}


<!-- README.md is generated from README.Rmd. Please edit that file -->

# ph2rand <img src='man/figures/ph2rand.png' align="right" height="139" />

*Design of Randomized Comparative Phase II Oncology Trials with a Binary
Primary Outcome*

## Description

**ph2rand** provides a suite of functions to assist in the design of
randomized comparative phase II oncology trials with a binary primary
outcome variable. Specifically, support is provided to: (a) perform a
sample size calculation when using one of several published designs
(Jung, 2008; Jung and Sargent, 2014; Kepner, 2010; Litwin *et al*, 2017,
Shan *et al*, 2013), (b) evaluate the operating characteristics of a
given design, and (c) produce informative plots.

## Getting started

You can install the the latest development version of **ph2rand**,
available from [GitHub](https://github.com/), with:

``` r
devtools::install_github("mjg211/ph2rand")
```

An introductory example of how to make use of the package’s core
functionality can be found below. For further help, please contact
Michael Grayling at <michael.grayling@newcastle.ac.uk>.

## Example

Find a two-stage design from Jung (2008) for the default parameters:

``` r
des <- des_two_stage()
```

Examine its required sample size in each arm, in each stage:

``` r
des$nC
#> [1] 17 17
des$nE
#> [1] 17 17
```

Next, look at its operating characteristics:

``` r
des$opchar
#> # A tibble: 2 x 13
#>     piC   piE `P(pi)` `ESS(pi)` `SDSS(pi)` `MSS(pi)` `E1(pi)` `E2(pi)`
#>   <dbl> <dbl>   <dbl>     <dbl>      <dbl>     <dbl>    <dbl>    <dbl>
#> 1   0.1   0.1  0.0702      47.0       16.5        34        0   0.0702
#> 2   0.1   0.3  0.813       64.7       10.1        68        0   0.813 
#> # … with 5 more variables: `F1(pi)` <dbl>, `F2(pi)` <dbl>, `S1(pi)` <dbl>,
#> #   `S2(pi)` <dbl>, `max N` <int>
```

Compare this to the equivalent design from Litwin *et al* (2017):

``` r
des <- des_two_stage(type  = "single_double",
                     nCmax = 20L)
des$nC
#> [1] 10 10
des$nE
#> [1] 10 10
des$opchar
#> # A tibble: 2 x 13
#>     piC   piE `P(pi)` `ESS(pi)` `SDSS(pi)` `MSS(pi)` `E1(pi)` `E2(pi)`
#>   <dbl> <dbl>   <dbl>     <dbl>      <dbl>     <dbl>    <dbl>    <dbl>
#> 1   0.1   0.1  0.1000      25.2       8.79        20        0   0.1000
#> 2   0.1   0.3  0.804       36.9       7.20        40        0   0.804 
#> # … with 5 more variables: `F1(pi)` <dbl>, `F2(pi)` <dbl>, `S1(pi)` <dbl>,
#> #   `S2(pi)` <dbl>, `max N` <int>
```

## References

Jung SH (2008) Randomized phase II trials with a prospective control.
*Stat Med* 27(4):568–83. DOI:
[10.1002/sim.2961](https://doi.org/10.1002/sim.2961). PMID:
[17573688](https://www.ncbi.nlm.nih.gov/pubmed/17573688).

Jung SH, Sargent DJ (2014) Randomized phase II clinical trials. *J
Biopharm Stat* 24(4):802–16. DOI:
[10.1080/10543406.2014.901343](https://doi.org/10.1080/10543406.2014.901343).
PMID: [24697589](https://www.ncbi.nlm.nih.gov/pubmed/24697589).

Kepner JL (2010) On group sequential designs comparing two binomial
proportions. *J Biopharm Stat* 20(1):145–59. DOI:
[10.1080/10543400903280621](https://doi.org/10.1080/10543400903280621).
PMID: [20077254](https://www.ncbi.nlm.nih.gov/pubmed/20077254).

Litwin S, Basickes S, Ross EA (2017) Two-sample binary phase 2 trials
with low type I error and low sample size. *Stat Med* 36(9):1383–94.
DOI: [10.1002/sim.7226](http://doi.org/10.1002/sim.7226). PMID:
[28118686](https://www.ncbi.nlm.nih.gov/pubmed/28118686).

Shan G, Ma C, Hutson AD, Wilding GE (2013) Randomized two-stage phase II
clinical trial designs based on Barnard’s exact test. *J Biopharm Stat*
23(5):1081–90. DOI:
[10.1080/10543406.2013.813525](https://doi.org/10.1080/10543406.2013.813525).
PMID: [23957517](https://www.ncbi.nlm.nih.gov/pubmed/23957517).


<!-- README.md is generated from README.Rmd. Please edit that file -->

# ph2rand <img src='man/figures/ph2rand.png' align="right" height="139" />

*Design of randomized comparative phase II oncology trials*

## Description

**ph2rand** provides a suite of functions to assist with the design of
randomized comparative phase II oncology trials. Specifically, support
is provided to perform a sample size calculation for each of the most
popular randomized comparative phase II oncology trial designs (Jung,
2008; Jung and Sargent, 2014; Kepner, 2010; Litwin *et al*, 2017, Shan
*et al*, 2013), for both point and composite null hypotheses.

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

    des <- des_two_stage()

## References

Jung SH (2008) [Randomized phase II trials with a prospective
control](https://doi.org/10.1002/sim.2961). *Stat Med* 27(4):568-83.

Jung SH, Sargent DJ (2014) [Randomized phase II clinical
trials](https://doi.org/10.1080/10543406.2014.901343). *J Biopharm Stat*
24(4):802-16.

Kepner JL (2010) [On group sequential designs comparing two binomial
proportions](https://doi.org/10.1080/10543400903280621). *J Biopharm
Stat* 20(1):145-59.

Litwin S, Basickes S, Ross EA (2017) [Two-sample binary phase 2 trials
with low type I error and low sample
size](http://doi.org/10.1002/sim.7226). *Stat Med* 36(9):1383-94.

Shan G, Ma C, Hutson AD, Wilding GE (2013) [Randomized two-stage phase
II clinical trial designs based on Barnard’s exact
test](https://doi.org/10.1080/10543406.2013.813525). *J Biopharm Stat*
23(5):1081-90.

#include <Rcpp.h>
using namespace Rcpp;

#ifndef DBINOM_TWO_STAGE_H
#define DBINOM_TWO_STAGE_H

NumericMatrix dbinom_two_stage(NumericVector pi, NumericVector nC,
                               NumericVector nE);

#endif
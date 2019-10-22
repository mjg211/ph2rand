#include <Rcpp.h>
using namespace Rcpp;

#ifndef DBINOM_ONE_STAGE_H
#define DBINOM_ONE_STAGE_H

NumericMatrix dbinom_one_stage(NumericVector pi, int nC, int nE);

#endif
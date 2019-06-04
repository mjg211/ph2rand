#include <Rcpp.h>
using namespace Rcpp;

#ifndef DBINOM_DES_TWO_STAGE_H
#define DBINOM_DES_TWO_STAGE_H

NumericMatrix dbinom_des_two_stage(NumericMatrix dbinom1, double pi_typeI,
                                   double pi_power, double delta, int n01,
                                   int n02, int n11, int n12);

#endif
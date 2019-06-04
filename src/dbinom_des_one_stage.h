#include <Rcpp.h>
using namespace Rcpp;

#ifndef DBINOM_DES_ONE_STAGE_H
#define DBINOM_DES_ONE_STAGE_H

NumericMatrix dbinom_des_one_stage(double pi_typeI, double pi_power,
                                   double delta, int n0, int n1);

#endif
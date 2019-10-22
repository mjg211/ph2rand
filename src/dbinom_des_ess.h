#include <Rcpp.h>
using namespace Rcpp;

#ifndef DBINOM_DES_ESS_H
#define DBINOM_DES_ESS_H

NumericMatrix dbinom_des_ess(NumericMatrix dbinom, double pi_typeI,
                             double pi_power, double delta, double pi_ess,
                             int nC, int nE);

#endif


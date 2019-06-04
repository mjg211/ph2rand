#include <Rcpp.h>
#include "dbinom_des_ess.h"
#include "dbinom_des_one_stage.h"
#include "dbinom_des_two_stage.h"
#include "dbinom_one_stage.h"
#include "dbinom_two_stage.h"
#include "pi_power_finder.h"
#include "pi_typeI_finder.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix dbinom_des_ess(NumericMatrix dbinom1, double pi_typeI,
                             double pi_power, double delta, double pi_ess,
                             int n0, int n1) {
  NumericMatrix dbinom_ess(3, max(n0, n1) + 1);
  if (pi_ess == pi_typeI) {
    dbinom_ess(0, _)      = dbinom1(0, _);
    dbinom_ess(1, _)      = dbinom1(2, _);
  }
  else if (pi_ess == pi_power) {
    dbinom_ess(0, _)      = dbinom1(1, _);
    if (n1 == n0) {
      dbinom_ess(1, _)    = dbinom_ess(0, _);
    }
    else {
      for (int x1 = 0; x1 <= n1; x1++) {
        dbinom_ess(1, x1) = R::dbinom(x1, n1, pi_ess, 0);
      }
    }
  }
  else if (pi_ess == pi_power + delta) {
    dbinom_ess(1, _)      = dbinom1(3, _);
    if (n0 == n1) {
      dbinom_ess(0, _)    = dbinom_ess(1, _);
    }
    else {
      for (int x0 = 0; x0 <= n0; x0++) {
        dbinom_ess(0, x0) = R::dbinom(x0, n0, pi_ess, 0);
      }
    }
  }
  else {
    for (int x0 = 0; x0 <= n0; x0++) {
      dbinom_ess(0, x0)   = R::dbinom(x0, n0, pi_ess, 0);
    }
    if (n0 == n1) {
      dbinom_ess(1, _)    = dbinom_ess(0, _);
    }
    else {
      for (int x1 = 0; x1 <= n1; x1++) {
        dbinom_ess(1, x1) = R::dbinom(x1, n1, pi_ess, 0);
      }
    }
  }
  if (pi_ess + delta == pi_typeI) {
    dbinom_ess(2, _)      = dbinom1(2, _);
  }
  else if (pi_ess + delta == pi_power) {
    if (n1 == n0) {
      dbinom_ess(2, _)    = dbinom1(1, _);
    }
    else {
      for (int x1 = 0; x1 <= n1; x1++) {
        dbinom_ess(2, x1) = R::dbinom(x1, n1, pi_power, 0);
      }
    }
  }
  else if (pi_ess == pi_power) {
    dbinom_ess(2, _)      = dbinom1(3, _);
  }
  else {
    for (int x1 = 0; x1 <= n1; x1++) {
      dbinom_ess(2, x1)   = R::dbinom(x1, n1, pi_ess + delta, 0);
    }
  }
  return dbinom_ess;
}

// [[Rcpp::export]]
NumericMatrix dbinom_des_one_stage(double pi_typeI, double pi_power,
                                   double delta, int n0, int n1) {
  NumericMatrix dbinom(4, max(n0, n1) + 1);
  for (int x0 = 0; x0 <= n0; x0++) {
    dbinom(0, x0)   = R::dbinom(x0, n0, pi_typeI, 0);
  }
  if (pi_typeI == pi_power) {
    dbinom(1, _)    = dbinom(0, _);
  }
  else {
    for (int x0 = 0; x0 <= n0; x0++) {
      dbinom(1, x0) = R::dbinom(x0, n0, pi_power, 0);
    }
  }
  if (n0 == n1) {
    dbinom(2, _)    = dbinom(0, _);
  }
  else {
    for (int x1 = 0; x1 <= n1; x1++) {
      dbinom(2, x1) = R::dbinom(x1, n1, pi_typeI, 0);
    }
  }
  if (pi_typeI == pi_power + delta) {
    dbinom(3, _)    = dbinom(2, _);
  }
  else {
    for (int x1 = 0; x1 <= n1; x1++) {
      dbinom(3, x1) = R::dbinom(x1, n1, pi_power + delta, 0);
    }
  }
  return dbinom;
}

// [[Rcpp::export]]
NumericMatrix dbinom_des_two_stage(NumericMatrix dbinom1, double pi_typeI,
                                   double pi_power, double delta, int n01,
                                   int n02, int n11, int n12) {
  NumericMatrix dbinom2(4, max(n02, n12) + 1);
  if (n01 == n02) {
    dbinom2(0, _)     = dbinom1(0, _);
  }
  else if (n11 == n02) {
    dbinom2(0, _)     = dbinom1(2, _);
  }
  else {
    for (int x02 = 0; x02 <= n02; x02++) {
      dbinom2(0, x02) = R::dbinom(x02, n02, pi_typeI, 0);
    }
  }
  if (n01 == n02) {
    dbinom2(1, _)     = dbinom1(1, _);
  }
  else if (pi_typeI == pi_power) {
    dbinom2(1, _)     = dbinom2(0, _);
  }
  else {
    for (int x02 = 0; x02 <= n02; x02++) {
      dbinom2(1, x02) = R::dbinom(x02, n02, pi_power, 0);
    }
  }
  if (n12 == n11) {
    dbinom2(2, _)     = dbinom1(2, _);
  }
  else if (n12 == n01) {
    dbinom2(2, _)     = dbinom1(0, _);
  }
  else if (n12 == n02) {
    dbinom2(2, _)     = dbinom2(0, _);
  }
  else {
    for (int x12 = 0; x12 <= n12; x12++) {
      dbinom2(2, x12) = R::dbinom(x12, n12, pi_typeI, 0);
    }
  }
  if (n12 == n11) {
    dbinom2(3, _)     = dbinom1(3, _);
  }
  else if (pi_typeI == pi_power + delta) {
    if (n01 == n12) {
      dbinom2(3, _) = dbinom1(0, _);
    }
    else if (n02 == n12) {
      dbinom2(3, _) = dbinom2(0, _);
    }
  }
  else {
    for (int x12 = 0; x12 <= n12; x12++) {
      dbinom2(3, x12) = R::dbinom(x12, n12, pi_power + delta, 0);
    }
  }
  return dbinom2;
}

// [[Rcpp::export]]
NumericMatrix dbinom_one_stage(NumericVector pi, int n0, int n1) {
  NumericMatrix dbinom(2, max(n0, n1) + 1);
  for (int x0 = 0; x0 <= n0; x0++) {
    dbinom(0, x0)   = R::dbinom(x0, n0, pi[0], 0);
  }
  if ((pi[1] == pi[0]) && (n0 == n1)) {
    dbinom(1, _)    = dbinom(0, _);
  }
  else {
    for (int x1 = 0; x1 <= n1; x1++) {
      dbinom(1, x1) = R::dbinom(x1, n1, pi[1], 0);
    }
  }
  return dbinom;
}

// [[Rcpp::export]]
NumericMatrix dbinom_two_stage(NumericVector pi, NumericVector n0,
                               NumericVector n1) {
  NumericMatrix dbinom(4, max(NumericVector::create(n0[0], n0[1],
                                                    n1[0], n1[1])) + 1);
  for (int x01 = 0; x01 <= n0[0]; x01++) {
    dbinom(0, x01)     = R::dbinom(x01, n0[0], pi[0], 0);
  }
  if (n0[0] == n0[1]) {
    dbinom(1, _)       = dbinom(0, _);
  }
  else {
    for (int x02 = 0; x02 <= n0[1]; x02++) {
      dbinom(1, x02)   = R::dbinom(x02, n0[1], pi[0], 0);
    }
  }
  if (pi[0] == pi[1]) {
    if (n1[0] == n0[0]) {
      dbinom(2, _)     = dbinom(0, _);
    }
    else if (n1[0] == n0[1]) {
      dbinom(2, _)     = dbinom(1, _);
    }
    else {
      for (int x11 = 0; x11 <= n1[0]; x11++) {
        dbinom(2, x11) = R::dbinom(x11, n1[0], pi[1], 0);
      }
    }
    if (n1[1] == n1[0]) {
      dbinom(3, _)     = dbinom(2, _);
    }
    else if (n1[1] == n0[0]) {
      dbinom(3, _)     = dbinom(0, _);
    }
    else if (n1[1] == n0[1]) {
      dbinom(3, _)     = dbinom(1, _);
    }
    else {
      for (int x12 = 0; x12 <= n1[1]; x12++) {
        dbinom(3, x12) = R::dbinom(x12, n1[1], pi[1], 0);
      }
    }
  }
  else {
    for (int x11 = 0; x11 <= n1[0]; x11++) {
      dbinom(2, x11)   = R::dbinom(x11, n1[0], pi[1], 0);
    }
    if (n1[1] == n1[0]) {
      dbinom(3, _)     = dbinom(2, _);
    }
    else {
      for (int x12 = 0; x12 <= n1[1]; x12++) {
        dbinom(3, x12) = R::dbinom(x12, n1[1], pi[1], 0);
      }
    }
  }
  return dbinom;
}

// [[Rcpp::export]]
double pi_power_finder(int point_alt, NumericVector pi_alt, double delta) {
  double pi;
  if (point_alt == 1) {
    pi   = pi_alt[0];
  }
  else {
    if ((pi_alt[0] <= 0.5 - 0.5*delta) && (pi_alt[1] >= 0.5 - 0.5*delta)) {
      pi = 0.5 - 0.5*delta;
    }
    else {
      pi = 0.5*(pi_alt[0] + pi_alt[1]);
    }
  }
  return pi;
}

// [[Rcpp::export]]
double pi_typeI_finder(int point_null, NumericVector pi_null) {
  double pi;
  if (point_null == 1) {
    pi   = pi_null[0];
  }
  else {
    if ((pi_null[0] <= 0.5) && (pi_null[1] >= 0.5)) {
      pi = 0.5;
    }
    else {
      pi = 0.5*(pi_null[0] + pi_null[1]);
    }
  }
  return pi;
}
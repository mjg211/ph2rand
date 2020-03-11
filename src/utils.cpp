#include <Rcpp.h>
#include "dbinom_des_ess.h"
#include "dbinom_des_one_stage.h"
#include "dbinom_des_two_stage.h"
#include "dbinom_one_stage.h"
#include "dbinom_two_stage.h"
#include "message_cpp.h"
#include "pi_power_finder.h"
#include "pi_typeI_finder.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix dbinom_des_ess(NumericMatrix dbinom1, double pi_typeI,
                             double pi_power, double delta, double pi_ess,
                             int nC, int nE) {
  NumericMatrix dbinom_ess(3, max(nC, nE) + 1);
  if (pi_ess == pi_typeI) {
    dbinom_ess(0, _)      = dbinom1(0, _);
    dbinom_ess(1, _)      = dbinom1(2, _);
  }
  else if (pi_ess == pi_power) {
    dbinom_ess(0, _)      = dbinom1(1, _);
    if (nC == nE) {
      dbinom_ess(1, _)    = dbinom_ess(0, _);
    }
    else {
      for (int xE = 0; xE <= nE; xE++) {
        dbinom_ess(1, xE) = R::dbinom(xE, nE, pi_ess, 0);
      }
    }
  }
  else if (pi_ess == pi_power + delta) {
    dbinom_ess(1, _)      = dbinom1(3, _);
    if (nC == nE) {
      dbinom_ess(0, _)    = dbinom_ess(1, _);
    }
    else {
      for (int xC = 0; xC <= nC; xC++) {
        dbinom_ess(0, xC) = R::dbinom(xC, nC, pi_ess, 0);
      }
    }
  }
  else {
    for (int xC = 0; xC <= nC; xC++) {
      dbinom_ess(0, xC)   = R::dbinom(xC, nC, pi_ess, 0);
    }
    if (nC == nE) {
      dbinom_ess(1, _)    = dbinom_ess(0, _);
    }
    else {
      for (int xE = 0; xE <= nE; xE++) {
        dbinom_ess(1, xE) = R::dbinom(xE, nE, pi_ess, 0);
      }
    }
  }
  if (pi_ess + delta == pi_typeI) {
    dbinom_ess(2, _)      = dbinom1(2, _);
  }
  else if (pi_ess + delta == pi_power) {
    if (nC == nE) {
      dbinom_ess(2, _)    = dbinom1(1, _);
    }
    else {
      for (int xE = 0; xE <= nE; xE++) {
        dbinom_ess(2, xE) = R::dbinom(xE, nE, pi_power, 0);
      }
    }
  }
  else if (pi_ess == pi_power) {
    dbinom_ess(2, _)      = dbinom1(3, _);
  }
  else {
    for (int xE = 0; xE <= nE; xE++) {
      dbinom_ess(2, xE)   = R::dbinom(xE, nE, pi_ess + delta, 0);
    }
  }
  return dbinom_ess;
}

// [[Rcpp::export]]
NumericMatrix dbinom_des_one_stage(double pi_typeI, double pi_power,
                                   double delta, int nC, int nE) {
  NumericMatrix dbinom(4, max(nC, nE) + 1);
  for (int xC = 0; xC <= nC; xC++) {
    dbinom(0, xC)   = R::dbinom(xC, nC, pi_typeI, 0);
  }
  if (pi_typeI == pi_power) {
    dbinom(1, _)    = dbinom(0, _);
  }
  else {
    for (int xC = 0; xC <= nC; xC++) {
      dbinom(1, xC) = R::dbinom(xC, nC, pi_power, 0);
    }
  }
  if (nC == nE) {
    dbinom(2, _)    = dbinom(0, _);
  }
  else {
    for (int xE = 0; xE <= nE; xE++) {
      dbinom(2, xE) = R::dbinom(xE, nE, pi_typeI, 0);
    }
  }
  if (pi_typeI == pi_power + delta) {
    dbinom(3, _)    = dbinom(2, _);
  }
  else {
    for (int xE = 0; xE <= nE; xE++) {
      dbinom(3, xE) = R::dbinom(xE, nE, pi_power + delta, 0);
    }
  }
  return dbinom;
}

// [[Rcpp::export]]
NumericMatrix dbinom_des_two_stage(NumericMatrix dbinom1, double pi_typeI,
                                   double pi_power, double delta, int n1C,
                                   int n2C, int n1E, int n2E) {
  NumericMatrix dbinom2(4, max(n2C, n2E) + 1);
  if (n1C == n2C) {
    dbinom2(0, _)     = dbinom1(0, _);
  }
  else if (n1E == n2C) {
    dbinom2(0, _)     = dbinom1(2, _);
  }
  else {
    for (int xC2 = 0; xC2 <= n2C; xC2++) {
      dbinom2(0, xC2) = R::dbinom(xC2, n2C, pi_typeI, 0);
    }
  }
  if (n1C == n2C) {
    dbinom2(1, _)     = dbinom1(1, _);
  }
  else if (pi_typeI == pi_power) {
    dbinom2(1, _)     = dbinom2(0, _);
  }
  else {
    for (int xC2 = 0; xC2 <= n2C; xC2++) {
      dbinom2(1, xC2) = R::dbinom(xC2, n2C, pi_power, 0);
    }
  }
  if (n2E == n1E) {
    dbinom2(2, _)     = dbinom1(2, _);
  }
  else if (n2E == n1C) {
    dbinom2(2, _)     = dbinom1(0, _);
  }
  else if (n2E == n2C) {
    dbinom2(2, _)     = dbinom2(0, _);
  }
  else {
    for (int xE2 = 0; xE2 <= n2E; xE2++) {
      dbinom2(2, xE2) = R::dbinom(xE2, n2E, pi_typeI, 0);
    }
  }
  if (n2E == n1E) {
    dbinom2(3, _)     = dbinom1(3, _);
  }
  else if (pi_typeI == pi_power + delta) {
    if (n1C == n2E) {
      dbinom2(3, _) = dbinom1(0, _);
    }
    else if (n2C == n2E) {
      dbinom2(3, _) = dbinom2(0, _);
    }
  }
  else {
    for (int xE2 = 0; xE2 <= n2E; xE2++) {
      dbinom2(3, xE2) = R::dbinom(xE2, n2E, pi_power + delta, 0);
    }
  }
  return dbinom2;
}

// [[Rcpp::export]]
NumericMatrix dbinom_one_stage(NumericVector pi, int nC, int nE) {
  NumericMatrix dbinom(2, max(nC, nE) + 1);
  for (int xC = 0; xC <= nC; xC++) {
    dbinom(0, xC)   = R::dbinom(xC, nC, pi[0], 0);
  }
  if ((pi[1] == pi[0]) && (nC == nE)) {
    dbinom(1, _)    = dbinom(0, _);
  }
  else {
    for (int xE = 0; xE <= nE; xE++) {
      dbinom(1, xE) = R::dbinom(xE, nE, pi[1], 0);
    }
  }
  return dbinom;
}

// [[Rcpp::export]]
NumericMatrix dbinom_two_stage(NumericVector pi, NumericVector nC,
                               NumericVector nE) {
  NumericMatrix dbinom(4, max(NumericVector::create(nC[0], nC[1],
                                                    nE[0], nE[1])) + 1);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    dbinom(0, xC1)     = R::dbinom(xC1, nC[0], pi[0], 0);
  }
  if (nC[0] == nC[1]) {
    dbinom(1, _)       = dbinom(0, _);
  }
  else {
    for (int xC2 = 0; xC2 <= nC[1]; xC2++) {
      dbinom(1, xC2)   = R::dbinom(xC2, nC[1], pi[0], 0);
    }
  }
  if (pi[0] == pi[1]) {
    if (nE[0] == nC[0]) {
      dbinom(2, _)     = dbinom(0, _);
    }
    else if (nE[0] == nC[1]) {
      dbinom(2, _)     = dbinom(1, _);
    }
    else {
      for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
        dbinom(2, xE1) = R::dbinom(xE1, nE[0], pi[1], 0);
      }
    }
    if (nE[1] == nE[0]) {
      dbinom(3, _)     = dbinom(2, _);
    }
    else if (nE[1] == nC[0]) {
      dbinom(3, _)     = dbinom(0, _);
    }
    else if (nE[1] == nC[1]) {
      dbinom(3, _)     = dbinom(1, _);
    }
    else {
      for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
        dbinom(3, xE2) = R::dbinom(xE2, nE[1], pi[1], 0);
      }
    }
  }
  else {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      dbinom(2, xE1)   = R::dbinom(xE1, nE[0], pi[1], 0);
    }
    if (nE[1] == nE[0]) {
      dbinom(3, _)     = dbinom(2, _);
    }
    else {
      for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
        dbinom(3, xE2) = R::dbinom(xE2, nE[1], pi[1], 0);
      }
    }
  }
  return dbinom;
}

// [[Rcpp::export]]
void message_cpp(std::string text_1, std::string text_2) {
  Rcpp::Function msg("message");
  msg(std::string("  [from Rcpp: ") + std::string(text_1) +
    std::string(text_2) + std::string("]"));
}

// [[Rcpp::export]]
double pi_power_finder(int point_alt, NumericVector Pi1, double delta) {
  double pi;
  if (point_alt == 1) {
    pi   = Pi1[0];
  }
  else {
    if ((Pi1[0] <= 0.5 - 0.5*delta) && (Pi1[1] >= 0.5 - 0.5*delta)) {
      pi = 0.5 - 0.5*delta;
    }
    else {
      pi = 0.5*(Pi1[0] + Pi1[1]);
    }
  }
  return pi;
}

// [[Rcpp::export]]
double pi_typeI_finder(int point_null, NumericVector Pi0) {
  double pi;
  if (point_null == 1) {
    pi   = Pi0[0];
  }
  else {
    if ((Pi0[0] <= 0.5) && (Pi0[1] >= 0.5)) {
      pi = 0.5;
    }
    else {
      pi = 0.5*(Pi0[0] + Pi0[1]);
    }
  }
  return pi;
}
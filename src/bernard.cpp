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
NumericMatrix bernard_pmf_two_stage_cpp(NumericVector pi, NumericVector nC,
                                        NumericVector nE, double e1, double f1,
                                        double e2, NumericVector k) {
  
  int           counter                    = 0,
                sum_nC                     = sum(nC),
                sum_nE                     = sum(nE);
  double        fact,
                statistic;
  NumericMatrix prob_x1(nC[0] + 1, nE[0] + 1),
                poss_x2(sum_nC + 1, sum_nE + 1),
                prob_x2(sum_nC + 1, sum_nE + 1),
                pmf((nC[0] + 1)*(nC[1] + 1)*(nE[0] + 1)*(nE[1] + 1), 8),
                dbinom                     = dbinom_two_stage(pi, nC, nE);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      prob_x1(xC1, xE1)                    = dbinom(0, xC1)*dbinom(2, xE1);
      if (((xC1 == 0) && (xE1 == 0)) || ((xC1 == nC[0]) && (xE1 == nE[0]))) {
        statistic                          = 0;
      }
      else {
        fact                               = (xC1 + xE1)/(nC[0] + nE[0]);
        statistic                          =
          (xE1/nE[0] - xC1/nC[0])/sqrt(fact*(1 - fact)*(1/nC[0] + 1/nE[0]));
      }
      if (((statistic >= e1) || (statistic <= f1)) && (k[0] == 1)) {
        pmf(counter, _)                    =
          NumericVector::create(xC1, xE1, nC[0], nE[0], statistic,
                                (statistic >= e1), 1, prob_x1(xC1, xE1));
        counter++;
      }
      else if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
        for (int xC2 = 0; xC2 <= nC[1]; xC2++) {
          for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
            poss_x2(xC1 + xC2, xE1 + xE2) += 1;
            prob_x2(xC1 + xC2, xE1 + xE2) += prob_x1(xC1, xE1)*dbinom(1, xC2)*
              dbinom(3, xE2);
          }
        }
      }
    }
  }
  if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
    for (int xC = 0; xC <= sum_nC; xC++) {
      for (int xE = 0; xE <= sum_nE; xE++) {
        if (poss_x2(xC, xE) > 0) {
          if (((xC == 0) && (xE == 0)) || ((xC == sum_nC) && (xE == sum_nE))) {
            statistic                     = 0;
          }
          else {
            fact                          = (xC + xE)/(sum_nC + sum_nE);
            statistic                     =
              (xE/sum_nE - xC/sum_nC)/sqrt(fact*(1 - fact)*
              (1/sum_nC + 1/sum_nE));
          }
          pmf(counter, _)                 =
            NumericVector::create(xC, xE, sum_nC, sum_nE, statistic,
                                  (statistic >= e2), 2, prob_x2(xC, xE));
          counter++;
        }
      }
    }
  }
  NumericMatrix output                    = pmf(Range(0, counter - 1),
                                                Range(0, 7));
  return output;
}

// [[Rcpp::export]]
double bernard_power_one_stage(NumericVector pi, int nC, int nE, double e,
                               NumericMatrix poss_x, NumericMatrix poss_B) {
  double        power  = 0;
  NumericMatrix dbinom = dbinom_one_stage(pi, nC, nE);
  for (int o = 0; o <= (nC + 1)*(nE + 1) - 1; o++) {
    if (poss_B(poss_x(o, 0), poss_x(o, 1)) >= e) {
      power           += dbinom(0, poss_x(o, 0))*dbinom(1, poss_x(o, 1));
    }
  }
  return power;
}

// [[Rcpp::export]]
double bernard_power_two_stage(NumericVector pi, NumericVector nC,
                               NumericVector nE, NumericVector e,
                               NumericVector f, List poss_x, List poss_B) {
  double        prob_x1,
                power   = 0;
  NumericMatrix poss_B1 = poss_B[0],
                poss_B2 = poss_B[1],
                poss_x1 = poss_x[0],
                poss_x2 = poss_x[1],
                dbinom  = dbinom_two_stage(pi, nC, nE);
  if ((f[0] >= poss_B1(nC[0], 0)) && (e[0] <= poss_B1(0, nE[0]))) {
    for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
      if (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) >= e[0]) {
        power          += dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
      }
      else if (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) > f[0]) {
        prob_x1         = dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
        for (int o2 = 0; o2 <= (nC[1] + 1)*(nE[1] + 1) - 1; o2++) {
          if (poss_B2(poss_x1(o1, 0) + poss_x2(o1, 0),
                      poss_x1(o1, 1) + poss_x2(o1, 1)) >= e[1]) {
            power      += prob_x1*dbinom(1, poss_x2(o1, 0))*
                            dbinom(3, poss_x2(o1, 1));
          }
        }
      }
    }
  }
  else if (f[0] >= poss_B1(nC[0], 0)) {
    for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
      if (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) > f[0]) {
        prob_x1         = dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
        for (int o2 = 0; o2 <= (nC[1] + 1)*(nE[1] + 1) - 1; o2++) {
          if (poss_B2(poss_x1(o1, 0) + poss_x2(o1, 0),
                      poss_x1(o1, 1) + poss_x2(o1, 1)) >= e[1]) {
            power      += prob_x1*dbinom(1, poss_x2(o1, 0))*
                            dbinom(3, poss_x2(o1, 1));
          }
        }
      }
    }
  }
  else {
    for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
      if (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) >= e[0]) {
        power          += dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
      }
      else {
        prob_x1         = dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
        for (int o2 = 0; o2 <= (nC[1] + 1)*(nE[1] + 1) - 1; o2++) {
          if (poss_B2(poss_x1(o1, 0) + poss_x2(o1, 0),
                      poss_x1(o1, 1) + poss_x2(o1, 1)) >= e[1]) {
            power      += prob_x1*dbinom(1, poss_x2(o1, 0))*
                            dbinom(3, poss_x2(o1, 1));
          }
        }
      }
    }
  }
  return power;
}

// [[Rcpp::export]]
NumericMatrix bernard_terminal_two_stage_cpp(NumericVector nC, NumericVector nE,
                                             double e1, double f1, double e2,
                                             NumericVector k) {
  
  int           xC,
                xE,
                counter            = 0,
                sum_nC             = sum(nC),
                sum_nE             = sum(nE);
  double        fact,
                statistic;
  NumericMatrix x2_mat(sum_nC + 1, sum_nE + 1),
                terminal((nC[0] + 1)*(nC[1] + 1)*(nE[0] + 1)*(nE[1] + 1), 7);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      if (((xC1 == 0) && (xE1 == 0)) || ((xC1 == nC[0]) && (xE1 == nE[0]))) {
        statistic                  = 0;
      }
      else {
        fact                       = (xC1 + xE1)/(nC[0] + nE[0]);
        statistic                  =
          (xE1/nE[0] - xC1/nC[0])/sqrt(fact*(1 - fact)*(1/nC[0] + 1/nE[0]));
      }
      if (((statistic >= e1) || (statistic <= f1)) && (k[0] == 1)) {
        terminal(counter, _)       =
          NumericVector::create(xC1, xE1, nC[0], nE[0], statistic,
                                (statistic >= e1), 1);
        counter++;
      }
      else if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
        for (int xC2 = 0; xC2 <= nC[1]; xC2++) {
          for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
            xC                     = xC1 + xC2;
            xE                     = xE1 + xE2;
            if (x2_mat(xC, xE) == 0) {
              if (((xC == 0) && (xE == 0)) ||
                    ((xC == sum_nC) && (xE == sum_nE))) {
                statistic          = 0;
              }
              else {
                fact               = (xC + xE)/(sum_nC + sum_nE);
                statistic          =
                  (xE/sum_nE - xC/sum_nC)/
                    sqrt(fact*(1 - fact)*(1/sum_nC + 1/sum_nE));
              }
              terminal(counter, _) =
                NumericVector::create(xC, xE, sum_nC, sum_nE, statistic,
                                      (statistic >= e2), 1);
              x2_mat(xC, xE)      += 1;
              counter++;
            }
          }
        }
      }
    }
  }
  NumericMatrix output             = terminal(Range(0, counter - 1),
                                              Range(0, 6));
  return output;
}

// [[Rcpp::export]]
NumericVector bernard_max_typeI(int J, double alpha, NumericVector nC,
                                NumericVector nE, NumericVector e,
                                NumericVector f, List poss_x, List poss_B,
                                NumericVector Pi0, int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = Pi0[0],
                x_right = Pi0[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      -bernard_power_one_stage(NumericVector::create(v, v), nC[0], nE[0], e[0],
                               poss_x[0], poss_B[0]);
  }
  else {
    f_v                 =
      -bernard_power_two_stage(NumericVector::create(v, v), nC, nE, e, f,
                               poss_x, poss_B);
  }
  NumericVector output(4);
  if ((-f_v > alpha) && (check == 1)) {
    output              = NumericVector::create(v, -f_v, 2, iter);
    return output;
  }
  double f_z,
         z              = pi_typeI_finder(0, Pi0);
  if (J == 1) {
    f_z                 =
      -bernard_power_one_stage(NumericVector::create(z, z), nC[0], nE[0], e[0],
                               poss_x[0], poss_B[0]);
  }
  else {
    f_z                 =
      -bernard_power_two_stage(NumericVector::create(z, z), nC, nE, e, f,
                               poss_x, poss_B);
  }
  if ((-f_z > alpha) && (check == 1)) {
    output              = NumericVector::create(z, -f_z, 2, iter);
    return output;
  }
  double f_u,
         midpoint,
         u,
         w              = v,
         f_w            = f_v,
         w_left         = z - x_left,
         w_right        = x_right - z,
         tol            = 1e-4*z,
         two_tol        = 2*tol,
         d              = 0,
         ed             = 0,
         p              = 0,
         q              = 0,
         r              = 0;
  while (iter < 100) {
    midpoint            = 0.5*(x_left + x_right);
    if (abs(z - midpoint) <= two_tol - 0.5*(x_right - x_left)) {
      output            = NumericVector::create(z, -f_z, 0, iter);
      return output;
    }
    if (abs(ed) > tol) {
      r                 = (z - w)*(f_z - f_v);
      q                 = (z - v)*(f_z - f_w);
      p                 = (z - v)*q - (z - w)*r;
      q                 = 2*(q - r);
      if (q > 0) {
        p               = -p;
      }
      else {
        q               = -q;
      }
      r                 = ed;
      ed                = d;
    }
    if ((abs(p) < abs(0.5*q*r)) && (p < q*w_left) && (p < q*w_right)) {
      d                 = p/q;
      u                 = z + d;
      if ((u - x_left < two_tol) || (x_right - u < two_tol)) {
        d               = (z < midpoint ? tol : -tol);
      }
    }
    else {
      ed                = (z < midpoint ? x_right - z : -(z - x_left));
      d                 = golden*ed;
    }
    if (abs(d) >= tol) {
      u                 = z + d;
    }
    else {
      u                 = z + (d > 0 ? tol : -tol);
    }
    if (J == 1) {
      f_u               =
        -bernard_power_one_stage(NumericVector::create(u, u), nC[0], nE[0],
                                 e[0], poss_x[0], poss_B[0]);
    }
    else {
      f_u               =
        -bernard_power_two_stage(NumericVector::create(u, u), nC, nE, e, f,
                                 poss_x, poss_B);
    }
    if ((-f_u > alpha) && (check == 1)) {
      output            = NumericVector::create(u, -f_u, 2, iter);
      return output;
    }
    if (f_u <= f_z) {
      if (u < z) {
        x_right         = z;
      }
      else {
        x_left          = z;
      }
      v                 = w;
      f_v               = f_w;
      w                 = z;
      f_w               = f_z;
      z                 = u;
      f_z               = f_u;
    }
    else {
      if (u < z) {
        x_left          = u;
      }
      else {
        x_right         = u;
      }
      if ((f_u <= f_w) || (w == z)) {
        v               = w;
        f_v             = f_w;
        w               = u;
        f_w             = f_u;
      }
      else if ((f_u <= f_v) || (v == z) || (v == w)) {
        v               = u;
        f_v             = f_u;
      }
    }
    w_left              = z - x_left;
    w_right             = x_right - z;
    iter++;
  }
  output                = NumericVector::create(z, -f_z, 1, iter);
  return output;
}

// [[Rcpp::export]]
NumericVector bernard_min_power(int J, double beta, double delta,
                                NumericVector nC, NumericVector nE,
                                NumericVector e, NumericVector f, List poss_x,
                                List poss_B, NumericVector Pi1, int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = Pi1[0],
                x_right = Pi1[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      bernard_power_one_stage(NumericVector::create(v, v + delta), nC[0],
                              nE[0], e[0], poss_x[0], poss_B[0]);
  }
  else {
    f_v                 =
      bernard_power_two_stage(NumericVector::create(v, v + delta), nC, nE, e,
                              f, poss_x, poss_B);
  }
  NumericVector output(4);
  if (f_v < 1 - beta && check == 1) {
    output              = NumericVector::create(v, f_v, 2, iter);
    return output;
  }
  double f_z,
         z              = pi_power_finder(0, Pi1, delta);
  if (J == 1) {
    f_z                 =
      bernard_power_one_stage(NumericVector::create(z, z + delta), nC[0],
                              nE[0], e[0], poss_x[0], poss_B[0]);
  }
  else {
    f_z                 =
      bernard_power_two_stage(NumericVector::create(z, z + delta), nC, nE, e,
                              f, poss_x, poss_B);
  }
  if ((f_z < 1 - beta) && (check == 1)) {
    output              = NumericVector::create(z, f_z, 2, iter);
    return output;
  }
  double f_u,
         midpoint,
         u,
         w              = v,
         f_w            = f_v,
         w_left         = z - x_left,
         w_right        = x_right - z,
         tol            = 1e-4*z,
         two_tol        = 2*tol,
         d              = 0,
         ed             = 0,
         p              = 0,
         q              = 0,
         r              = 0;
  while (iter < 100) {
    midpoint            = 0.5*(x_left + x_right);
    if (abs(z - midpoint) <= two_tol - 0.5*(x_right - x_left)) {
      output            = NumericVector::create(z, f_z, 0, iter);
      return output;
    }
    if (abs(ed) > tol) {
      r                 = (z - w)*(f_z - f_v);
      q                 = (z - v)*(f_z - f_w);
      p                 = (z - v)*q - (z - w)*r;
      q                 = 2*(q - r);
      if (q > 0) {
        p               = -p;
      }
      else {
        q               = -q;
      }
      r                 = ed;
      ed                = d;
    }
    if ((abs(p) < abs(0.5*q*r)) && (p < q*w_left) && (p < q*w_right)) {
      d                 = p/q;
      u                 = z + d;
      if ((u - x_left < two_tol) || (x_right - u < two_tol)) {
        d               = (z < midpoint ? tol : -tol);
      }
    }
    else {
      ed                = (z < midpoint ? x_right - z : -(z - x_left));
      d                 = golden*ed;
    }
    if (abs(d) >= tol) {
      u                 = z + d;
    }
    else {
      u                 = z + (d > 0 ? tol : -tol);
    }
    if (J == 1) {
      f_u               =
        bernard_power_one_stage(NumericVector::create(u, u + delta), nC[0],
                                nE[0], e[0], poss_x[0], poss_B[0]);
    }
    else {
      f_u               =
        bernard_power_two_stage(NumericVector::create(u, u + delta), nC, nE, e,
                                f, poss_x, poss_B);
    }
    if ((f_u < 1 - beta) && (check == 1)) {
      output            = NumericVector::create(u, f_u, 2, iter);
      return output;
    }
    if (f_u <= f_z) {
      if (u < z) {
        x_right         = z;
      }
      else {
        x_left          = z;
      }
      v                 = w;
      f_v               = f_w;
      w                 = z;
      f_w               = f_z;
      z                 = u;
      f_z               = f_u;
    }
    else {
      if (u < z) {
        x_left          = u;
      }
      else {
        x_right         = u;
      }
      if ((f_u <= f_w) || (w == z)) {
        v               = w;
        f_v             = f_w;
        w               = u;
        f_w             = f_u;
      }
      else if ((f_u <= f_v) || (v == z) || (v == w)) {
        v               = u;
        f_v             = f_u;
      }
    }
    w_left              = z - x_left;
    w_right             = x_right - z;
    iter++;
  }
  output                = NumericVector::create(z, f_z, 1, iter);
  return output;
}

// [[Rcpp::export]]
double bernard_ess_two_stage(NumericVector pi, NumericVector nC,
                             NumericVector nE, double e1, double f1) {
  double        statistic,
  fact,
  S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, nC[0], nE[0]);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      if (((xC1 == 0) && (xE1 == 0)) || ((xC1 == nC[0]) && (xE1 == nE[0]))) {
        statistic       = 0;
      }
      else {
        fact            = (xC1 + xE1)/(nC[0] + nE[0]);
        statistic       = (xE1/nE[0] - xC1/nC[0])/
          sqrt(fact*(1 - fact)*(1/nC[0] + 1/nE[0]));
      }
      if ((statistic <= f1) | (statistic >= e1)) {
        S1             += dbinom1(0, xC1)*dbinom1(1, xE1);
      }
    }
  }
  double ess           = nC[0] + nE[0] + (1 - S1)*(nC[1] + nE[1]);
  return ess;
}

// [[Rcpp::export]]
double bernard_des_ess_two_stage(NumericVector pi, NumericVector nC,
                                 NumericVector nE, double e1, double f1,
                                 NumericMatrix poss_x, NumericMatrix poss_B) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, nC[0], nE[0]);
  for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
    if ((poss_B(poss_x(o1, 0), poss_x(o1, 1)) <= f1) ||
          (poss_B(poss_x(o1, 0), poss_x(o1, 1)) >= e1)) {
      S1               += dbinom1(0, poss_x(o1, 0))*dbinom1(1, poss_x(o1, 1));
    }
  }
  double ess            = nC[0] + nE[0] + (1 - S1)*(nC[1] + nE[1]);
  return ess;
}

// [[Rcpp::export]]
NumericVector bernard_max_ess_1d_two_stage(NumericVector nC, NumericVector nE,
                                           double e1, double f1,
                                           NumericMatrix poss_x1,
                                           NumericMatrix poss_B1) {
  int           iter    = 0;
  double        f_u,
                midpoint,
                u,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = 0,
                x_right = 1,
                v       = x_left + golden*(x_right - x_left),
                f_v     =
                  -bernard_des_ess_two_stage(NumericVector::create(v, v), nC,
                                             nE, e1, f1, poss_x1, poss_B1),
                z       = 0.5,
                f_z     =
                  -bernard_des_ess_two_stage(NumericVector::create(z, z), nC,
                                             nE, e1, f1, poss_x1, poss_B1),
                w       = v,
                f_w     = f_v,
                w_left  = z - x_left,
                w_right = x_right - z,
                tol     = 1e-4*z,
                two_tol = 2*tol,
                d       = 0,
                ed      = 0,
                p       = 0,
                q       = 0,
                r       = 0;
  NumericVector output(4);
  while (iter < 100) {
    midpoint            = 0.5*(x_left + x_right);
    if (abs(z - midpoint) <= two_tol - 0.5*(x_right - x_left)) {
      output            = NumericVector::create(z, -f_z, 0, iter);
      return output;
    }
    if (abs(ed) > tol) {
      r                 = (z - w)*(f_z - f_v);
      q                 = (z - v)*(f_z - f_w);
      p                 = (z - v)*q - (z - w)*r;
      q                 = 2*(q - r);
      if (q > 0) {
        p               = -p;
      }
      else {
        q               = -q;
      }
      r                 = ed;
      ed                = d;
    }
    if ((abs(p) < abs(0.5*q*r)) && (p < q*w_left) && (p < q*w_right)) {
      d                 = p/q;
      u                 = z + d;
      if ((u - x_left < two_tol) || (x_right - u < two_tol)) {
        d               = (z < midpoint ? tol : -tol);
      }
    }
    else {
      ed                = (z < midpoint ? x_right - z : -(z - x_left));
      d                 = golden*ed;
    }
    if (abs(d) >= tol) {
      u                 = z + d;
    }
    else {
      u                 = z + (d > 0 ? tol : -tol);
    }
    f_u                 =
      -bernard_des_ess_two_stage(NumericVector::create(u, u), nC, nE, e1, f1,
                                 poss_x1, poss_B1);
    if (f_u <= f_z) {
      if (u < z) {
        x_right         = z;
      }
      else {
        x_left          = z;
      }
      v                 = w;
      f_v               = f_w;
      w                 = z;
      f_w               = f_z;
      z                 = u;
      f_z               = f_u;
    }
    else {
      if (u < z) {
        x_left          = u;
      }
      else {
        x_right         = u;
      }
      if ((f_u <= f_w) || (w == z)) {
        v               = w;
        f_v             = f_w;
        w               = u;
        f_w             = f_u;
      }
      else if ((f_u <= f_v) || (v == z) || (v == w)) {
        v               = u;
        f_v             = f_u;
      }
    }
    w_left              = z - x_left;
    w_right             = x_right - z;
    iter++;
  }
  output                = NumericVector::create(z, -f_z, 1, iter);
  return output;
}

// [[Rcpp::export]]
NumericMatrix bernard_des_one_stage_cpp(double alpha, double beta, double delta,
                                        NumericVector poss_nC,
                                        NumericVector poss_nE, List poss_x,
                                        List poss_B, List unique_B,
                                        int point_null, NumericVector Pi0,
                                        int point_alt, NumericVector Pi1,
                                        int summary) {
  int           nC,
                nE,
                counter      = 0,
                interrupt    = 0,
                nCmax        = max(poss_nC);
  double        e1,
                f1,
                power,
                typeI,
                pi_power     = pi_power_finder(point_alt, Pi1, delta),
                pi_typeI     = pi_typeI_finder(point_null, Pi0);
  NumericMatrix feasible_designs(10000000, 6);
  for (int n = 1; n <= poss_nC.length(); n++) {
    nC                       = poss_nC[n - 1];
    nE                       = poss_nE[n - 1];
    if ((summary == 1) && (nC%10 == 0)) {
      string str = std::to_string(nC);
      message_cpp("currently analysing designs with nC = ", str);
    }
    NumericVector unique_B_n = unique_B[nC + nCmax*(nE - 1) - 1];
    NumericMatrix prob_x_power(nC + 1, nE + 1),
                  prob_x_typeI(nC + 1, nE + 1),
                  poss_B_n   = poss_B[nC + nCmax*(nE - 1) - 1],
                  poss_x_n   = poss_x[nC + nCmax*(nE - 1) - 1],
    dbinom                   = dbinom_des_one_stage(pi_typeI, pi_power, delta,
                                                    nC, nE);
    for (int o = 0; o <= (nC + 1)*(nE + 1) - 1; o++) {
      prob_x_power(poss_x_n(o, 0), poss_x_n(o, 1)) =
        dbinom(1, poss_x_n(o, 0))*dbinom(3, poss_x_n(o, 1));
      prob_x_typeI(poss_x_n(o, 0), poss_x_n(o, 1)) =
        dbinom(0, poss_x_n(o, 0))*dbinom(2, poss_x_n(o, 1));
    }
    for (int ei = 0; ei <= unique_B_n.length() - 1; ei++) {
      interrupt++;
      if (interrupt % 1000 == 0) {
        Rcpp::checkUserInterrupt();
      }
      e1                     = unique_B_n[ei];
      f1                     = e1;
      power                  = 0;
      typeI                  = 0;
      for (int o = 0; o <= (nC + 1)*(nE + 1) - 1; o++) {
        if (poss_B_n(poss_x_n(o, 0), poss_x_n(o, 1)) >= e1) {
          power             += prob_x_power(poss_x_n(o, 0), poss_x_n(o, 1));
          typeI             += prob_x_typeI(poss_x_n(o, 0), poss_x_n(o, 1));
        }
      }
      if (power < 1 - beta) {
        break;
      }
      if (typeI <= alpha) {
        if (point_null == 1) {
          if (point_alt == 1) {
            feasible_designs(counter, _)   =
              NumericVector::create(nC, e1, pi_typeI, typeI, pi_power, power);
            counter++;
          }
          else {
            NumericVector min_power        =
              bernard_min_power(1, beta, delta, NumericVector::create(nC),
                                NumericVector::create(nE),
                                NumericVector::create(e1),
                                NumericVector::create(f1),
                                List::create(poss_x_n), List::create(poss_B_n),
                                Pi1, 1);
            if (min_power[1] < 1 - beta) {
              break;
            }
            feasible_designs(counter, _)   =
              NumericVector::create(nC, e1, pi_typeI, typeI, min_power[0],
                                    min_power[1]);
            counter++;
          }
        }
        else {
          NumericVector max_typeI          =
            bernard_max_typeI(1, alpha, NumericVector::create(nC),
                              NumericVector::create(nE),
                              NumericVector::create(e1),
                              NumericVector::create(f1), List::create(poss_x_n),
                              List::create(poss_B_n), Pi0, 1);
          if (max_typeI[1] <= alpha) {
            if (point_alt == 1) {
              feasible_designs(counter, _) =
                NumericVector::create(nC, e1, max_typeI[0], max_typeI[1],
                                      pi_power, power);
              counter++;
            }
            else {
              NumericVector min_power      =
                bernard_min_power(1, beta, delta, NumericVector::create(nC),
                                  NumericVector::create(nE),
                                  NumericVector::create(e1),
                                  NumericVector::create(f1),
                                  List::create(poss_x_n),
                                  List::create(poss_B_n), Pi1, 1);
              if (min_power[1] < 1 - beta) {
                break;
              }
              feasible_designs(counter, _) =
                NumericVector::create(nC, e1, max_typeI[0], max_typeI[1],
                                      min_power[0], min_power[1]);
              counter++;
            }
          }
        }
      }
    }
  }
  NumericMatrix output                     =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 5));
  return output;
}

// [[Rcpp::export]]
NumericMatrix bernard_des_two_stage_cpp(double alpha, double beta, double delta,
                                        NumericVector poss_nC,
                                        NumericVector poss_nE, List poss_x,
                                        List poss_B, List unique_B,
                                        int point_null, NumericVector Pi0,
                                        int point_alt, NumericVector Pi1,
                                        int equal, int efficacy, int futility,
                                        double pi_ess, int summary) {
  int           len_B1,
                len_B2,
                nC1,
                nE1,
                nC2,
                nE2,
                counter                 = 0,
                interrupt               = 0,
                nCmax                   = max(poss_nC);
  double        e1,
                e2,
                f1,
                f2,
                power1,
                power2,
                S1_ess0,
                S1_ess1,
                typeI1,
                typeI2,
                typeII1,
                pi_power                = pi_power_finder(point_alt, Pi1,
                                                          delta),
                pi_typeI                = pi_typeI_finder(point_null, Pi0);
  NumericMatrix dbinom_ess,
                dbinom1,
                dbinom2,
                feasible_designs(10000000, 11);
  for (int n1 = 0; n1 <= poss_nC.length() - 1; n1++) {
    nC1                                 = poss_nC[n1];
    nE1                                 = poss_nE[n1];
    Rcout << "The value is A = " << n1 << std::endl;
    if (((equal != 1) && (nC1 < nCmax)) ||
          ((equal == 1) && (nC1 <= 0.5*nCmax))) {
      if ((summary == 1) && (nC1%10 == 0)) {
        string str = std::to_string(nC1);
        message_cpp("currently analysing designs with nC1 = ", str);
      }
      Rcout << "The value is B = " << n1 << std::endl;
      NumericVector unique_B1           = unique_B[nC1 + nCmax*(nE1 - 1) - 1];
      len_B1                            = unique_B1.length();
      Rcout << "The value is C = " << len_B1 << std::endl;
      NumericMatrix prob_x1_ess0(nC1 + 1, nE1 + 1),
                    prob_x1_ess1(nC1 + 1, nE1 + 1),
                    prob_x1_power(nC1 + 1, nE1 + 1),
                    prob_x1_typeI(nC1 + 1, nE1 + 1),
                    poss_B1             = poss_B[nC1 + nCmax*(nE1 - 1) - 1],
                    poss_x1             = poss_x[nC1 + nCmax*(nE1 - 1) - 1],
                    dbinom1             =
                      dbinom_des_one_stage(pi_typeI, pi_power, delta, nC1, nE1),
                    dbinom_ess          =
                      dbinom_des_ess(dbinom1, pi_typeI, pi_power, delta, pi_ess,
                                     nC1, nE1);
      Rcout << "The value is D = " << len_B1 << std::endl;
      for (int o1 = 0; o1 <= (nC1 + 1)*(nE1 + 1) - 1; o1++) {
        prob_x1_ess0(poss_x1(o1, 0), poss_x1(o1, 1))  =
          dbinom_ess(0, poss_x1(o1, 0))*dbinom_ess(1, poss_x1(o1, 1));
        prob_x1_ess1(poss_x1(o1, 0), poss_x1(o1, 1))  =
          dbinom_ess(0, poss_x1(o1, 0))*dbinom_ess(2, poss_x1(o1, 1));
        prob_x1_power(poss_x1(o1, 0), poss_x1(o1, 1)) =
          dbinom1(1, poss_x1(o1, 0))*dbinom1(3, poss_x1(o1, 1));
        prob_x1_typeI(poss_x1(o1, 0), poss_x1(o1, 1)) =
          dbinom1(0, poss_x1(o1, 0))*dbinom1(2, poss_x1(o1, 1));
      }
      Rcout << "The value is E = " << len_B1 << std::endl;
      for (int fi1 = (futility == 1 ? 1 : 0);
           fi1 <= (futility == 1 ?
                     (efficacy == 1 ? len_B1 - 4 : len_B1 - 3) : 0); fi1++) {
        f1                              = unique_B1[fi1];
        Rcout << "The value is F = " << f1 << std::endl;
        typeII1                         = 0;
        for (int o1 = 0; o1 <= (nC1 + 1)*(nE1 + 1) - 1; o1++) {
          if (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) <= f1) {
            typeII1                    += prob_x1_power(poss_x1(o1, 0),
                                                        poss_x1(o1, 1));
          }
        }
        Rcout << "The value is G = " << f1 << std::endl;
        if (typeII1 > beta) {
          break;
        }
        Rcout << "The value is H = " << f1 << std::endl;
        for (int ei1 = (efficacy == 1 ? fi1 + 2 : len_B1 - 1);
             ei1 <= (efficacy == 1 ? len_B1 - 2 : len_B1 - 1); ei1++) {
          e1                            = unique_B1[ei1];
          power1                        = 0;
          typeI1                        = 0;
          for (int o1 = 0; o1 <= (nC1 + 1)*(nE1 + 1) - 1; o1++) {
            if (poss_B1(poss_x1(o1, 0), poss_B1(o1, 1)) >= e1) {
              power1                   += prob_x1_power(poss_x1(o1, 0),
                                                        poss_x1(o1, 1));
              typeI1                   += prob_x1_typeI(poss_x1(o1, 0),
                                                        poss_x1(o1, 1));
            }
          }
          Rcout << "The value is I = " << f1 << std::endl;
          if (typeI1 < alpha) {
            S1_ess0                     = 0;
            S1_ess1                     = 0;
            for (int o1 = 0; o1 <= (nC1 + 1)*(nE1 + 1) - 1; o1++) {
              if ((poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) >= e1) ||
                    (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) <= f1)) {
                S1_ess0                += prob_x1_ess0(poss_x1(o1, 0),
                                                       poss_x1(o1, 1));
                S1_ess1                += prob_x1_ess1(poss_x1(o1, 0),
                                                       poss_x1(o1, 1));
              }
            }
            Rcout << "The value is J = " << f1 << std::endl;
            for (int n2 = (equal == 1 ? n1 : 0);
                 n2 <= (equal == 1 ? n1 : poss_nC.length() - 1); n2++) {
              nC2                       = poss_nC[n2];
              Rcout << "The value is K = " << nC2 << std::endl;
              if (nC1 + nC2 <= nCmax) {
                nE2                     = poss_nE[n2];
                Rcout << "The value is nE2 = " << nE2 << std::endl;
                NumericVector unique_B2 = unique_B[nC1 + nC2 +
                  nCmax*(nE1 + nE2 - 1) - 1];
                len_B2                  = unique_B2.length();
                NumericMatrix prob_x_power(nC1 + nC2 + 1, nE1 + nE2 + 1),
                              prob_x_typeI(nC1 + nC2 + 1, nE1 + nE2 + 1),
                              poss_B2   = poss_B[nC1 + nC2 +
                                                   nCmax*(nE1 + nE2 - 1) - 1],
                              poss_x2   = poss_x[nC2 + nCmax*(nE2 - 1) - 1];
                dbinom2                 =
                  dbinom_des_two_stage(dbinom1, pi_typeI, pi_power, delta, nC1,
                                       nC2, nE1, nE2);
                Rcout << "The value is L = " << f1 << std::endl;
                for (int o1 = 0; o1 <= (nC1 + 1)*(nE1 + 1) - 1; o1++) {
                  if ((poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) > f1) &&
                      (poss_B1(poss_x1(o1, 0), poss_x1(o1, 1)) < e1)) {
                    for (int o2 = 0; o2 <= (nC2 + 1)*(nE2 + 1) - 1; o2++) {
                      prob_x_power(poss_x1(o1, 0) + poss_x2(o2, 0),
                                   poss_x1(o1, 1) + poss_x2(o2, 1)) +=
                        prob_x1_power(poss_x1(o1, 0), poss_x1(o1, 1))*
                        dbinom2(1, poss_x2(o2, 0))*dbinom2(3, poss_x2(o2, 1));
                      prob_x_typeI(poss_x1(o1, 0) + poss_x2(o2, 0),
                                   poss_x1(o1, 1) + poss_x2(o2, 1)) +=
                        prob_x1_typeI(poss_x1(o1, 0), poss_x1(o1, 1))*
                        dbinom2(0, poss_x2(o2, 0))*dbinom2(2, poss_x2(o2, 1));
                    }
                  }
                }
                Rcout << "The value is M = " << f1 << std::endl;
                for (int ei2 = 1; ei2 <= len_B2; ei2++) {
                  interrupt++;
                  if (interrupt % 1000 == 0) {
                    Rcpp::checkUserInterrupt();
                  }
                  e2                    = unique_B2[ei2];
                  f2                    = e2;
                  power2                = 0;
                  typeI2                = 0;
                  Rcout << "The value is f1 = " << f1 << std::endl;
                  Rcout << "The value is e1 = " << e1 << std::endl;
                  Rcout << "The value is e2 = " << e2 << std::endl;
                  Rcout << "The value is nC1 = " << nC1 << std::endl;
                  Rcout << "The value is nC2 = " << nC2 << std::endl;
                  Rcout << "The value is nE1 = " << nE1 << std::endl;
                  Rcout << "The value is nE2 = " << nE2 << std::endl;
                  for (int xC = 0; xC <= nC1 + nC2; xC++) {
                    for (int xE = nE1 + nE2; xE >= 0; xE--) {
                      if (poss_B2(xC, xE) >= e2) {
                        power2         += prob_x_power(xC, xE);
                        typeI2         += prob_x_typeI(xC, xE);
                      }
                      else {
                        break;
                      }
                    }
                  }
                  Rcout << "The value is power1 = " << power1 << std::endl;
                  Rcout << "The value is power2 = " << power2 << std::endl;
                  Rcout << "The value is typeI1 = " << typeI1 << std::endl;
                  Rcout << "The value is typeI2 = " << typeI2 << std::endl;
                  if (power1 + power2 < 1 - beta) {
                    break;
                  }
                  if (typeI1 + typeI2 <= alpha) {
                    if (point_null == 1) {
                      if (point_alt == 1) {
                        Rcout << "The value is P = " << f1 << std::endl;
                        feasible_designs(counter, _)   =
                          NumericVector::create(nC1, nC2, e1, e2, f1, pi_typeI,
                                                typeI1 + typeI2, pi_power,
                                                power1 + power2,
                                                nC1 + nE1 +
                                                  (1 - S1_ess0)*(nC2 + nE2),
                                                nC1 + nE1 +
                                                  (1 - S1_ess1)*(nC2 + nE2));
                        counter++;
                        Rcout << "The value is Q = " << f1 << std::endl;
                      }
                      else {
                        NumericVector min_power        =
                          bernard_min_power(2, beta, delta,
                                            NumericVector::create(nC1, nC2),
                                            NumericVector::create(nE1, nE2),
                                            NumericVector::create(e1, e2),
                                            NumericVector::create(f1, f2),
                                            List::create(poss_x1, poss_x2),
                                            List::create(poss_B1, poss_B2),
                                            Pi1, 1);
                        if (min_power[1] < 1 - beta) {
                          break;
                        }
                        feasible_designs(counter, _)   =
                          NumericVector::create(nC1, nC2, e1, e2, f1, pi_typeI,
                                                typeI1 + typeI2, min_power[0],
                                                min_power[1],
                                                nC1 + nE1 +
                                                  (1 - S1_ess0)*(nC2 + nE2),
                                                nC1 + nE1 +
                                                  (1 - S1_ess1)*(nC2 + nE2));
                        counter++;
                      }
                    }
                    else {
                      NumericVector max_typeI          =
                        bernard_max_typeI(2, alpha,
                                          NumericVector::create(nC1, nC2),
                                          NumericVector::create(nE1, nE2),
                                          NumericVector::create(e1, e2),
                                          NumericVector::create(f1, f2),
                                          List::create(poss_x1, poss_x2),
                                          List::create(poss_B1, poss_B2),
                                          Pi0, 1);
                      if (max_typeI[1] <= alpha) {
                        if (point_alt == 1) {
                          feasible_designs(counter, _) =
                            NumericVector::create(nC1, nC2, e1, e2, f1,
                                                  max_typeI[0], max_typeI[1],
                                                  pi_power, power1 + power2,
                                                  nC1 + nE1 +
                                                    (1 - S1_ess0)*(nC2 + nE2),
                                                  nC1 + nE1 +
                                                    (1 - S1_ess1)*(nC2 + nE2));
                          counter++;
                        }
                        else {
                          NumericVector min_power      =
                            bernard_min_power(2, beta, delta,
                                              NumericVector::create(nC1, nC2),
                                              NumericVector::create(nE1, nE2),
                                              NumericVector::create(e1, e2),
                                              NumericVector::create(f1, f2),
                                              List::create(poss_x1, poss_x2),
                                              List::create(poss_B1, poss_B2),
                                              Pi1, 1);
                          if (min_power[1] < 1 - beta) {
                            break;
                          }
                          feasible_designs(counter, _) =
                            NumericVector::create(nC1, nC2, e1, e2, f1,
                                                  max_typeI[0], max_typeI[1],
                                                  min_power[0], min_power[1],
                                                  nC1 + nE1 +
                                                    (1 - S1_ess0)*(nC2 + nE2),
                                                  nC1 + nE1 +
                                                    (1 - S1_ess1)*(nC2 + nE2));
                          counter++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  NumericMatrix output                  =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 10));
  return output;
}
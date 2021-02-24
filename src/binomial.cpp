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
NumericMatrix binomial_pmf_two_stage_cpp(NumericVector pi, NumericVector nC,
                                         NumericVector nE, double e1, double f1,
                                         double e2, NumericVector k) {
  int           counter                    = 0,
                sum_nC                     = sum(nC),
                sum_nE                     = sum(nE);
  NumericMatrix prob_x1(nC[0] + 1, nE[0] + 1),
                poss_x2(sum_nC + 1, sum_nE + 1),
                prob_x2(sum_nC + 1, sum_nE + 1),
                pmf((nC[0] + 1)*(nC[1] + 1)*(nE[0] + 1)*(nE[1] + 1), 8),
                dbinom                     = dbinom_two_stage(pi, nC, nE);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      prob_x1(xC1, xE1) = dbinom(0, xC1)*dbinom(2, xE1);
      if (((xE1 - xC1 <= f1) || (xE1 - xC1 >= e1)) && (k[0] == 1)) {
        pmf(counter, _)                    =
          NumericVector::create(xC1, xE1, nC[0], nE[0], xE1 - xC1,
                                (xE1 - xC1 >= e1), 1, prob_x1(xC1, xE1));
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
          pmf(counter, _)                  =
            NumericVector::create(xC, xE, sum_nC, sum_nE, xE - xC,
                                  (xE - xC >= e2), 2, prob_x2(xC, xE));
          counter++;
        }
      }
    }
  }
  NumericMatrix output = pmf(Range(0, counter - 1), Range(0, 7));
  return output;
}

// [[Rcpp::export]]
double binomial_power_one_stage(NumericVector pi, int nC, int nE, double e,
                                NumericMatrix poss_x, NumericVector poss_y) {
  double        power  = 0;
  NumericMatrix dbinom = dbinom_one_stage(pi, nC, nE);
  for (int o = 0; o <= (nC + 1)*(nE + 1) - 1; o++) {
    if (poss_y[o] >= e) {
      power           += dbinom(0, poss_x(o, 0))*dbinom(1, poss_x(o, 1));
    }
  }
  return power;
}

// [[Rcpp::export]]
double binomial_power_two_stage(NumericVector pi, NumericVector nC,
                                NumericVector nE, NumericVector e,
                                NumericVector f, List poss_x, List poss_y) {
  double        power                 = 0;
  NumericVector prob_x1(nC[0] + nE[0] + 1),
                poss_y1               = poss_y[0],
                poss_y2               = poss_y[1];
  NumericMatrix dbinom                = dbinom_two_stage(pi, nC, nE),
                poss_x1               = poss_x[0],
                poss_x2               = poss_x[1];
  if ((f[0] >= -nC[0]) && (e[0] <= nE[0])) {
    for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
      if (poss_y1[o1] >= e[0]) {
        power                        += dbinom(0, poss_x1(o1, 0))*
                                          dbinom(2, poss_x1(o1, 1));
      }
      else if (poss_y1[o1] > f[0]) {
        prob_x1[poss_y1[o1] + nC[0]] += dbinom(0, poss_x1(o1, 0))*
                                          dbinom(2, poss_x1(o1, 1));
      }
    }
  }
  else if (f[0] >= -nC[0]) {
    for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
      if (poss_y1[o1] > f[0]) {
        prob_x1[poss_y1[o1] + nC[0]] += dbinom(0, poss_x1(o1, 0))*
                                          dbinom(2, poss_x1(o1, 1));
      }
    }
  }
  else {
    for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
      if (poss_y1[o1] >= e[0]) {
        power                        += dbinom(0, poss_x1(o1, 0))*
                                          dbinom(2, poss_x1(o1, 1));
      }
      else {
        prob_x1[poss_y1[o1] + nC[0]] += dbinom(0, poss_x1(o1, 0))*
                                          dbinom(2, poss_x1(o1, 1));
      }
    }
  }
  for (int y = f[0] + 1; y <= e[0] - 1; y++) {
    for (int o2 = 0; o2 <= (nC[1] + 1)*(nE[1] + 1) - 1; o2++) {
      if (y + poss_y2[o2] >= e[1]) {
        power                        += prob_x1[y + nC[0]]*
                                          dbinom(1, poss_x2(o2, 0))*
                                          dbinom(3, poss_x2(o2, 1));
      }
    }
  }
  return power;
}

// [[Rcpp::export]]
NumericMatrix binomial_terminal_two_stage_cpp(NumericVector nC,
                                              NumericVector nE, double e1,
                                              double f1, double e2,
                                              NumericVector k) {
  int           counter                     = 0,
                sum_nC                      = sum(nC),
                sum_nE                      = sum(nE);
  NumericMatrix x2_mat(sum_nC + 1, sum_nE + 1),
  terminal((nC[0] + 1)*(nC[1] + 1)*(nE[0] + 1)*(nE[1] + 1), 7);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      if (k[0] == 1) {
        if ((xE1 - xC1 <= f1) || (xE1 - xC1 >= e1)) {
          terminal(counter, _)                =
            NumericVector::create(xC1, xE1, nC[0], nE[0], xE1 - xC1,
                                  (xE1 - xC1 >= e1) + 1, 1);
        }
        else {
          terminal(counter, _)                =
            NumericVector::create(xC1, xE1, nC[0], nE[0], xE1 - xC1, 3, 1);
        }
        counter++;
      }
      if ((xE1 - xC1 > f1) && (xE1 - xC1 < e1) &&
            ((k[0] == 2) || (k[k.length() - 1] == 2))) {
        for (int xC2 = 0; xC2 <= nC[1]; xC2++) {
          for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
            if (x2_mat(xC1 + xC2, xE1 + xE2) == 0) {
              terminal(counter, _)          =
                NumericVector::create(xC1 + xC2, xE1 + xE2, sum_nC, sum_nE,
                                      xE1 + xE2 - xC1 - xC2,
                                      (xE1 + xE2 - xC1 - xC2 >= e2) + 1, 2);
              x2_mat(xC1 + xC2, xE1 + xE2) += 1;
              counter++;
            }
          }
        }
      }
    }
  }
  NumericMatrix output                      = terminal(Range(0, counter - 1),
                                                       Range(0, 6));
  return output;
}

// [[Rcpp::export]]
NumericVector binomial_max_typeI(int J, double alpha, NumericVector nC,
                                 NumericVector nE, NumericVector e,
                                 NumericVector f, List poss_x, List poss_y,
                                 NumericVector Pi0, int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = Pi0[0],
                x_right = Pi0[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      -binomial_power_one_stage(NumericVector::create(v, v), nC[0], nE[0], e[0],
                                poss_x[0], poss_y[0]);
  }
  else {
    f_v                 =
      -binomial_power_two_stage(NumericVector::create(v, v), nC, nE, e, f,
                                poss_x, poss_y);
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
      -binomial_power_one_stage(NumericVector::create(z, z), nC[0], nE[0], e[0],
                                poss_x[0], poss_y[0]);
  }
  else {
    f_z                 =
      -binomial_power_two_stage(NumericVector::create(z, z), nC, nE, e, f,
                                poss_x, poss_y);
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
        -binomial_power_one_stage(NumericVector::create(u, u), nC[0], nE[0],
                                  e[0], poss_x[0], poss_y[0]);
    }
    else {
      f_u               =
        -binomial_power_two_stage(NumericVector::create(u, u), nC, nE, e, f,
                                  poss_x, poss_y);
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
NumericVector binomial_min_power(int J, double beta, double delta,
                                 NumericVector nC, NumericVector nE,
                                 NumericVector e, NumericVector f, List poss_x,
                                 List poss_y, NumericVector pi_alt, int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = pi_alt[0],
                x_right = pi_alt[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      binomial_power_one_stage(NumericVector::create(v, v + delta), nC[0],
                               nE[0], e[0], poss_x[0], poss_y[0]);
  }
  else {
    f_v                 =
      binomial_power_two_stage(NumericVector::create(v, v + delta), nC, nE, e,
                               f, poss_x, poss_y);
  }
  NumericVector output(4);
  if ((f_v < 1 - beta) && (check == 1)) {
    output              = NumericVector::create(v, f_v, 2, iter);
    return output;
  }
  double f_z,
         z              = pi_power_finder(0, pi_alt, delta);
  if (J == 1) {
    f_z                 =
      binomial_power_one_stage(NumericVector::create(z, z + delta), nC[0],
                               nE[0], e[0], poss_x[0], poss_y[0]);
  }
  else {
    f_z                 =
      binomial_power_two_stage(NumericVector::create(z, z + delta), nC, nE, e,
                               f, poss_x, poss_y);
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
        binomial_power_one_stage(NumericVector::create(u, u + delta), nC[0],
                                 nE[0], e[0], poss_x[0], poss_y[0]);
    }
    else {
      f_u               =
        binomial_power_two_stage(NumericVector::create(u, u + delta), nC, nE, e,
                                 f, poss_x, poss_y);
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
double binomial_ess_two_stage(NumericVector pi, NumericVector nC,
                              NumericVector nE, int e1, int f1) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, nC[0], nE[0]);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      if ((xE1 - xC1 <= f1) || (xE1 - xC1 >= e1)) {
        S1             += dbinom1(0, xC1)*dbinom1(1, xE1);
      }
    }
  }
  double ess            = nC[0] + nE[0] + (1 - S1)*(nC[1] + nE[1]);
  return ess;
}

// [[Rcpp::export]]
double binomial_des_ess_two_stage(NumericVector pi, NumericVector nC,
                                  NumericVector nE, int e1, int f1,
                                  NumericMatrix poss_x1,
                                  NumericVector poss_y1) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, nC[0], nE[0]);
  for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
    if ((poss_y1[o1] <= f1) || (poss_y1[o1] >= e1)) {
      S1               += dbinom1(0, poss_x1(o1, 0))*dbinom1(1, poss_x1(o1, 1));
    }
  }
  double ess            = nC[0] + nE[0] + (1 - S1)*(nC[1] + nE[1]);
  return ess;
}

// [[Rcpp::export]]
NumericVector binomial_max_ess_1d_two_stage(NumericVector nC, NumericVector nE,
                                            int e1, int f1,
                                            NumericMatrix poss_x1,
                                            NumericVector poss_y1) {
  int           iter    = 0;
  double        f_u,
                midpoint,
                u,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = 0,
                x_right = 1,
                v       = x_left + golden*(x_right - x_left),
                f_v     =
                  -binomial_des_ess_two_stage(NumericVector::create(v, v), nC,
                                              nE, e1, f1, poss_x1, poss_y1),
                z       = 0.5,
                f_z     =
                  -binomial_des_ess_two_stage(NumericVector::create(z, z), nC,
                                              nE, e1, f1, poss_x1, poss_y1),
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
        -binomial_des_ess_two_stage(NumericVector::create(u, u), nC, nE, e1, f1,
                                    poss_x1, poss_y1);
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
NumericMatrix binomial_des_one_stage_cpp(double alpha, double beta,
                                         double delta, NumericVector poss_nC,
                                         NumericVector poss_nE, List poss_x,
                                         List poss_y, int point_null,
                                         NumericVector Pi0, int point_alt,
                                         NumericVector pi_alt, int summary) {
  int           f1,
                nC,
                nE,
                counter                    = 0,
                interrupt                  = 0,
                nCmax                      = max(poss_nC);
  double        power,
                typeI,
                pi_power                   = pi_power_finder(point_alt, pi_alt,
                                                             delta),
                pi_typeI                   = pi_typeI_finder(point_null,
                                                             Pi0);
  NumericMatrix dbinom,
                feasible_designs(10000000, 6);
  for (int n = 0; n <= poss_nC.length() - 1; n++) {
    nC                                     = poss_nC[n];
    nE                                     = poss_nE[n];
    if ((summary == 1) && (nC%10 == 0)) {
      Rprintf("currently analysing designs with nC = %i \n", nC);
    }
    NumericVector prob_y_power(nC + nE + 1),
                  prob_y_typeI(nC + nE + 1);
    NumericVector poss_y_n                  = poss_y[nC + nCmax*(nE - 1) - 1];
    NumericMatrix poss_x_n                  = poss_x[nC + nCmax*(nE - 1) - 1];
    dbinom                                  =
      dbinom_des_one_stage(pi_typeI, pi_power, delta, nC, nE);
    for (int o = 0; o <= (nC + 1)*(nE + 1) - 1; o++) {
      prob_y_power[poss_y_n[o] + nC]       += dbinom(1, poss_x_n(o, 0))*
                                                dbinom(3, poss_x_n(o, 1));
      prob_y_typeI[poss_y_n[o] + nC]       += dbinom(0, poss_x_n(o, 0))*
                                                dbinom(2, poss_x_n(o, 1));
    }
    for (int e1 = -nC + 1; e1 <= nE; e1++) {
      interrupt++;
      if (interrupt % 1000 == 0) {
        Rcpp::checkUserInterrupt();
      }
      f1                                   = e1 - 1;
      power                                = 0;
      typeI                                = 0;
      for (int y = e1; y <= nE; y++) {
        power                             += prob_y_power[y + nC];
        typeI                             += prob_y_typeI[y + nC];
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
              binomial_min_power(1, beta, delta, NumericVector::create(nC),
                                 NumericVector::create(nE),
                                 NumericVector::create(e1),
                                 NumericVector::create(f1),
                                 List::create(poss_x_n), List::create(poss_y_n),
                                 pi_alt, 1);
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
            binomial_max_typeI(1, alpha, NumericVector::create(nC),
                               NumericVector::create(nE),
                               NumericVector::create(e1),
                               NumericVector::create(f1),
                               List::create(poss_x_n), List::create(poss_y_n),
                               Pi0, 1);
          if (max_typeI[1] <= alpha) {
            if (point_alt == 1) {
              feasible_designs(counter, _) =
                NumericVector::create(nC, e1, max_typeI[0], max_typeI[1],
                                      pi_power, power);
              counter++;
            }
            else {
              NumericVector min_power      =
                binomial_min_power(1, beta, delta, NumericVector::create(nC),
                                   NumericVector::create(nE),
                                   NumericVector::create(e1),
                                   NumericVector::create(f1),
                                   List::create(poss_x_n),
                                   List::create(poss_y_n), pi_alt, 1);
              if (min_power[1] < 1 - beta) {
                break;
              }
            feasible_designs(counter, _)   =
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
NumericMatrix binomial_des_two_stage_cpp(double alpha, double beta,
                                         double delta, NumericVector poss_nC,
                                         NumericVector poss_nE, List poss_x,
                                         List poss_y, int point_null,
                                         NumericVector Pi0, int point_alt,
                                         NumericVector pi_alt, int equal,
                                         int efficacy, int futility,
                                         double pi_ess, int summary) {
  int           f2,
                n1C,
                n1E,
                n2C,
                n2E,
                counter                   = 0,
                interrupt                 = 0,
                nCmax                     = max(poss_nC);
  double        power1,
                power2,
                S2_ess0,
                S2_ess1,
                typeI1,
                typeI2,
                typeII1,
                pi_power                  =
                  pi_power_finder(point_alt, pi_alt, delta),
                pi_typeI                  =
                  pi_typeI_finder(point_null, Pi0);
  NumericMatrix dbinom_ess,
                dbinom1,
                dbinom2,
                feasible_designs(10000000, 11);
  for (int n1 = 0; n1 <= poss_nC.length() - 1; n1++) {
    n1C                                   = poss_nC[n1];
    n1E                                   = poss_nE[n1];
    if (((equal != 1) && (n1C < nCmax)) ||
          ((equal == 1) && (n1C <= 0.5*nCmax))) {
      if ((summary == 1) && (n1C%10 == 0)) {
        Rprintf("currently analysing designs with n1C = %i \n", n1C);
      }
      NumericVector prob_y1_ess0(n1C + n1E + 1),
                    prob_y1_ess1(n1C + n1E + 1),
                    prob_y1_power(n1C + n1E + 1),
                    prob_y1_typeI(n1C + n1E + 1),
                    poss_y1               = poss_y[n1C + nCmax*(n1E - 1) - 1];
      NumericMatrix poss_x1               = poss_x[n1C + nCmax*(n1E - 1) - 1];
      dbinom1                             =
        dbinom_des_one_stage(pi_typeI, pi_power, delta, n1C, n1E);
      dbinom_ess                          =
        dbinom_des_ess(dbinom1, pi_typeI, pi_power, delta, pi_ess, n1C, n1E);
      for (int o1 = 0; o1 <= (n1C + 1)*(n1E + 1) - 1; o1++) {
        prob_y1_ess0[poss_y1[o1] + n1C]  += dbinom_ess(0, poss_x1(o1, 0))*
                                              dbinom_ess(1, poss_x1(o1, 1));
        prob_y1_ess1[poss_y1[o1] + n1C]  += dbinom_ess(0, poss_x1(o1, 0))*
                                              dbinom_ess(2, poss_x1(o1, 1));
        prob_y1_power[poss_y1[o1] + n1C] += dbinom1(1, poss_x1(o1, 0))*
                                              dbinom1(3, poss_x1(o1, 1));
        prob_y1_typeI[poss_y1[o1] + n1C] += dbinom1(0, poss_x1(o1, 0))*
                                              dbinom1(2, poss_x1(o1, 1));
      }
      for (int f1 = (futility == 1 ? -n1C : -n1C - 1);
           f1 <= (futility == 1 ?
                    (efficacy == 1 ? n1E - 2 : n1E - 1) : -n1C - 1); f1++) {
        typeII1                           = 0;
        for (int y1 = -n1C; y1 <= f1; y1++) {
          typeII1                        += prob_y1_power[y1 + n1C];
        }
        if (typeII1 > beta) {
          break;
        }
        for (int e1 = (efficacy == 1 ? f1 + 2 : n1E + 1);
             e1 <= (efficacy == 1 ? n1E : n1E + 1); e1++) {
          power1                          = 0;
          typeI1                          = 0;
          for (int y1 = e1; y1 <= n1E; y1++) {
            power1                       += prob_y1_power[y1 + n1C];
            typeI1                       += prob_y1_typeI[y1 + n1C];
          }
          if (typeI1 < alpha) {
            S2_ess0                       = 0;
            S2_ess1                       = 0;
            for (int y1 = f1 + 1; y1 <= e1 - 1; y1++) {
              S2_ess0                    += prob_y1_ess0[y1 + n1C];
              S2_ess1                    += prob_y1_ess1[y1 + n1C];
            }
            for (int n2 = (equal == 1 ? n1 : 0);
                 n2 <= (equal == 1 ? n1 : poss_nC.length() - 1); n2++) {
              n2C                         = poss_nC[n2];
              if (n1C + n2C <= nCmax) {
                n2E                       = poss_nE[n2];
                NumericVector prob_y_power(n1C + n1E + n2C + n2E + 1),
                              prob_y_typeI(n1C + n1E + n2C + n2E + 1),
                              poss_y2     = poss_y[n2C + nCmax*(n2E - 1) - 1];
                NumericMatrix poss_x2     = poss_x[n2C + nCmax*(n2E - 1) - 1];
                dbinom2                   =
                  dbinom_des_two_stage(dbinom1, pi_typeI, pi_power, delta, n1C,
                                       n2C, n1E, n2E);
                for (int y1 = f1 + 1; y1 <= e1 - 1; y1++) {
                  for (int o2 = 0; o2 <= (n2C + 1)*(n2E + 1) - 1; o2++) {
                    prob_y_power[y1 + poss_y2[o2] + n1C + n2C] +=
                      prob_y1_power[y1 + n1C]*dbinom2(1, poss_x2(o2, 0))*
                      dbinom2(3, poss_x2(o2, 1));
                    prob_y_typeI[y1 + poss_y2[o2] + n1C + n2C] +=
                      prob_y1_typeI[y1 + n1C]*dbinom2(0, poss_x2(o2, 0))*
                      dbinom2(2, poss_x2(o2, 1));
                  }
                }
                for (int e2 = f1 + 2 - n2C; e2 <= e1 - 1 + n2E; e2++) {
                  interrupt++;
                  if (interrupt % 1000 == 0) {
                    Rcpp::checkUserInterrupt();
                  }
                  f2                      = e2 - 1;
                  power2                  = 0;
                  typeI2                  = 0;
                  for (int y_index = n1C + n2C + e2;
                       y_index <= n1C + n1E + n2C + n2E; y_index++) {
                    power2               += prob_y_power[y_index];
                    typeI2               += prob_y_typeI[y_index];
                  }
                  if (power1 + power2 < 1 - beta) {
                    break;
                  }
                  if (typeI1 + typeI2 <= alpha) {
                    if (point_null == 1) {
                      if (point_alt == 1) {
                        feasible_designs(counter, _)   =
                          NumericVector::create(n1C, n2C, e1, e2, f1, pi_typeI,
                                                typeI1 + typeI2, pi_power,
                                                power1 + power2,
                                                n1C + n1E + S2_ess0*(n2C + n2E),
                                                n1C + n1E +
                                                  S2_ess1*(n2C + n2E));
                        counter++;
                      }
                      else {
                        NumericVector min_power        =
                          binomial_min_power(2, beta, delta,
                                             NumericVector::create(n1C, n2C),
                                             NumericVector::create(n1E, n2E),
                                             NumericVector::create(e1, e2),
                                             NumericVector::create(f1, f2),
                                             List::create(poss_x1, poss_x2),
                                             List::create(poss_y1, poss_y2),
                                             pi_alt, 1);
                        if (min_power[1] < 1 - beta) {
                          break;
                        }
                        feasible_designs(counter, _)   =
                          NumericVector::create(n1C, n2C, e1, e2, f1, pi_typeI,
                                                typeI1 + typeI2, min_power[0],
                                                min_power[1],
                                                n1C + n1E + S2_ess0*(n2C + n2E),
                                                n1C + n1E +
                                                  S2_ess1*(n2C + n2E));
                        counter++;
                      }
                    }
                    else {
                      NumericVector max_typeI          =
                        binomial_max_typeI(2, alpha,
                                           NumericVector::create(n1C, n2C),
                                           NumericVector::create(n1E, n2E),
                                           NumericVector::create(e1, e2),
                                           NumericVector::create(f1, f2),
                                           List::create(poss_x1, poss_x2),
                                           List::create(poss_y1, poss_y2),
                                           Pi0, 1);
                      if (max_typeI[1] <= alpha) {
                        if (point_alt == 1) {
                          feasible_designs(counter, _) =
                            NumericVector::create(n1C, n2C, e1, e2, f1,
                                                  max_typeI[0], max_typeI[1],
                                                  pi_power, power1 + power2,
                                                  n1C + n1E +
                                                    S2_ess0*(n2C + n2E),
                                                  n1C + n1E + 
                                                    S2_ess1*(n2C + n2E));
                          counter++;
                        }
                        else {
                          NumericVector min_power      =
                            binomial_min_power(2, beta, delta,
                                               NumericVector::create(n1C, n2C),
                                               NumericVector::create(n1E, n2E),
                                               NumericVector::create(e1, e2),
                                               NumericVector::create(f1, f2),
                                               List::create(poss_x1, poss_x2),
                                               List::create(poss_y1, poss_y2),
                                               pi_alt, 1);
                          if (min_power[1] < 1 - beta) {
                            break;
                          }
                          feasible_designs(counter, _) =
                            NumericVector::create(n1C, n2C, e1, e2, f1,
                                                  max_typeI[0], max_typeI[1],
                                                  min_power[0], min_power[1],
                                                  n1C + n1E +
                                                    S2_ess0*(n2C + n2E),
                                                  n1C + n1E +
                                                    S2_ess1*(n2C + n2E));
                          counter++;
                        }
                      }
                    }
                  }
                }
              }
              else {
                break;
              }
            }
          }
        }
      }
    }
  }
  NumericMatrix output                    =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 10));
  return output;
}

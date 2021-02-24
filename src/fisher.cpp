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
NumericMatrix fisher_pmf_two_stage_cpp(NumericVector pi, NumericVector nC,
                                       NumericVector nE, NumericVector e1,
                                       NumericVector f1, NumericMatrix e2,
                                       NumericVector k) {
  
  int           counter     = 0,
                sum_nC      = sum(nC),
                sum_nE      = sum(nE);
  NumericMatrix prob_x1(nC[0] + 1, nE[0] + 1),
                poss_x2(sum_nC + 1, sum_nE + 1),
                prob_x2(sum_nC + 1, sum_nE + 1),
                pmf((nC[0] + 1)*(nC[1] + 1)*(nE[0] + 1)*(nE[1] + 1), 12),
                dbinom      = dbinom_two_stage(pi, nC, nE);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      prob_x1(xC1, xE1)     = dbinom(0, xC1)*dbinom(2, xE1);
      if (((xE1 - xC1 <= f1[xC1 + xE1]) || (xE1 - xC1 >= e1[xC1 + xE1])) &&
          (k[0] == 1)) {
        pmf(counter, _)     =
          NumericVector::create(xC1, xE1, -1, -1, nC[0], nE[0], xC1 + xE1, -1,
                                xE1 - xC1, (xE1 - xC1 >= e1[xC1 + xE1]), 1,
                                prob_x1(xC1, xE1));
        counter++;
      }
      else if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
        for (int xC2 = 0; xC2 <= nC[1]; xC2++) {
          for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
            pmf(counter, _) =
              NumericVector::create(xC1, xE1, xC2, xE2, sum_nC, sum_nE,
                                    xC1 + xE1, xC2 + xE2, xE1 + xE2 - xC1 - xC2,
                                    (xE1 + xE2 - xC1 - xC2 >=
                                      e2(xC1 + xE1, xC2 + xE2)), 2,
                                      prob_x1(xC1, xE1)*dbinom(1, xC2)*
                                        dbinom(3, xE2));
            counter++;
          }
        }
      }
    }
  }
  NumericMatrix output      = pmf(Range(0, counter - 1), Range(0, 11));
  return output;
}

// [[Rcpp::export]]
double fisher_power_one_stage(NumericVector pi, int nC, int nE, NumericVector e,
                              NumericMatrix poss_x, NumericVector poss_y,
                              NumericVector poss_z) {
  double        power  = 0;
  NumericMatrix dbinom = dbinom_one_stage(pi, nC, nE);
  for (int o = 0; o <= (nC + 1)*(nE + 1) - 1; o++) {
    if (poss_y[o] >= e[poss_z[o]]) {
      power           += dbinom(0, poss_x(o, 0))*dbinom(1, poss_x(o, 1));
    }
  }
  return power;
}

// [[Rcpp::export]]
double fisher_power_two_stage(NumericVector pi, NumericVector nC,
                              NumericVector nE, NumericVector e1,
                              NumericVector f1, NumericMatrix e2,
                              List poss_x, List poss_y, List poss_z) {
  double        prob_xo1,
                power   = 0;
  NumericVector poss_y1 = poss_y[0],
                poss_y2 = poss_y[1],
                poss_z1 = poss_z[0],
                poss_z2 = poss_z[1];
  NumericMatrix poss_x1 = poss_x[0],
                poss_x2 = poss_x[1],
                dbinom  = dbinom_two_stage(pi, nC, nE);
  for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
    if (poss_y1[o1] >= e1[poss_z1[o1]]) {
      power            += dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
    }
    else if (poss_y1[o1] > f1[poss_z1[o1]]) {
      prob_xo1          = dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
      for (int o2 = 0; o2 <= (nC[1] + 1)*(nE[1] + 1) - 1; o2++) {
        if (poss_y1[o1] + poss_y2[o2] >= e2(poss_z1[o1], poss_z2[o2])) {
          power        += prob_xo1*dbinom(1, poss_x2(o2, 0))*
                            dbinom(3, poss_x2(o2, 1));
        }
      }
    }
  }
  return power;
}

// [[Rcpp::export]]
NumericMatrix fisher_terminal_two_stage_cpp(NumericVector nC,
                                            NumericVector nE, NumericVector e1,
                                            NumericMatrix e2, NumericVector f1,
                                            NumericVector k) {
  
  int           counter          = 0,
                sum_nC           = sum(nC),
                sum_nE           = sum(nE);
  NumericMatrix terminal(100000000, 11);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xE1++) {
      if (k[0] == 1) {
        if ((xE1 - xC1 <= f1[xC1 + xE1]) || (xE1 - xC1 >= e1[xC1 + xE1])) {
          terminal(counter, _)     =
            NumericVector::create(xC1, xE1, -1, -1, nC[0], nE[0], xC1 + xE1, -1,
                                  xE1 - xC1, (xE1 - xC1 >= e1[xC1 + xE1]) + 1,
                                  1);
        }
        else {
          terminal(counter, _)     =
            NumericVector::create(xC1, xE1, -1, -1, nC[0], nE[0], xC1 + xE1, -1,
                                  xE1 - xC1, 3, 1);
        }
        counter++;
      }
      if ((xE1 - xC1 > f1[xC1 + xE1]) && (xE1 - xC1 < e1[xC1 + xE1]) &&
            ((k[0] == 2) || (k[k.length() - 1] == 2))) {
        for (int xC2 = 0; xC2 <= nC[1]; xC2++) {
          for (int xE2 = 0; xE2 <= nE[1]; xE2++) {
            terminal(counter, _) =
              NumericVector::create(xC1, xE1, xC2, xE2, sum_nC, sum_nE,
                                    xC1 + xE1, xC2 + xE2, xE1 + xE2 - xC1 - xC2,
                                    (xE1 + xE2 - xC1 - xC2 >=
                                      e2(xC1 + xE1, xC2 + xE2)) + 1, 2);
            counter++;
          }
        }
      }
    }
  }
  NumericMatrix output           = terminal(Range(0, counter - 1),
                                            Range(0, 10));
  return output;
}

// [[Rcpp::export]]
NumericVector fisher_max_typeI(int J, double alpha, NumericVector nC,
                               NumericVector nE, NumericVector e1,
                               NumericVector f1, NumericMatrix e2, List poss_x,
                               List poss_y, List poss_z, NumericVector pi_null,
                               int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = pi_null[0],
                x_right = pi_null[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 = -fisher_power_one_stage(NumericVector::create(v, v),
                                                  nC[0], nE[0], e1, poss_x[0],
                                                  poss_y[0], poss_z[0]);
  }
  else {
    f_v                 = -fisher_power_two_stage(NumericVector::create(v, v),
                                                  nC, nE, e1, f1, e2, poss_x,
                                                  poss_y, poss_z);
  }
  NumericVector output(4);
  if ((-f_v > alpha) && (check == 1)) {
    output              = NumericVector::create(v, -f_v, 2, iter);
    return output;
  }
  double f_z,
         z              = pi_typeI_finder(0, pi_null);
  if (J == 1) {
    f_z                 = -fisher_power_one_stage(NumericVector::create(z, z),
                                                  nC[0], nE[0], e1, poss_x[0],
                                                  poss_y[0], poss_z[0]);
  }
  else {
    f_z                 = -fisher_power_two_stage(NumericVector::create(z, z),
                                                  nC, nE, e1, f1, e2, poss_x,
                                                  poss_y, poss_z);
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
      if ((u - x_left < two_tol) | (x_right - u < two_tol)) {
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
      f_u               = -fisher_power_one_stage(NumericVector::create(u, u),
                                                  nC[0], nE[0], e1, poss_x[0],
                                                  poss_y[0], poss_z[0]);
    }
    else {
      f_u               = -fisher_power_two_stage(NumericVector::create(u, u),
                                                  nC, nE, e1, f1, e2, poss_x,
                                                  poss_y, poss_z);
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
      if (f_u <= f_w || w == z) {
        v               = w;
        f_v             = f_w;
        w               = u;
        f_w             = f_u;
      }
      else if (f_u <= f_v || v == z || v == w) {
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
NumericVector fisher_min_power(int J, double beta, double delta,
                               NumericVector nC, NumericVector nE,
                               NumericVector e1, NumericVector f1,
                               NumericMatrix e2, List poss_x, List poss_y,
                               List poss_z, NumericVector pi_alt, int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = pi_alt[0],
                x_right = pi_alt[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      fisher_power_one_stage(NumericVector::create(v, v + delta), nC[0], nE[0],
                             e1, poss_x[0], poss_y[0], poss_z[0]);
  }
  else {
    f_v                 =
      fisher_power_two_stage(NumericVector::create(v, v + delta), nC, nE, e1,
                             f1, e2, poss_x, poss_y, poss_z);
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
      fisher_power_one_stage(NumericVector::create(z, z + delta), nC[0], nE[0],
                             e1, poss_x[0], poss_y[0], poss_z[0]);
  }
  else {
    f_z                 =
      fisher_power_two_stage(NumericVector::create(z, z + delta), nC, nE, e1,
                             f1, e2, poss_x, poss_y, poss_z);
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
    if (abs(p) < abs(0.5*q*r) && p < q*w_left && p < q*w_right) {
      d                 = p/q;
      u                 = z + d;
      if ((u - x_left < two_tol) | (x_right - u < two_tol)) {
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
      f_u                 =
        fisher_power_one_stage(NumericVector::create(u, u + delta), nC[0],
                               nE[0], e1, poss_x[0], poss_y[0], poss_z[0]);
    }
    else {
      f_u                 =
        fisher_power_two_stage(NumericVector::create(u, u + delta), nC, nE, e1,
                               f1, e2, poss_x, poss_y, poss_z);
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
double fisher_ess_two_stage(NumericVector pi, NumericVector nC,
                            NumericVector nE, NumericVector e1,
                            NumericVector f1) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, nC[0], nE[0]);
  for (int xC1 = 0; xC1 <= nC[0]; xC1++) {
    for (int xE1 = 0; xE1 <= nE[0]; xC1++) {
      if ((xE1 - xC1 >= e1[xC1 + xE1]) || (xE1 - xC1 <= f1[xC1 + xE1])) {
        S1             += dbinom1(0, xC1)*dbinom1(1, xE1);
      }
    }
  }
  double ess            = nC[0] + nE[0] + (1 - S1)*(nC[1] + nE[1]);
  return ess;
}

// [[Rcpp::export]]
double fisher_des_ess_two_stage(NumericVector pi, NumericVector nC,
                                NumericVector nE, NumericVector e1,
                                NumericVector f1, NumericMatrix poss_x1,
                                NumericVector poss_y1, NumericVector poss_z1) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, nC[0], nE[0]);
  for (int o1 = 0; o1 <= (nC[0] + 1)*(nE[0] + 1) - 1; o1++) {
    if ((poss_y1[o1] >= e1[poss_z1[o1]]) ||
          (poss_x1(o1, 1) <= f1[poss_z1[o1]])) {
      S1               += dbinom1(0, poss_x1(o1, 0))*dbinom1(1, poss_x1(o1, 1));
    }
  }
  double ess            = nC[0] + nE[0] + (1 - S1)*(nC[1] + nE[1]);
  return ess;
}

// [[Rcpp::export]]
NumericVector fisher_max_ess_1d_two_stage(NumericVector nC, NumericVector nE,
                                          NumericVector e1, NumericVector f1,
                                          NumericMatrix poss_x1,
                                          NumericVector poss_y1,
                                          NumericVector poss_z1) {
  int           iter    = 0;
  double        f_u,
                midpoint,
                u,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = 0,
                x_right = 1,
                v       = x_left + golden*(x_right - x_left),
                f_v     =
                  -fisher_des_ess_two_stage(NumericVector::create(v, v), nC, nE,
                                            e1, f1, poss_x1, poss_y1, poss_z1),
                z       = 0.5,
                f_z     =
                  -fisher_des_ess_two_stage(NumericVector::create(z, z), nC, nE,
                                            e1, f1, poss_x1, poss_y1, poss_z1),
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
      -fisher_des_ess_two_stage(NumericVector::create(u, u), nC, nE, e1, f1,
                                poss_x1, poss_y1, poss_z1);
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
NumericMatrix fisher_des_one_stage_cpp(double alpha, double beta, double delta,
                                       NumericVector poss_nC,
                                       NumericVector poss_nE, List poss_x,
                                       List poss_y, List poss_z,
                                       NumericMatrix choose_mat, int point_null,
                                       NumericVector pi_null, int point_alt,
                                       NumericVector pi_alt, int summary) {
  int           counter2,
                len_x,
                m_minus,
                m_plus,
                nC,
                nE,
                counter     = 0,
                interrupt   = 0,
                nCmax       = max(poss_nC),
                nEmax       = max(poss_nE);
  double        power,
                typeI,
                pi_power    = pi_power_finder(point_alt, pi_alt, delta),
                pi_typeI    = pi_typeI_finder(point_null, pi_null),
                theta_power = (1 - pi_power)*(pi_power + delta)/
                                (pi_power*(1 - pi_power - delta));
  NumericVector B,
                f_z,
                xC,
                xE,
                feasible(6 + nCmax + nEmax);
  NumericMatrix dbinom,
                feasible_designs(poss_nC.length(), 5 + nCmax + nEmax + 1);
  for (int n = 0; n <= poss_nC.length() - 1; n++) {
    nC                      = poss_nC[n];
    nE                      = poss_nE[n];
    if ((summary == 1) && (nC%10 == 0)) {
      Rprintf("currently analysing designs with nC = %i \n", nC);
    }
    NumericVector e_z(nC + nE + 1),
                  g_power(nC + nE + 1),
                  g_typeI(nC + nE + 1),
                  power_z(nC + nE + 1),
                  typeI_z(nC + nE + 1);
    dbinom                  = dbinom_des_one_stage(pi_typeI, pi_power, delta,
                                                   nC, nE);
    for (int z = 0; z <= nC + nE; z++) {
      interrupt++;
      if (interrupt % 1000 == 0) {
        Rcpp::checkUserInterrupt();
      }
      m_minus               = max(0, z - nC);
      m_plus                = min(z, nE);
      xE                    = seq(m_minus, m_plus);
      xC                    = z - xE;
      len_x                 = m_plus - m_minus + 1;
      B                     = xE - xC;
      NumericVector h_z_typeI(len_x),
                    h_z_power(len_x);
      for (int i = 0; i < len_x; i++) {
        h_z_typeI[i]        = choose_mat(nC - 1, xC[i])*
                                choose_mat(nE - 1, xE[i]);
        h_z_power[i]        = h_z_typeI[i]*pow(theta_power, xE[i]);
        g_power[z]         += dbinom(1, xC[i])*dbinom(3, xE[i]);
        g_typeI[z]         += dbinom(0, xC[i])*dbinom(2, xE[i]);
      }
      h_z_typeI             = h_z_typeI/sum(h_z_typeI);
      h_z_power             = h_z_power/sum(h_z_power);
      if (h_z_typeI[len_x - 1] > alpha) {
        e_z[z]              = B[len_x - 1] + 1;
      }
      else {
        counter2            = len_x - 1;
        while (typeI_z[z] <= alpha) {
          typeI_z[z]       += h_z_typeI[counter2];
          counter2--;
        }
        counter2++;
        typeI_z[z]         -= h_z_typeI[counter2];
        e_z[z]              = B[counter2 + 1];
        power_z[z]          = sum(h_z_power[Range(counter2 + 1, len_x - 1)]);
      }
    }
    f_z                     = e_z - 1;
    power                   = sum(power_z*g_power);
    if (power >= 1 - beta) {
      if (point_null == 1) {
        typeI               = sum(typeI_z*g_typeI);
        if (point_alt == 1) {
          feasible[0]                                             = nC;
          feasible[Range(1, nC + nE + 1)]                         = e_z;
          feasible[Range(2 + nCmax + nEmax, 5 + nCmax + nEmax)]   =
            NumericVector::create(pi_typeI, typeI, pi_power, power);
          feasible_designs(counter, _)                            = feasible;
          counter++;
        }
        else {
          NumericVector min_power                                 =
            fisher_min_power(1, beta, delta, NumericVector::create(nC),
                             NumericVector::create(nE), e_z, f_z,
                             NumericMatrix(0, 0),
                             List::create(poss_x[nC + nCmax*(nE - 1) - 1]),
                             List::create(poss_y[nC + nCmax*(nE - 1) - 1]),
                             List::create(poss_z[nC + nCmax*(nE - 1) - 1]),
                             pi_alt, 1);
          if (min_power[1] >= 1 - beta) {
            feasible[0]                                           = nC;
            feasible[Range(1, nC + nE + 1)]                       = e_z;
            if (nC + nE < nCmax + nEmax) {
              feasible[Range(nC + nE + 2, 1 + nCmax + nEmax)]     =
                NumericVector(nCmax + nEmax - nC - nE);
            }
            feasible[Range(2 + nCmax + nEmax, 5 + nCmax + nEmax)] =
              NumericVector::create(pi_typeI, typeI, pi_power, power);
            feasible_designs(counter, _)                          = feasible;
            counter++;
          }
        }
      }
      else {
        NumericVector max_typeI                                   =
          fisher_max_typeI(1, alpha, NumericVector::create(nC),
                           NumericVector::create(nE), e_z, f_z,
                           NumericMatrix(0, 0),
                           List::create(poss_x[nC + nCmax*(nE - 1) - 1]),
                           List::create(poss_y[nC + nCmax*(nE - 1) - 1]),
                           List::create(poss_z[nC + nCmax*(nE - 1) - 1]),
                           pi_null, 1);
        if (point_alt == 1) {
          feasible[0]                                             = nC;
          feasible[Range(1, nC + nE + 1)]                         = e_z;
          if (nC + nE < nCmax + nEmax) {
            feasible[Range(nC + nE + 2, 1 + nCmax + nEmax)]       =
              NumericVector(nCmax + nEmax - nC - nE);
          }
          feasible[Range(2 + nCmax + nEmax, 5 + nCmax + nEmax)]   =
            NumericVector::create(max_typeI[0], max_typeI[1], pi_power, power);
          feasible_designs(counter, _)                            = feasible;
          counter++;
        }
        else {
          NumericVector min_power                                 =
            fisher_min_power(1, beta, delta, NumericVector::create(nC),
                             NumericVector::create(nE), e_z, f_z,
                             NumericMatrix(0, 0),
                             List::create(poss_x[nC + nCmax*(nE - 1) - 1]),
                             List::create(poss_y[nC + nCmax*(nE - 1) - 1]),
                             List::create(poss_z[nC + nCmax*(nE - 1) - 1]),
                             pi_alt, 1);
          if (min_power[1] >= 1 - beta) {
            feasible[0]                                           = nC;
            feasible[Range(1, nC + nE + 1)]                       = e_z;
            if (nC + nE < nCmax + nEmax) {
              feasible[Range(nC + nE + 2, 1 + nCmax + nEmax)]     =
                NumericVector(nCmax + nEmax - nC - nE);
            }
            feasible[Range(2 + nCmax + nEmax, 5 + nCmax + nEmax)] =
              NumericVector::create(max_typeI[0], max_typeI[1], min_power[0], min_power[1]);
            feasible_designs(counter, _)                          = feasible;
            counter++;
          }
        }
      }
    }
  }
  NumericMatrix output                                            =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 5 + nCmax + nEmax));
  return output;
}

// [[Rcpp::export]]
List fisher_des_two_stage_cpp(double alpha, double beta, double delta,
                              NumericVector poss_nC, NumericVector poss_nE,
                              List poss_x, List poss_y, List poss_z,
                              NumericMatrix choose_mat, int point_null,
                              NumericVector pi_null, int point_alt,
                              NumericVector pi_alt, int equal,
                              int efficacy_type, double efficacy_param,
                              int futility_type, double futility_param,
                              double pi_ess, int summary) {
  int           counter2,
                len_x2,
                m_minus1,
                m_minus2,
                m_plus1,
                m_plus2,
                n1C,
                n2C,
                n1E,
                n2E,
                rows,
                counter     = 0,
                interrupt   = 0,
                nCmax       = max(poss_nC),
                nEmax       = max(poss_nE);
  double        power,
                power1,
                power2,
                S1_ess0,
                S1_ess1,
                typeI,
                typeI1,
                typeI2,
                typeII1,
                pi_power    = pi_power_finder(point_alt, pi_alt, delta),
                pi_typeI    = pi_typeI_finder(point_null, pi_null),
                theta_ess1  = (1 - pi_ess)*(pi_ess + delta)/
                                (pi_ess*(1 - pi_ess - delta)),
                theta_power = (1 - pi_power)*(pi_power + delta)/
                                (pi_power*(1 - pi_power - delta));
  NumericVector B2,
                xC1,
                xC2,
                xE1,
                xE2;
  NumericMatrix dbinom_ess,
                dbinom1,
                dbinom2;
  if (equal == 1) {
    rows                    = poss_nC.length();
  }
  else {
    rows                    = poss_nC.length()*poss_nE.length();
  }
  NumericMatrix feasible_designs_info(rows, 8),
                feasible_designs_e_z1(rows, nCmax + nEmax + 1),
                feasible_designs_f_z1(rows, nCmax + nEmax + 1),
                feasible_designs_e_z(rows, (nCmax + nEmax + 1)*
                                             (2*(nCmax + nEmax) + 1));
  int           check;
  for (int n1 = 0; n1 <= poss_nC.length() - 1; n1++) {
    check                   = 0;
    n1C                     = poss_nC[n1];
    n1E                     = poss_nE[n1];
    if (((equal != 1) && (n1C < nCmax)) ||
          ((equal == 1) && (n1C <= 0.5*nCmax))) {
      NumericVector poss_y1 = poss_y[n1C + nCmax*(n1E - 1) - 1],
                    poss_z1 = poss_z[n1C + nCmax*(n1E - 1) - 1];
      NumericMatrix poss_x1 = poss_x[n1C + nCmax*(n1E - 1) - 1];
      if ((summary == 1) && (n1C%10 == 0)) {
        Rprintf("currently analysing designs with n1C = %i \n", n1C);
      }
      NumericVector e_z1(nCmax + nEmax + 1),
                    ess0_z1(n1C + n1E + 1),
                    ess1_z1(n1C + n1E + 1),
                    f_z1(nCmax + nEmax + 1),
                    g_ess01(n1C + n1E + 1),
                    g_ess11(n1C + n1E + 1),
                    g_power1(n1C + n1E + 1),
                    g_typeI1(n1C + n1E + 1),
                    len_x1(n1C + n1E + 1),
                    power_z1(n1C + n1E + 1),
                    typeI_z1(n1C + n1E + 1);
      NumericMatrix B1(n1C + n1E + 1, n1C + n1E + 1),
                    h_typeI_z1(n1C + n1E + 1, n1C + n1E + 1),
                    h_power_z1(n1C + n1E + 1, n1C + n1E + 1),
                    h_ess1_z1(n1C + n1E + 1, n1C + n1E + 1);
      dbinom1               = dbinom_des_one_stage(pi_typeI, pi_power, delta,
                                                   n1C, n1E),
      dbinom_ess            = dbinom_des_ess(dbinom1, pi_typeI, pi_power, delta,
                                             pi_ess, n1C, n1E);
      for (int z1 = 0; z1 <= n1C + n1E; z1++) {
        interrupt++;
        if (interrupt % 1000 == 0) {
          Rcpp::checkUserInterrupt();
        }
        m_minus1                          = max(0, z1 - n1C);
        m_plus1                           = min(z1, n1E);
        len_x1[z1]                        = m_plus1 - m_minus1 + 1;
        xE1                               = seq(m_minus1, m_plus1);
        xC1                               = z1 - xE1;
        if (len_x1[z1] < n1C + n1E + 1) {
          NumericVector B1_z1(n1C + n1E + 1);
          B1_z1[Range(0, len_x1[z1] - 1)] = xE1 - xC1;
          B1(z1, _)                       = B1_z1;
        }
        else {
          B1(z1, _)                       = xE1 - xC1;
        }
        for (int i = 0; i <= len_x1[z1] - 1; i++) {
          h_typeI_z1(z1, i)               = choose_mat(n1C - 1, xC1[i])*
                                              choose_mat(n1E - 1, xE1[i]);
          h_power_z1(z1, i)               = h_typeI_z1(z1, i)*
                                              pow(theta_power, xE1[i]);
          h_ess1_z1(z1, i)                = h_typeI_z1(z1, i)*
                                              pow(theta_ess1, xE1[i]);
          g_ess01[z1]                    += dbinom_ess(0, xC1[i])*
                                              dbinom_ess(1, xE1[i]);
          g_ess11[z1]                    += dbinom_ess(0, xC1[i])*
                                              dbinom_ess(2, xE1[i]);
          g_power1[z1]                   += dbinom1(1, xC1[i])*
                                              dbinom1(3, xE1[i]);
          g_typeI1[z1]                   += dbinom1(0, xC1[i])*
                                              dbinom1(2, xE1[i]);
        }
        h_typeI_z1(z1, _)                 = h_typeI_z1(z1, _)/
                                              sum(h_typeI_z1(z1, _));
        h_power_z1(z1, _)                 = h_power_z1(z1, _)/
                                              sum(h_power_z1(z1, _));
        h_ess1_z1(z1, _)                  = h_ess1_z1(z1, _)/
                                              sum(h_ess1_z1(z1, _));
        counter2                          = len_x1[z1] - 1;
        if (efficacy_type == 0) {
          e_z1[z1]                        = z1 + 1;
        }
        else if (efficacy_type == 1) {
          if (efficacy_param == -0.5) {
            e_z1[z1]                      = round(0.5*(n1C + n1E)*delta) + 1;
          }
          else {
            e_z1[z1]                      = efficacy_param;
          }
          while ((B1(z1, counter2) >= e_z1[z1]) && (counter2 >= 0)) {
            ess0_z1[z1]                  += h_typeI_z1(z1, counter2);
            ess1_z1[z1]                  += h_ess1_z1(z1, counter2);
            typeI_z1[z1]                 += h_typeI_z1(z1, counter2);
            power_z1[z1]                 += h_power_z1(z1, counter2);
            counter2--;
          }
          if (typeI_z1[z1] > alpha) {
            check = 1;
            break;
          }
        }
        else {
          typeI_z1[z1]                    = h_typeI_z1(z1, counter2);
          if (typeI_z1[z1] > efficacy_param) {
            typeI_z1[z1]                  = 0;
            e_z1[z1]                      = z1 + 1;
          }
          else {
            ess0_z1[z1]                   = h_typeI_z1(z1, counter2);
            ess1_z1[z1]                   = h_ess1_z1(z1, counter2);
            power_z1[z1]                  = h_power_z1(z1, counter2);
            while ((typeI_z1[z1] <= efficacy_param) && (counter2 >= 1)) {
              counter2--;
              ess0_z1[z1]                += h_typeI_z1(z1, counter2);
              ess1_z1[z1]                += h_ess1_z1(z1, counter2);
              typeI_z1[z1]               += h_typeI_z1(z1, counter2);
              power_z1[z1]               += h_power_z1(z1, counter2);
            }
            ess0_z1[z1]                  -= h_typeI_z1(z1, counter2);
            ess1_z1[z1]                  -= h_ess1_z1(z1, counter2);
            typeI_z1[z1]                 -= h_typeI_z1(z1, counter2);
            power_z1[z1]                 -= h_power_z1(z1, counter2);
            e_z1[z1]                      = B1(z1, counter2 + 1);
          }
        }
        counter2                          = 0;
        typeII1                           = 0;
        if (futility_type == 0) {
          f_z1[z1]                        = -z1 - 1;
        }
        else if (futility_type == 1) {
          f_z1[z1]                        = futility_param;
          while ((B1(z1, counter2) <= futility_param) &&
                   (counter2 <= len_x1[z1] - 1)) {
            ess0_z1[z1]                  += h_typeI_z1(z1, counter2);
            ess1_z1[z1]                  += h_ess1_z1(z1, counter2);
            counter2++;
          }
        }
        else {
          typeII1                         = h_power_z1(z1, counter2);
          if (typeII1 > futility_param) {
            typeII1                       = 0;
            f_z1[z1]                      = -n1C - 1;
          }
          else {
            ess0_z1[z1]                   += h_typeI_z1(z1, counter2);
            ess1_z1[z1]                   += h_ess1_z1(z1, counter2);
            while ((typeII1 <= futility_param) && (counter2 < len_x1[z1] - 1)) {
              counter2++;
              typeII1                    += h_power_z1(z1, counter2);
              ess0_z1[z1]                += h_typeI_z1(z1, counter2);
              ess1_z1[z1]                += h_ess1_z1(z1, counter2);
            }
            ess0_z1[z1]                  -= h_typeI_z1(z1, counter2);
            ess1_z1[z1]                  -= h_ess1_z1(z1, counter2);
            f_z1[z1]                      = B1(z1, counter2 - 1);
          }
        }
        if (f_z1[z1] > e_z1[z1]) {
          check                           = 1;
          break;
        }
      }
      if (check == 0) {
        S1_ess0                           = sum(ess0_z1*g_ess01);
        S1_ess1                           = sum(ess1_z1*g_ess11);
        power1                            = sum(g_power1*power_z1);
        typeI1                            = sum(g_typeI1*typeI_z1);
        for (int n2 = (equal == 1 ? n1 : 0);
             n2 <= (equal == 1 ? n1 : poss_nC.length() - 1); n2++) {
          n2C                             = poss_nC[n2];
          if (n1C + n2C <= nCmax) {
            n2E                           = poss_nE[n2];
            NumericVector poss_y2         = poss_y[n2C + nCmax*(n2E - 1) - 1],
                          poss_z2         = poss_z[n2C + nCmax*(n2E - 1) - 1];
            NumericMatrix poss_x2         = poss_x[n2C + nCmax*(n2E - 1) - 1];
            dbinom2                       =
              dbinom_des_two_stage(dbinom1, pi_typeI, pi_power, delta, n1C, n2C,
                                   n1E, n2E);
            NumericVector g_power2(n2C + n2E + 1),
                          g_typeI2(n2C + n2E + 1);
            NumericMatrix e_z(nCmax + nEmax + 1, 2*(nCmax + nEmax) + 1),
                          power_z(n1C + n1E + 1, n2C + n2E + 1),
                          typeI_z(n1C + n1E + 1, n2C + n2E + 1);
            for (int z1 = 0; z1 <= nCmax + nEmax; z1++) {
              for (int z2 = 0; z2 <= 2*(nCmax + nEmax); z2++) {
                e_z(z1, z2)               = -0.5;
              }
            }
            typeI2                        = 0;
            power2                        = 0;
            for (int z2 = 0; z2 <= n2C + n2E; z2++) {
              m_minus2                    = max(0, z2 - n2C);
              m_plus2                     = min(z2, n2E);
              len_x2                      = m_plus2 - m_minus2 + 1;
              xE2                         = seq(m_minus2, m_plus2);
              xC2                         = z2 - xE2;
              B2                          = xE2 - xC2;
              NumericVector h_typeI_z2(len_x2),
                            h_power_z2(len_x2);
              for (int i = 0; i <= len_x2 - 1; i++) {
                h_typeI_z2[i]             = choose_mat(n2C - 1, xC2[i])*
                                              choose_mat(n2E - 1, xE2[i]);
                h_power_z2[i]             = h_typeI_z2[i]*
                                              pow(theta_power, xE2[i]);
                g_power2[z2]             += dbinom2(1, xC2[i])*
                                              dbinom2(3, xE2[i]);
                g_typeI2[z2]             += dbinom2(0, xC2[i])*
                                              dbinom2(2, xE2[i]);
              }
              h_typeI_z2                  = h_typeI_z2/sum(h_typeI_z2);
              h_power_z2                  = h_power_z2/sum(h_power_z2);
              for (int z1 = 0; z1 <= n1C + n1E; z1++) {
                NumericVector h_typeI(n1C + n2C + n1E + n2E + 1),
                              h_power(n1C + n2C + n1E + n2E + 1);
                for (int B1i = 0; B1i < len_x1[z1]; B1i++) {
                  if (B1(z1, B1i) > f_z1[z1] && B1(z1, B1i) < e_z1[z1]) {
                    for (int B2i = 0; B2i < len_x2; B2i++) {
                      h_typeI[B1(z1, B1i) + B2[B2i] + n1C + n2C] +=
                        h_typeI_z1(z1, B1i)*h_typeI_z2[B2i];
                      h_power[B1(z1, B1i) + B2[B2i] + n1C + n2C] +=
                        h_power_z1(z1, B1i)*h_power_z2[B2i];
                    }
                  }
                }
                counter2                  = n1C + n2C + n1E + n2E;
                if (h_typeI[counter2] > alpha - typeI_z1[z1]) {
                  e_z(z1, z2)             = z1 + z2 + 1;
                }
                else {
                  typeI_z(z1, z2)         = h_typeI[counter2];
                  power_z(z1, z2)         = h_power[counter2];
                  while ((typeI_z(z1, z2) <= alpha - typeI_z1[z1]) &&
                           (counter2 > 0)) {
                    counter2--;
                    typeI_z(z1, z2)      += h_typeI[counter2];
                    power_z(z1, z2)      += h_power[counter2];
                  }
                  typeI_z(z1, z2)        -= h_typeI[counter2];
                  power_z(z1, z2)        -= h_power[counter2];
                  e_z(z1, z2)             = counter2 + 1 - n1C - n2C;
                }
                typeI2                   += g_typeI1[z1]*g_typeI2[z2]*
                                              typeI_z(z1, z2);
                power2                   += g_power1[z1]*g_power2[z2]*
                                              power_z(z1, z2);
              }
            }
            power                         = power1 + power2;
            if (power >= 1 - beta) {
              if (point_null == 1) {
                typeI                     = typeI1 + typeI2;
                if (point_alt == 1) {
                  feasible_designs_info(counter, _) =
                    NumericVector::create(n1C, n2C, pi_typeI, typeI, pi_power,
                                          power, n1C + n1E +
                                                   (1 - S1_ess0)*(n2C + n2E),
                                          n1C + n1E +
                                            (1 - S1_ess1)*(n2C + n2E));
                  feasible_designs_e_z1(counter, _) = e_z1;
                  feasible_designs_f_z1(counter, _) = f_z1;
                  NumericVector e_z_vec((nCmax + nEmax + 1)*
                                          (2*(nCmax + nEmax) + 1));
                  counter2                          = 0;
                  for (int z1 = 0; z1 <= nCmax + nEmax; z1++) {
                    for (int z2 = 0; z2 <= 2*(nCmax + nEmax); z2++) {
                      e_z_vec[counter2]             = e_z(z1, z2);
                      counter2++;
                    }
                  }
                  feasible_designs_e_z(counter, _)  = e_z_vec;
                  counter++;
                }
                else {
                  NumericVector min_power             =
                    fisher_min_power(2, beta, delta,
                                     NumericVector::create(n1C, n2C),
                                     NumericVector::create(n1E, n2E), e_z1,
                                     f_z1, e_z, List::create(poss_x1, poss_x2),
                                     List::create(poss_y1, poss_y2),
                                     List::create(poss_z1, poss_z2), pi_alt, 1);
                  if (min_power[1] >= 1 - beta) {
                    feasible_designs_info(counter, _) =
                      NumericVector::create(n1C, n2C, pi_typeI, typeI,
                                            min_power[0], min_power[1],
                                            n1C + n1E +
                                              (1 - S1_ess0)*(n2C + n2E),
                                            n1C + n1E +
                                              (1 - S1_ess1)*(n2C + n2E));
                    feasible_designs_e_z1(counter, _) = e_z1;
                    feasible_designs_f_z1(counter, _) = f_z1;
                    NumericVector e_z_vec((nCmax + nEmax + 1)*
                                            (2*(nCmax + nEmax) + 1));
                    counter2                          = 0;
                    for (int z1 = 0; z1 <= nCmax + nEmax; z1++) {
                      for (int z2 = 0; z2 < 2*(nCmax + nEmax); z2++) {
                        e_z_vec[counter2]             = e_z(z1, z2);
                        counter2++;
                      }
                    }
                    feasible_designs_e_z(counter, _)  = e_z_vec;
                    counter++;
                  }
                }
              }
              else {
                NumericVector max_typeI               =
                  fisher_max_typeI(2, alpha, NumericVector::create(n1C, n2C),
                                   NumericVector::create(n1E, n2E), e_z1, f_z1,
                                   e_z,
                                 List::create(poss_x[n1C +
                                                         nCmax*(n1E - 1) - 1],
                                                poss_x[n2C +
                                                         nCmax*(n2E - 1) - 1]),
                                   List::create(poss_y[n1C +
                                                         nCmax*(n1E - 1) - 1],
                                                poss_y[n2C +
                                                         nCmax*(n2E - 1) - 1]),
                                   List::create(poss_z[n1C + 
                                                         nCmax*(n1E - 1) - 1],
                                                poss_z[n2C +
                                                         nCmax*(n2E - 1) - 1]),
                                   pi_null, 1);
                if (point_alt == 1) {
                  feasible_designs_info(counter, _)   =
                    NumericVector::create(n1C, n2C, max_typeI[0], max_typeI[1],
                                          pi_power, power,
                                          n1C + n1E + (1 - S1_ess0)*(n2C + n2E),
                                          n1C + n1E +
                                            (1 - S1_ess1)*(n2C + n2E));
                  feasible_designs_e_z1(counter, _)   = e_z1;
                  feasible_designs_f_z1(counter, _)   = f_z1;
                  NumericVector e_z_vec((nCmax + nEmax + 1)*
                                          (2*(nCmax + nEmax) + 1));
                  counter2                            = 0;
                  for (int z1 = 0; z1 <= nCmax + nEmax; z1++) {
                    for (int z2 = 0; z2 < 2*(nCmax + nEmax); z2++) {
                      e_z_vec[counter2]               = e_z(z1, z2);
                      counter2++;
                    }
                  }
                  feasible_designs_e_z(counter, _)    = e_z_vec;
                  counter++;
                }
                else {
                  NumericVector min_power             =
                    fisher_min_power(2, beta, delta,
                                     NumericVector::create(n1C, n2C),
                                     NumericVector::create(n1E, n2E), e_z1,
                                     f_z1, e_z,
                                     List::create(poss_x[n1C +
                                                           nCmax*(n1E - 1) - 1],
                                                  poss_x[n2C + nCmax*(n2E - 1) -
                                                           1]),
                                     List::create(poss_y[n1C +
                                                           nCmax*(n1E - 1) - 1],
                                                  poss_y[n2C + nCmax*(n2E - 1) -
                                                           1]),
                                     List::create(poss_z[n1C +
                                                           nCmax*(n1E - 1) - 1],
                                                  poss_z[n2C + nCmax*(n2E - 1) -
                                                           1]), pi_alt, 1);
                  if (min_power[1] >= 1 - beta) {
                    feasible_designs_info(counter, _) =
                      NumericVector::create(n1C, n2C, max_typeI[0],
                                            max_typeI[1], min_power[0],
                                            min_power[1],
                                            n1C + n1E +
                                              (1 - S1_ess0)*(n2C + n2E),
                                            n1C + n1E +
                                              (1 - S1_ess1)*(n2C + n2E));
                    feasible_designs_e_z1(counter, _) = e_z1;
                    feasible_designs_f_z1(counter, _) = f_z1;
                    NumericVector e_z_vec((nCmax + nEmax + 1)*
                                            (2*(nCmax + nEmax) + 1));
                    counter2                          = 0;
                    for (int z1 = 0; z1 <= nCmax + nEmax; z1++) {
                      for (int z2 = 0; z2 <= 2*(nCmax + nEmax); z2++) {
                        e_z_vec[counter2]             = e_z(z1, z2);
                        counter2++;
                      }
                    }
                    feasible_designs_e_z(counter, _)  = e_z_vec;
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
  List output                             =
    List::create(feasible_designs_info(
                   Range(0, 0 + (counter > 0 ? counter - 1 : 0)), Range(0, 7)),
                 feasible_designs_e_z1(
                   Range(0, 0 + (counter > 0 ? counter - 1 : 0)), _),
                 feasible_designs_f_z1(
                   Range(0, 0 + (counter > 0 ? counter - 1 : 0)), _),
                 feasible_designs_e_z(
                   Range(0, 0 + (counter > 0 ? counter - 1 : 0)), _));
  return output;
}
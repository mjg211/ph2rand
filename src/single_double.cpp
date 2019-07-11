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
NumericMatrix single_double_pmf_two_stage_cpp(NumericVector pi,
                                              NumericVector n0,
                                              NumericVector n1, double eS1,
                                              double eT1, double fS1,
                                              double fT1, double eS2,
                                              double eT2, NumericVector k) {
  
  int           counter                    = 0,
    sum_n0                     = sum(n0),
    sum_n1                     = sum(n1);
  NumericMatrix x1_prob(n0[0] + 1, n1[0] + 1),
  x2_poss(sum_n0 + 1, sum_n1 + 1),
  x2_prob(sum_n0 + 1, sum_n1 + 1),
  pmf((n0[0] + 1)*(n0[1] + 1)*(n1[0] + 1)*(n1[1] + 1), 9),
  dbinom                     = dbinom_two_stage(pi, n0, n1);
  for (int x01 = 0; x01 <= n0[0]; x01++) {
    for (int x11 = 0; x11 <= n1[0]; x11++) {
      x1_prob(x01, x11)                    = dbinom(0, x01)*dbinom(2, x11);
      if (((x11 <= fS1) || (x11 - x01 <= fT1) ||
          ((x11 >= eS1) && (x11 - x01 >= eT1))) && (k[0] == 1)) {
        pmf(counter, _)                    =
          NumericVector::create(x01, x11, n0[0], n1[0], x11, x11 - x01,
                                ((x11 >= eS1) && (x11 - x01 >= eT1)), 1,
                                x1_prob(x01, x11));
        counter++;
      }
      else if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
        for (int x02 = 0; x02 <= n0[1]; x02++) {
          for (int x12 = 0; x12 <= n1[1]; x12++) {
            x2_poss(x01 + x02, x11 + x12) += 1;
            x2_prob(x01 + x02, x11 + x12) += x1_prob(x01, x11)*dbinom(1, x02)*
              dbinom(3, x12);
          }
        }
      }
    }
  }
  if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
    for (int x0 = 0; x0 <= sum_n0; x0++) {
      for (int x1 = 0; x1 <= sum_n1; x1++) {
        if (x2_poss(x0, x1) > 0) {
          pmf(counter, _)                  =
            NumericVector::create(x0, x1, sum_n0, sum_n1, x1, x1 - x0,
                                  ((x1 >= eS2) && (x1 - x0 >= eT2)), 2,
                                  x2_prob(x0, x1));
          counter++;
        }
      }
    }
  }
  NumericMatrix output                     = pmf(Range(0, counter - 1),
                                                 Range(0, 8));
  return output;
}

// [[Rcpp::export]]
double single_double_power_one_stage(NumericVector pi, int n0, int n1, int eS,
                                     int eT, NumericMatrix poss_x,
                                     NumericVector poss_y) {
  double        power  = 0;
  NumericMatrix dbinom = dbinom_one_stage(pi, n0, n1);
  for (int o = 0; o <= (n0 + 1)*(n1 + 1) - 1; o++) {
    if ((poss_y[o] >= eT) && (poss_x(o, 1) >= eS)) {
      power           += dbinom(0, poss_x(o, 0))*dbinom(1, poss_x(o, 1));
    }
  }
  return power;
}

// [[Rcpp::export]]
double single_double_power_two_stage(NumericVector pi, NumericVector n0,
                                     NumericVector n1, NumericVector eS,
                                     NumericVector eT, NumericVector fS,
                                     NumericVector fT, List poss_x,
                                     List poss_y) {
  double        prob_x1,
                power   = 0;
  NumericVector poss_y1 = poss_y[0],
                poss_y2 = poss_y[1];
  NumericMatrix poss_x1 = poss_x[0],
                poss_x2 = poss_x[1],
                dbinom  = dbinom_two_stage(pi, n0, n1);
  for (int o1 = 0; o1 <= (n0[0] + 1)*(n1[0] + 1) - 1; o1++) {
    if ((poss_y1[o1] >= eT[0]) && (poss_x1(o1, 1) >= eS[0])) {
      power            += dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
    }
    else if ((poss_y1[o1] > fT[0]) && (poss_x1(o1, 1) > fS[0])) {
      prob_x1           = dbinom(0, poss_x1(o1, 0))*dbinom(2, poss_x1(o1, 1));
      for (int o2 = 0; o2 <= (n0[1] + 1)*(n1[1] + 1) - 1; o2++) {
        if ((poss_y1[o1] + poss_y2[o2] >= eT[1]) &&
              (poss_x1(o1, 1) + poss_x2(o2, 1) >= eS[1])) {
          power        += prob_x1*dbinom(1, poss_x2(o2, 0))*
                            dbinom(3, poss_x2(o2, 1));
        }
      }
    }
  }
  return power;
}

// [[Rcpp::export]]
NumericMatrix single_double_terminal_two_stage_cpp(NumericVector pi,
                                                   NumericVector n0,
                                                   NumericVector n1, double eS1,
                                                   double eT1, double fS1,
                                                   double fT1, double eS2,
                                                   double eT2,
                                                   NumericVector k) {
  
  int           x0,
                x1,
                counter            = 0,
                sum_n0             = sum(n0),
                sum_n1             = sum(n1);
  NumericMatrix x2_mat(sum_n0 + 1, sum_n1 + 1),
                terminal((n0[0] + 1)*(n0[1] + 1)*(n1[0] + 1)*(n1[1] + 1), 8);
  for (int x01 = 0; x01 <= n0[0]; x01++) {
    for (int x11 = 0; x11 <= n1[0]; x11++) {
      if (((x11 <= fS1) || (x11 - x01 <= fT1) ||
          ((x11 >= eS1) && (x11 - x01 >= eT1))) && (k[0] == 1)) {
        terminal(counter, _)       =
          NumericVector::create(x01, x11, n0[0], n1[0], x11, x11 - x01,
                                ((x11 >= eS1) && (x11 - x01 >= eT1)), 1);
        counter++;
      }
      else if ((k[0] == 2) || (k[k.length() - 1] == 2)) {
        for (int x02 = 0; x02 <= n0[1]; x02++) {
          for (int x12 = 0; x12 <= n1[1]; x12++) {
            if (x2_mat(x01 + x02, x11 + x12) == 0) {
              x0                   = x01 + x02;
              x1                   = x11 + x12;
              terminal(counter, _) =
                NumericVector::create(x0, x1, sum_n0, sum_n1, x1, x0 + x1,
                                      ((x1 >= eS2) && (x1  - x0 >= eT2)),
                                      2);
              x2_mat(x0, x1)      += 1;
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
NumericVector single_double_max_typeI(int J, double alpha, NumericVector n0,
                                      NumericVector n1, NumericVector eS,
                                      NumericVector eT, NumericVector fS,
                                      NumericVector fT, List poss_x,
                                      List poss_y, NumericVector pi_null,
                                      int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = pi_null[0],
                x_right = pi_null[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      -single_double_power_one_stage(NumericVector::create(v, v), n0[0], n1[0],
                                     eS[0], eT[0], poss_x[0], poss_y[0]);
  }
  else {
    f_v                 =
      -single_double_power_two_stage(NumericVector::create(v, v), n0, n1, eS,
                                     eT, fS, fT, poss_x, poss_y);
  }
  NumericVector output(4);
  if ((-f_v > alpha) && (check == 1)) {
    output              = NumericVector::create(v, -f_v, 2, iter);
    return output;
  }
  double f_z,
         z              = pi_typeI_finder(0, pi_null);
  if (J == 1) {
    f_z                 =
      -single_double_power_one_stage(NumericVector::create(z, z), n0[0], n1[0],
                                     eS[0], eT[0], poss_x[0], poss_y[0]);
  }
  else {
    f_z                 =
      -single_double_power_two_stage(NumericVector::create(z, z), n0, n1, eS,
                                     eT, fS, fT, poss_x, poss_y);
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
        -single_double_power_one_stage(NumericVector::create(u, u), n0[0],
                                       n1[0], eS[0], eT[0], poss_x[0],
                                       poss_y[0]);
    }
    else {
      f_u               =
        -single_double_power_two_stage(NumericVector::create(u, u), n0, n1, eS,
                                       eT, fS, fT, poss_x, poss_y);
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
NumericVector single_double_min_power(int J, double beta, double delta,
                                      NumericVector n0, NumericVector n1,
                                      NumericVector eS, NumericVector eT,
                                      NumericVector fS, NumericVector fT,
                                      List poss_x, List poss_y,
                                      NumericVector pi_alt, int check) {
  int           iter    = 0;
  double        f_v,
                golden  = 0.5*(3 - pow(5, 0.5)),
                x_left  = pi_alt[0],
                x_right = pi_alt[1],
                v       = x_left + golden*(x_right - x_left);
  if (J == 1) {
    f_v                 =
      single_double_power_one_stage(NumericVector::create(v, v + delta), n0[0],
                                    n1[0], eS[0], eT[0], poss_x[0], poss_y[0]);
  }
  else {
    f_v                 =
      single_double_power_two_stage(NumericVector::create(v, v + delta), n0, n1,
                                    eS, eT, fS, fT, poss_x, poss_y);
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
      single_double_power_one_stage(NumericVector::create(z, z + delta), n0[0],
                                    n1[0], eS[0], eT[0], poss_x[0], poss_y[0]);
  }
  else {
    f_z                 =
      single_double_power_two_stage(NumericVector::create(z, z + delta), n0, n1,
                                    eS, eT, fS, fT, poss_x, poss_y);
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
        single_double_power_one_stage(NumericVector::create(u, u + delta),
                                      n0[0], n1[0], eS[0], eT[0], poss_x[0],
                                      poss_y[0]);
    }
    else {
      f_u               =
        single_double_power_two_stage(NumericVector::create(u, u + delta), n0,
                                      n1, eS, eT, fS, fT, poss_x, poss_y);
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
double single_double_ess_two_stage(NumericVector pi, NumericVector n0,
                                   NumericVector n1, int eS1, int eT1, int fS1,
                                   int fT1) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, n0[0], n1[0]);
  for (int x01 = 0; x01 <= n0[0]; x01++) {
    for (int x11 = 0; x11 <= n1[0]; x11++) {
      if (((x11 - x01 >= eT1) && (x11 >= eS1)) || (x11 - x01 <= fT1) ||
          (x11 <= fS1)) {
        S1             += dbinom1(0, x01)*dbinom1(1, x11);
      }
    }
  }
  double ess            = n0[0] + n1[0] + (1 - S1)*(n0[1] + n1[1]);
  return ess;
}

// [[Rcpp::export]]
double single_double_des_ess_two_stage(NumericVector pi, NumericVector n0,
                                       NumericVector n1, int eS1, int eT1,
                                       int fS1, int fT1, NumericMatrix poss_x1,
                                       NumericVector poss_y1) {
  double        S1      = 0;
  NumericMatrix dbinom1 = dbinom_one_stage(pi, n0[0], n1[0]);
  for (int o1 = 0; o1 <= (n0[0] + 1)*(n1[0] + 1) - 1; o1++) {
    if ((poss_y1[o1] >= eT1 & poss_x1(o1, 1) >= eS1) ||
          (poss_x1(o1, 1) <= fS1) || (poss_y1[o1] <= fT1)) {
      S1               += dbinom1(0, poss_x1(o1, 0))*dbinom1(1, poss_x1(o1, 1));
    }
  }
  double ess            = n0[0] + n1[0] + (1 - S1)*(n0[1] + n1[1]);
  return ess;
}

// [[Rcpp::export]]
NumericVector single_double_max_ess_1d_two_stage(NumericVector n0,
                                                 NumericVector n1, int eS1,
                                                 int eT1, int fS1, int fT1,
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
                  -single_double_des_ess_two_stage(NumericVector::create(v, v),
                                                   n0, n1, eS1, eT1, fS1, fT1,
                                                   poss_x1, poss_y1),
                z       = 0.5,
                f_z     =
                  -single_double_des_ess_two_stage(NumericVector::create(z, z),
                                                   n0, n1, eS1, eT1, fS1, fT1,
                                                   poss_x1, poss_y1),
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
    f_u     =
      -single_double_des_ess_two_stage(NumericVector::create(u, u), n0, n1, eS1,
                                       eT1, fS1, fT1, poss_x1, poss_y1);
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
NumericMatrix single_double_des_one_stage_cpp(double alpha, double beta,
                                              double delta,
                                              NumericVector poss_n0,
                                              NumericVector poss_n1,
                                              List poss_x, List poss_y,
                                              int point_null,
                                              NumericVector pi_null,
                                              int point_alt,
                                              NumericVector pi_alt,
                                              int summary) {
  int           n0,
                n1,
                counter    = 0,
                interrupt  = 0,
                n0max      = max(poss_n0);
  double        power,
                typeI,
                pi_power   = pi_power_finder(point_alt, pi_alt, delta),
                pi_typeI   = pi_typeI_finder(point_null, pi_null);
  NumericMatrix dbinom,
                feasible_designs(10000000, 7);
  for (int n = 0; n <= poss_n0.length() - 1; n++) {
    n0                     = poss_n0[n];
    n1                     = poss_n1[n];
    if ((summary == 1) && (n0%10 == 0)) {
      Rcpp::Rcout << "  currently analysing designs with n0 = " << n0 <<
        std::endl;
    }
    NumericVector poss_y_n = poss_y[n0 + n0max*(n1 - 1) - 1];
    NumericMatrix prob_x_power(n0 + 1, n1 + 1),
                  prob_x_typeI(n0 + 1, n1 + 1),
                  poss_x_n = poss_x[n0 + n0max*(n1 - 1) - 1];
    dbinom                 = dbinom_des_one_stage(pi_typeI, pi_power, delta, n0,
                                                  n1);
    for (int o = 0; o <= (n0 + 1)*(n1 + 1) - 1; o++) {
      prob_x_power(poss_x_n(o, 0), poss_x_n(o, 1)) =
        dbinom(1, poss_x_n(o, 0))*dbinom(3, poss_x_n(o, 1));
      prob_x_typeI(poss_x_n(o, 0), poss_x_n(o, 1)) =
        dbinom(0, poss_x_n(o, 0))*dbinom(2, poss_x_n(o, 1));
    }
    for (int eT = -n0 + 1; eT <= n1; eT++) {
      for (int eS = 1; eS <= n1; eS++) {
        interrupt++;
        if (interrupt % 1000 == 0) {
          Rcpp::checkUserInterrupt();
        }
        power              = 0;
        typeI              = 0;
        for (int x1 = eS; x1 <= n1; x1++) {
          for (int x0 = 0; x0 <= n0; x0++) {
            if (x1 - x0 >= eT) {
              power       += prob_x_power(x0, x1);
              typeI       += prob_x_typeI(x0, x1);
            }
          }
        }
        if (power < 1 - beta) {
          break;
        }
        if (typeI <= alpha) {
          if (point_null == 1) {
            if (point_alt == 1) {
              feasible_designs(counter, _)   =
                NumericVector::create(n0, eS, eT, pi_typeI, typeI, pi_power,
                                      power);
              counter++;
            }
            else {
              NumericVector min_power        =
                single_double_min_power(1, beta, delta,
                                        NumericVector::create(n0),
                                        NumericVector::create(n1),
                                        NumericVector::create(eS),
                                        NumericVector::create(eT),
                                        NumericVector::create(eS - 1),
                                        NumericVector::create(eT - 1),
                                        List::create(poss_x_n),
                                        List::create(poss_y_n), pi_alt, 1);
              if (min_power[1] < 1 - beta) {
                break;
              }
              feasible_designs(counter, _)   =
                NumericVector::create(n0, eS, eT, pi_typeI, typeI, min_power[0],
                                      min_power[1]);
              counter++;
            }
          }
          else {
            NumericVector max_typeI          =
              single_double_max_typeI(1, alpha, NumericVector::create(n0),
                                      NumericVector::create(n1),
                                      NumericVector::create(eS),
                                      NumericVector::create(eT),
                                      NumericVector::create(eS - 1),
                                      NumericVector::create(eT - 1),
                                      List::create(poss_x_n),
                                      List::create(poss_y_n), pi_null, 1);
            if (max_typeI[1] <= alpha) {
              if (point_alt == 1) {
                feasible_designs(counter, _) =
                  NumericVector::create(n0, eS, eT, max_typeI[0], max_typeI[1],
                                        pi_power, power);
                counter++;
              }
              else {
                NumericVector min_power      =
                  single_double_min_power(1, beta, delta,
                                          NumericVector::create(n0),
                                          NumericVector::create(n1),
                                          NumericVector::create(eS),
                                          NumericVector::create(eT),
                                          NumericVector::create(eS - 1),
                                          NumericVector::create(eT - 1),
                                          List::create(poss_x_n),
                                          List::create(poss_y_n), pi_alt, 1);
                if (min_power[1] < 1 - beta) {
                  break;
                }
                feasible_designs(counter, _) =
                  NumericVector::create(n0, eS, eT, max_typeI[0], max_typeI[1],
                                        min_power[0], min_power[1]);
                counter++;
              }
            }
          }
        }
      }
    }
  }
  NumericMatrix output                       =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 6));
  return output;
}

// [[Rcpp::export]]
NumericMatrix single_double_des_two_stage_cpp(double alpha, double beta,
                                              double delta,
                                              NumericVector poss_n0,
                                              NumericVector poss_n1,
                                              List poss_x, List poss_y,
                                              int point_null,
                                              NumericVector pi_null,
                                              int point_alt,
                                              NumericVector pi_alt, int equal,
                                              int efficacy, int futility,
                                              double pi_ess, int summary) {
  int           fT2,
                fS2,
                n01,
                n11,
                n02,
                n12,
                counter                   = 0,
                interrupt                 = 0,
                n0max                     = max(poss_n0);
  double        power1,
                power2,
                S1_ess0,
                S1_ess1,
                typeI1,
                typeI2,
                typeII1,
                pi_power                  = pi_power_finder(point_alt, pi_alt,
                                                            delta),
                pi_typeI                  = pi_typeI_finder(point_null,
                                                            pi_null);
  NumericMatrix dbinom_ess,
                dbinom1,
                dbinom2,
                feasible_designs(10000000, 14);
  for (int n1 = 0; n1 <= poss_n0.length() - 1; n1++) {
    n01                                   = poss_n0[n1];
    n11                                   = poss_n1[n1];
    if (((equal != 1) && (n01 < n0max)) || ((equal == 1) && (n01 <= 0.5*n0max))) {
      if ((summary == 1) && (n01%10 == 0)) {
        Rcpp::Rcout << "  currently analysing designs with n01 = " << n01 <<
          std::endl;
      }
      NumericVector poss_y1               = poss_y[n01 + n0max*(n11 - 1) - 1];
      NumericMatrix prob_x1_ess0(n01 + 1, n11 + 1),
                    prob_x1_ess1(n01 + 1, n11 + 1),
                    prob_x1_power(n01 + 1, n11 + 1),
                    prob_x1_typeI(n01 + 1, n11 + 1),
                    poss_x1               = poss_x[n01 + n0max*(n11 - 1) - 1],
                    dbinom1               =
                      dbinom_des_one_stage(pi_typeI, pi_power, delta, n01, n11),
                    dbinom_ess            =
                      dbinom_des_ess(dbinom1, pi_typeI, pi_power, delta,
                                      pi_ess, n01, n11);
      for (int o1 = 0; o1 <= (n01 + 1)*(n11 + 1) - 1; o1++) {
        prob_x1_ess0(poss_x1(o1, 0), poss_x1(o1, 1))  =
          dbinom_ess(0, poss_x1(o1, 0))*dbinom_ess(1, poss_x1(o1, 1));
        prob_x1_ess1(poss_x1(o1, 0), poss_x1(o1, 1))  =
          dbinom_ess(0, poss_x1(o1, 0))*dbinom_ess(2, poss_x1(o1, 1));
        prob_x1_power(poss_x1(o1, 0), poss_x1(o1, 1)) =
          dbinom1(1, poss_x1(o1, 0))*dbinom1(3, poss_x1(o1, 1));
        prob_x1_typeI(poss_x1(o1, 0), poss_x1(o1, 1)) =
          dbinom1(0, poss_x1(o1, 0))*dbinom1(2, poss_x1(o1, 1));
      }
      for (int fT1 = (futility == 1 ? -n01 : -n01 - 1);
           fT1 <= (futility == 1 ? (efficacy == 1 ? n11 - 2 : n11 - 1) :
                     -n01 - 1); fT1++) {
        for (int fS1 = (futility == 1 ? 0 : -1);
             fS1 <= (futility == 1 ? (efficacy == 1 ? n11 - 2 : n11 - 1) : -1);
             fS1++) {
          typeII1                         = 0;
          for (int x01 = 0; x01 <= n01; x01++) {
            for (int x11 = 0; x11 <= n11; x11++) {
              if ((x11 - x01 <= fT1) || (x11 <= fS1)) {
                typeII1                  += prob_x1_power(x01, x11);
              }
            }
          }
          if (typeII1 > beta) {
            break;
          }
          for (int eT1 = (efficacy == 1 ? fT1 + 2 : n11 + 1);
               eT1 <= (efficacy == 1 ? n11 : n11 + 1); eT1++) {
            for (int eS1 = (efficacy == 1 ? fS1 + 2 : n11 + 1);
                 eS1 <= (efficacy == 1 ? n11 : n11 + 1); eS1++) {
              power1                      = 0;
              typeI1                      = 0;
              for (int x11 = eS1; x11 <= n11; x11++) {
                for (int x01 = 0; x01 <= n01; x01++) {
                  if (x11 - x01 >= eT1) {
                    power1               += prob_x1_power(x01, x11);
                    typeI1               += prob_x1_typeI(x01, x11);
                  }
                }
              }
              if (typeI1 < alpha) {
                S1_ess0                   = 0;
                S1_ess1                   = 0;
                for (int o1 = 0; o1 <= (n01 + 1)*(n11 + 1) - 1; o1++) {
                  if (((poss_y1[o1] >= eT1) && (poss_x1(o1, 1) >= eS1)) ||
                      (poss_y1[o1] <= fT1) || (poss_x1(o1, 1) <= fS1)) {
                    S1_ess0              += prob_x1_ess0(poss_x1(o1, 0),
                                                         poss_x1(o1, 1));
                    S1_ess1              += prob_x1_ess1(poss_x1(o1, 0),
                                                         poss_x1(o1, 1));
                  }
                }
                for (int n2 = (equal == 1 ? n1 : 0);
                     n2 <= (equal == 1 ? n1 : poss_n0.length() - 1); n2++) {
                  n02                     = poss_n0[n2];
                  if (n01 + n02 <= n0max) {
                    n12                   = poss_n1[n2];
                    NumericVector poss_y2 = poss_y[n02 + n0max*(n12 - 1) - 1];
                    NumericMatrix prob_x_power(n01 + n02 + 1, n11 + n12 + 1),
                                  prob_x_typeI(n01 + n02 + 1, n11 + n12 + 1),
                                  poss_x2 = poss_x[n02 + n0max*(n12 - 1) - 1];
                    dbinom2               =
                      dbinom_des_two_stage(dbinom1, pi_typeI, pi_power, delta,
                                           n01, n02, n11, n12);
                    for (int o1 = 0; o1 <= (n01 + 1)*(n11 + 1) - 1; o1++) {
                      if ((poss_y1[o1] > fT1) && (poss_x1(o1, 1) > fS1) &&
                          ((poss_y1[o1] < eT1) || (poss_x1(o1, 1) < eS1))) {
                        for (int o2 = 0; o2 <= (n02 + 1)*(n12 + 1) - 1; o2++) {
                          prob_x_power(poss_x1(o1, 0) + poss_x2(o2, 0),
                                       poss_x1(o1, 1) + poss_x2(o2, 1)) +=
                            prob_x1_power(poss_x1(o1, 0), poss_x1(o1, 1))*
                            dbinom2(1, poss_x2(o2, 0))*
                            dbinom2(3, poss_x2(o2, 1));
                          prob_x_typeI(poss_x1(o1, 0) + poss_x2(o2, 0),
                                       poss_x1(o1, 1) + poss_x2(o2, 1)) +=
                            prob_x1_typeI(poss_x1(o1, 0), poss_x1(o1, 1))*
                            dbinom2(0, poss_x2(o2, 0))*
                            dbinom2(2, poss_x2(o2, 1));
                        }
                      } 
                    }
                    for (int eT2 = fT1 + 2 - n02; eT2 <= eT1 - 1 + n12; eT2++) {
                      fT2                 = eT2 - 1;
                      for (int eS2 = fS1 + 2; eS2 <= eS1 - 1 + n12; eS2++) {
                        interrupt++;
                        if (interrupt % 1000 == 0) {
                          Rcpp::checkUserInterrupt();
                        }
                        fS2               = eS2 - 1;
                        power2            = 0;
                        typeI2            = 0;
                        for (int x1 = eS2; x1 <= n11 + n12; x1++) {
                          for (int x0 = 0; x0 <= n01 + n02; x0++) {
                            if (x1 - x0 >= eT2) {
                              power2     += prob_x_power(x0, x1);
                              typeI2     += prob_x_typeI(x0, x1);
                            }
                          }
                        }
                        if (power1 + power2 < 1 - beta) {
                          break;
                        }
                        if (typeI1 + typeI2 <= alpha) {
                          if (point_null == 1) {
                            if (point_alt == 1) {
                              feasible_designs(counter, _)   =
                                NumericVector::create(n01, n02, eS1, eS2, eT1,
                                                      eT2, fS1, fT1, pi_typeI,
                                                      typeI1 + typeI2, pi_power,
                                                      power1 + power2,
                                                      n01 + n11 + (1 - S1_ess0)*
                                                                    (n02 + n12),
                                                      n01 + n11 +
                                                        (1 - S1_ess1)*
                                                        (n02 + n12));
                              counter++;
                            }
                            else {
                              NumericVector min_power        =
                                single_double_min_power(
                                  2, beta, delta,
                                  NumericVector::create(n01, n02),
                                  NumericVector::create(n11, n12),
                                  NumericVector::create(eS1, eS2),
                                  NumericVector::create(eT1, eT2),
                                  NumericVector::create(fS1, fS2),
                                  NumericVector::create(fT1, fT2),
                                  List::create(poss_x1, poss_x2),
                                  List::create(poss_y1, poss_y2), pi_alt, 1);
                              if (min_power[1] < 1 - beta) {
                                break;
                              }
                              feasible_designs(counter, _)   =
                                NumericVector::create(
                                  n01, n02, eS1, eS2, eT1, eT2, fS1, fT1,
                                  pi_typeI, typeI1 + typeI2, min_power[0],
                                  min_power[1],
                                  n01 + n11 + (1 - S1_ess0)*(n02 + n12),
                                  n01 + n11 + (1 - S1_ess1)*(n02 + n12));
                              counter++;
                            }
                          }
                          else {
                            NumericVector max_typeI          =
                              single_double_max_typeI(
                                2, alpha, NumericVector::create(n01, n02),
                                NumericVector::create(n11, n12),
                                NumericVector::create(eS1, eS2),
                                NumericVector::create(eT1, eT2),
                                NumericVector::create(fS1, fS2),
                                NumericVector::create(fT1, fT2),
                                List::create(poss_x1, poss_x2),
                                List::create(poss_y1, poss_y2), pi_null, 1);
                            if (max_typeI[1] <= alpha) {
                              if (point_alt == 1) {
                                feasible_designs(counter, _) =
                                  NumericVector::create(n01, n02, eS1, eS2, eT1,
                                                        eT2, fS1, fT1,
                                                        max_typeI[0],
                                                        max_typeI[1],
                                                        pi_power,
                                                        power1 + power2,
                                                        n01 + n11 +
                                                          (1 - S1_ess0)*
                                                          (n02 + n12),
                                                        n01 + n11 +
                                                          (1 - S1_ess0)*
                                                          (n02 + n12));
                                counter++;
                              }
                              else {
                                NumericVector min_power      =
                                  single_double_min_power(
                                    2, beta, delta,
                                    NumericVector::create(n01, n02),
                                    NumericVector::create(n11, n12),
                                    NumericVector::create(eS1, eS2),
                                    NumericVector::create(eT1, eT2),
                                    NumericVector::create(fS1, fS2),
                                    NumericVector::create(fT1, fT2),
                                    List::create(poss_x1, poss_x2),
                                    List::create(poss_y1, poss_y2), pi_alt, 1);
                                if (min_power[1] < 1 - beta) {
                                  break;
                                }
                                feasible_designs(counter, _) =
                                  NumericVector::create(n01, n02, eS1, eS2, eT1,
                                                        eT2, fS1, fT1,
                                                        max_typeI[0],
                                                        max_typeI[1],
                                                        min_power[0],
                                                        min_power[1],
                                                        n01 + n11 +
                                                          (1 - S1_ess0)*
                                                          (n02 + n12),
                                                        n01 + n11 +
                                                          (1 - S1_ess1)*
                                                          (n02 + n12));
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
      }
    }
  }
  NumericMatrix output                    =
    feasible_designs(Range(0, 0 + (counter > 0 ? counter - 1 : 0)),
                     Range(0, 13));
  return output;
}
#include <algorithm>
#include <cmath>
#include <glpk.h>
#include <iostream>
#include <sstream>
#include <string>

#include "PASAQ.h"
#include "lin_prog.h"

using std::cout;
using std::endl;
using std::to_string;
using std::string;

// #define DEBUG

/*
 * Expected Attacker utility on attacking target i, given defender strategy x.
 */
double U_a(const size_t i, const strategy &x, const PayoffMatrix &Pm) {
  return x[i] * Pm.P_a[i] + (1 - x[i]) * Pm.R_a[i];
}

/*
 * Expected defender utility on if the adversary attacks target i, given
 * strategy x. 
 */
double U_d(const size_t i, const strategy &x, const PayoffMatrix &Pm) {
  return x[i] * Pm.R_d[i] + (1 - x[i]) * Pm.P_d[i];
}

double q_i(const size_t i, const strategy &s, const PayoffMatrix Pm,
           const double lambda) {
  const size_t T = Pm.P_a.size();
  double sum = 0;
  for (size_t j = 0; j < T; j++)
    sum += exp(lambda * U_a(j, s, Pm));
  return exp(lambda * U_a(i, s, Pm)) / sum;
}

/* 
 * Expected defender utility. 
 */
double UD(const strategy &x, const PayoffMatrix &Pm, const double lambda) {
  double sum = 0;
  for (size_t i = 0; i < x.size(); i++)
    sum += q_i(i, x, Pm, lambda) * U_d(i, x, Pm);
  return sum;
}

/* 
 * Expected Attacker utility.
 */
double UA(const strategy &x, const PayoffMatrix &Pm, const double lambda) {
  double sum = 0;
  for (size_t i = 0; i < x.size(); i++)
    sum += q_i(i, x, Pm, lambda) * U_a(i, x, Pm);
  return sum;
}

/* 
 * Symbols for targets in SSG. 
 */
double theta(const int i, const PayoffMatrix &Pm, const double lambda) {
  return exp(lambda * Pm.R_a[i]);
}

double alpha(const int i, const PayoffMatrix &Pm, const double lambda) {
  return (Pm.R_d[i] - Pm.P_d[i]);
}

double beta(const int i, const PayoffMatrix &Pm, const double lambda) {
  return lambda * (Pm.R_a[i] - Pm.P_a[i]);
}

double f1(int i, double x, const PayoffMatrix &Pm, const double lambda) {
  return exp(-beta(i, Pm, lambda) * x);
}

double f2(int i, double x, const PayoffMatrix &Pm, const double lambda) {
  return x * exp(-beta(i, Pm, lambda) * x);
}

/* 
 * Estimate the upper and lower bound of utility the defender can achieve. The
 * lower bound is set to the expected utiltiy of a uniform strategy. The upper
 * bound is set to the sum of all rewards.
 */
pair<double, double> EstimateBounds(int numRes, const PayoffMatrix &Pm,
                                    const double lambda) {
  const int T = Pm.P_a.size();
  pair<double, double> result = {0, 0};
  strategy uniform_strategy(
      T, std::max(numRes * (1.0 / static_cast<double>(T)), 1.0));
  result.first = UD(uniform_strategy, Pm, lambda);
  for (const auto &reward : Pm.R_d) 
    result.second += reward;
  return result;
}

/* Helper functions that return column index of the variables in the LP. */

/*
 * Index of variable x_ik
 */
inline int x_index(const int T, const int K, const int i, const int k) {
  return (i-1)* K + k;
}

/*
 * Index of variable z_ik
 */
inline int z_index(const int T, const int K, const int i, const int k) {
  return (T * K) + (i - 1) * K + k;
}

/* 
 * Index of variable a_j
 */
inline int a_index(const int T, const int K, const int j) {
  return (T * K) * 2 + j;
}

/* 
 * Helper functions to set PASAQ constraints as defined in paper. 
 *
 * Constraint (11): SUM x_ik < M 
 * Constraint (12): EACH 0 < x_ik < 1/K 
 * Constraint (13): EACH zik * (1/k) < xik => zik * (1/k) - xik < 0 
 * Constraint (14): EACH x_ik+1 <= z_ik => EACH xx_ik+1 - z_ik <= 0. 
 * Constraint (15): EACH z_ik is either 0 or 1
 * Constraint (16): SUM x_ik = SUM a_j * A_ij. 
 * Constraint (17): SUM a_j = 1
 * Constraint (18): EACH a_j is between 0 and 1
 */
void set_pasaq_constraint_11(lin_prog &LP, const size_t T, const size_t K,
                               const int num_res);
void set_pasaq_constraint_12(lin_prog &LP, const size_t T, const size_t K);
void set_pasaq_constraint_13(lin_prog &LP, const size_t T, const size_t K);
void set_pasaq_constraint_14(lin_prog &LP, const size_t T, const size_t K);
void set_pasaq_constraint_15(lin_prog &LP, const size_t T, const size_t K);

// PASAQ with assignment constraints.
void set_pasaq_constraint_16(lin_prog &LP, const size_t T, const size_t K,
                             const vector<vector<double>> &A);
void set_pasaq_constraint_17(lin_prog &LP, const size_t T, const size_t K,
                             const vector<vector<double>> &A);
void set_pasaq_constraint_18(lin_prog &LP, const size_t T, const size_t K,
                             const vector<vector<double>> &A);

size_t set_pasaq_constraints(const size_t T, const size_t K);

/*
 * Set objective function for a PASAQ problem with constraints within a binary
 * search method
 *
 * Arguments:
 * lp - problem object
 * r - 
 */
void set_pasaq_obj(lin_prog &lp, const double r, const PayoffMatrix &Pm,
                   const double lambda, const int K);


void print_lp_result(int result);

/*
 * Generate CF-OPT and solve it using GPLK, to check that a strategy is feasible
 * and return such a strategy. We do this by creating a linear program defined
 * by PASAQ with assignment constraints.
 *
 * A = Probability matrix.
 */
pair<bool, vector<double>> CheckFeasibility(const double r, const int num_res,
                                            const PayoffMatrix &Pm,
                                            const vector<vector<double>> &A,
                                            const double K, double lambda) {
  const size_t T = Pm.P_a.size();
  cout << "CheckFeasibility(" << r << ");" << endl;
  pair<bool, vector<double>> result;
  std::string name = "Feasibility check r = " + to_string(r);
  cout << "\tT = " << T << " K=" << K << " A=" << A.size() << "x" << A[1].size()
       << endl;

  lin_prog LP("Check Feasibility r = " + to_string(r));
  LP.declare_variables("x", T*K);
  LP.declare_variables("z", T*K);
  LP.declare_variables("a", A[1].size());  

  set_pasaq_obj(LP, r, Pm, lambda, K);
  set_pasaq_constraint_11(LP, T, K, num_res);
  set_pasaq_constraint_12(LP, T, K);
  set_pasaq_constraint_13(LP, T, K);
  set_pasaq_constraint_14(LP, T, K);
  set_pasaq_constraint_15(LP, T, K);
  set_pasaq_constraint_16(LP, T, K, A);
  set_pasaq_constraint_17(LP, T, K, A);
  set_pasaq_constraint_18(LP, T, K, A);
  
  glp_iocp parm;
  parm.presolve = GLP_ON;
  LP.run(&parm);

  double obj_val = LP.get_obj_val();
  cout << "obj value = " << obj_val << endl;

  result.first = obj_val == 0;
  result.second = vector<double>(T);
  
  std::cout << "\nVariable x values:" << "\n";
  for (size_t i = 1; i <= T; i++) {
    double sum = 0;
    for (int k = 1; k <= K; k++) {
      auto x = LP.get_var_val("x",  (i - 1) * K + k);
      // cout << "x_" << i << k << "=" << x << endl;
      sum += x;
    }
    result.second[i] = sum;
    cout << "x_" << i << "=" << result.second[i] << (i % 5 == 0 ? "\n" : " ");
  }

  std::cout << "\nVariable z values:" << "\n";
  for (size_t i = 1; i <= T; i++) {
    for (size_t k = 1; k <= K; k++) {
      auto z = LP.get_var_val("z",  (i - 1) * K + k);
      cout << "z_{" << i << "," << k << "}=" << z << " ";
    }
    cout << endl;
  }

  // cout << "PRINTING AS: " << A[1].size() << endl;
  std::cout << "\nVariable a values:" << "\n";
  for (size_t j = 1; j <= A[1].size(); j++) {
    auto a =  LP.get_var_val("a", j);
    cout << "a_" << j << "=" << a << (j % 5 == 0 ? "\n" : " ");
  }
  cout << endl;
  return result;
}

// Main algorithm for finding strategy.
pair<double, vector<double>>
BinarySearchMethod(const double e, const int numRes, const PayoffMatrix &Pm,
                   const vector<vector<double>> &A, const double lambda,
                   const double K) {
  cout << "BinarySearchMethod(" << e << ", " << numRes << ")" << endl;
  const auto pair = EstimateBounds(numRes, Pm, lambda);
  auto L = pair.first;
  auto U = pair.second;
  vector<double> x;
  cout << "U = " << U << " L=" << L << endl;
  while (U - L > e) {
    double r = (U + L) / 2;
    cout << "U = " << U << " L=" << L << " r = " << r << endl;
    const auto f_x_pair = CheckFeasibility(r, numRes, Pm, A, K, lambda);
     x = f_x_pair.second;
    if (f_x_pair.first) {
      L = r;
    } else {
      U = r;
    }
  }
  return std::pair<double, vector<double>>(L, x);
}

void set_pasaq_obj(lin_prog  &LP, const double r, const PayoffMatrix &Pm,
                   const double lambda, const int K) {
  const double T = Pm.P_a.size();
#ifdef DEBUG
  cout << "OBJECTIVE:";
#endif
  LP.set_max();
  for (int i = 1; i <= T; i++) {
    const double theta_ = theta(i, Pm, lambda);
    const double alpha_ = alpha(i, Pm, lambda);
    const double coef = theta_ * (r - Pm.P_d[i]);
    for (int k = 1; k <= K; k++) {
      const double k_ = static_cast<double>(k);
      const double left = (k_ - 1.0) / (double)K;
      const double right = k_ / (double)K;
      const double y_ik =
          (f1(i, right, Pm, lambda) - f1(i, left, Pm, lambda)) / (right - left);
      const double u_ik =
          (f2(i, right, Pm, lambda) - f2(i, left, Pm, lambda)) / (right - left);
      double coef_val = coef * y_ik;
      coef_val = // coef_val
        0 - (theta_ * alpha_ * u_ik);
      LP.set_objective_var("x", ((i - 1) * K) + k, coef_val);
#ifdef DEBUG
      cout << (k > 1 ? " " : "\n");
      cout << "x_{" << i << "," << k << "} = " << coef_val;
      cout << "y_" << i << k << " = " << y_ik << endl;
      cout << "u_" << i << k << " = " << u_ik << endl;
      cout << "r = " << r << endl;
      cout << "theta_i=" << theta_ << endl;
      cout << "alpha_i=" << alpha_ << endl;
#endif
    }
#ifdef DEBUG
    cout << endl;
#endif
  }
}

void set_pasaq_constraint_11(lin_prog &LP, const size_t T, const size_t K,
                               const int num_res) {
  // glp_set_row_name(lp, 1, "11");
  LP.add_row("(11)");
  LP.set_row_bnd(GLP_UP, 0, num_res);
  // glp_set_row_bnds(lp, 1, GLP_UP, 0, num_res);
  for (size_t i = 1; i <= T; i++) {
    for (size_t k = 1; k <= K; k++) {
      LP.add_constraint("x", ((i - 1) * K) + k, 1);
      // cm.row_index.push_back(1);
      // cm.col_index.push_back(((i - 1) * K) + k);
      // cm.value.push_back(1);
    }
  }
}

void set_pasaq_constraint_12(lin_prog& LP, const size_t T, const size_t K) {
  for (size_t i = 1; i <= T; i++) {
    for (size_t k = 1; k <= K; k++) {  
    std::ostringstream col_name;
      col_name << "x_" << i << k;
      LP.set_var_bnd("x", ((i - 1) * K) + k, GLP_DB, 0,
                     1.0 / static_cast<double>(K));
      // glp_set_col_name(lp, ((i - 1) * K) + k, col_name.str().c_str());
      // glp_set_col_bnds(lp, ((i - 1) * K) + k, GLP_DB, 0,
      //                  1.0 / static_cast<double>(K));
    }
  }
}

/*
 * This constraint requires each z_{ik} / K <= x_{ik} (for all i and K). Thus,
 * this constraint must add T*K rows to our lp.
 *
 * Constraint 13: EACH zik * (1/k) < xik => zik * (1/k) - xik < 0 
 */
void set_pasaq_constraint_13(lin_prog &LP, const size_t T, const size_t K) {
  for (size_t i = 1; i <= T; i++) {
    for (size_t k = 1; k <= K; k++) {
      LP.add_row("13-" + to_string(i) + " " + to_string(k));
      // string row_name = "13-" + to_string(i) + " " + to_string(k);
      // glp_set_row_name(lp, ++last_row_used, row_name.c_str());
      LP.set_row_bnd(GLP_UP, 0, 0);
      LP.add_constraint("x", ((i - 1) * K) + k,-1);
      LP.add_constraint("z", ((i - 1) * K) + k, (1 / static_cast<double>(K)));
      // cm.col_index.push_back(z_offset + ((i - 1) * K) + k);
      // glp_set_row_bnds(lp, last_row_used, GLP_UP, 0, 0);
      // cm.row_index.push_back(last_row_used);
      // cm.col_index.push_back(((i - 1) * K) + k);
      // cm.value.push_back(-1); 
      // cm.row_index.push_back(last_row_used);
      // cm.col_index.push_back(z_offset + ((i - 1) * K) + k);
      // cm.value.push_back(1 / (double)K);
    }
  }
}

void set_pasaq_constraint_14(lin_prog &LP, const size_t T, const size_t K) {
  for (size_t i = 1; i <= T; i++) {
    for (size_t k = 1; k <= K - 1; k++) {
      LP.add_row("14-"+ std::to_string(i) + std::to_string(k));
      LP.set_row_bnd(GLP_UP, 0, 0);
      LP.add_constraint("x", ((i - 1) * K) + k + 1, 1);
      LP.add_constraint("z", ((i - 1) * K) + k, -1);
      // LP.add_row("14-" + std::to_string(i) + std::to_string(k));
      // glp_set_row_bnds(lp, row_index++, GLP_UP, 0, 0);
      // cm.row_index.push_back(row_index);
      // cm.col_index.push_back(x_offset + ((i - 1) * K) + k);
      // cm.value.push_back(1);
      // cm.row_index.push_back(row_index);
      // cm.col_index.push_back(z_offset + ((i - 1) * K) + k);
      // cm.value.push_back(-1);
    }
  }
}

void set_pasaq_constraint_15(lin_prog &LP, const size_t T, const size_t K) {
  for (size_t i = 1; i <= T; i++) {
    for (size_t k = 1; k <= K; k++) {
      // Each z is a binary variables, 0 or 1.
      LP.set_var_kind("z", ((i - 1) * K) + k,GLP_BV);
      // the follwoing should be unnecessary, as the var kind enforces this.
      LP.add_row("15-"+ std::to_string(i) + std::to_string(k));
      LP.add_constraint("z", ((i - 1) * K) + k, 1);
      LP.set_row_bnd(GLP_DB, 0, 1);
    }
  }
}

void set_pasaq_constraint_16(lin_prog &LP, const size_t T, const size_t K,
                             const vector<vector<double>> &A) {
  for (size_t i = 1; i < T; i++) {
    LP.add_row("16-" + std::to_string(i));
    LP.set_row_bnd(GLP_FX, 0, 0);
    for (size_t k = 1; k <= K; k++) {
      LP.add_constraint("x", ((i - 1) * K + k), 1);
    }
    for (size_t j = 1; j <= A[1].size(); j++) {
      LP.add_constraint("a", j, -A[i][j]);
    }
  }
}

void set_pasaq_constraint_17(lin_prog &LP, const size_t T, const size_t K,
                             const vector<vector<double>> &A) {
  // Set bounds for a_j
  LP.add_row("(17)");
  LP.set_row_bnd(GLP_UP,0,1);
  for (size_t  j = 1; j <= A[1].size(); j++)
    LP.add_constraint("a", j, 1);
}

void set_pasaq_constraint_18(lin_prog &LP, const size_t T, const size_t K,
                             const vector<vector<double>> &A) {
  for (size_t j = 1; j < A[1].size(); j++) {
    LP.set_var_bnd("a", j, GLP_DB, 0, 1);
 }
}

void print_lp_result(int result) {
  // glp_write_lp(lp, NULL, "logs/log.txt");
  cout << "MIP GLPK RESULT " << result << endl;
  switch (result) {
  case 0:
    cout << "MILP SUCCESS" << endl;
    break;
  case GLP_EBOUND:
    cout << "MILP EBOUND" << endl;
    break;
  case GLP_EROOT:
    cout << "MILP EROOT" << endl;
    break;
  case GLP_ENOPFS:
    cout << "MILP ENOPFS" << endl;
    break;
  case GLP_ENODFS:
    cout << "MILP ENODFS" << endl;
    break;
  case GLP_EFAIL:
    cout << "MILP EFAIL" << endl;
    break;
  case GLP_EMIPGAP:
    cout << "MILP EMIPGAP" << endl;
    break;
  case GLP_ETMLIM:
    cout << "MILP ETMLIM" << endl;
    break;
  case GLP_ESTOP:
    cout << "MILP ESTOP" << endl;
    break;
  default:
    cout << "MILP: UNKOWN" << endl;
    break;
  }  
}

#ifndef PASAQ_H
#define PASAQ_H

#include <utility>
#include <vector>

using std::vector;
using std::pair;

typedef vector<double> QuantalResponse;
typedef vector<double> strategy;
typedef vector<int> Payoff;

// Payoff Matrix for a game.
struct PayoffMatrix {
  Payoff R_d; // Defender reward.
  Payoff P_d; // Defender penalty.
  Payoff R_a; // Attacker reward.
  Payoff P_a; // Attacker penalty.
  PayoffMatrix(const Payoff &attacker_reward, const Payoff &attacker_penalty,
               const Payoff &defender_reward, const Payoff &defender_penalty)
      : R_d(defender_reward), P_d(defender_penalty), R_a(attacker_reward),
        P_a(attacker_penalty) {}
};

double q_i(const int i, const strategy &s, const double lambda);

pair<double, vector<double>>
BinarySearchMethod(const double e, const int numRes, const PayoffMatrix &Pm,
                   const vector<vector<double>> &A, const double lambda,
                   const double K);

#endif /* PASAQ_H */

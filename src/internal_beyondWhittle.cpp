#include <Rcpp.h>
using namespace Rcpp;

//' C++ function for generating  p from v in Stick Breaking DP representation
//' @keywords internal
// [[Rcpp::export]]
NumericVector pFromV(NumericVector v) {
  unsigned L = v.size();
  NumericVector p(L + 1);
  double currentProduct = 1.0;
  double pSum = 0.0;
  for (unsigned l = 0; l < L; ++l) {
    p[l + 1] = currentProduct * v[l];
    currentProduct *= (1.0 - v[l]);
    pSum += p[l + 1];
  }
  p[0] = std::max(1.0 - pSum, 0.0); // account for numerical instabilities
  return p;
}

//' C++ function for computing mixture weights of Bernstein-Mixtures given the probabilities p, values w, and degree k.
//' @keywords internal
// [[Rcpp::export]]
NumericVector mixtureWeight(NumericVector p, NumericVector w, unsigned k) {
  typedef std::pair<double, double> wpType;
  std::vector<wpType> wp;
  for (unsigned l = 0; l < p.size(); ++l) {
    wp.push_back(wpType(w[l], p[l]));
  }
  std::sort(wp.begin(), wp.end());
  NumericVector weight(k);
  unsigned l = 0;
  for (unsigned j = 1; j <= k; ++j) {
    weight[j-1] = 0;
    double wMax = j / (double)k;
    while (l < wp.size() && wp[l].first <= wMax) {
      weight[j-1] += wp[l].second;
      l += 1;
    }
  }
  return weight;
}

//' C++ function for building a density mixture, given mixture weights and functions.
//' @keywords internal
// [[Rcpp::export]]
NumericVector densityMixture(NumericVector weights, NumericMatrix densities) {
  if (weights.size() != densities.nrow()) {
    return(NumericVector());
  }
  const unsigned n = densities.ncol();
  NumericVector res(n);
  for (unsigned omega = 0; omega < n; ++omega) {
    res[omega] = 0.0;
  }
  for (unsigned j = 0; j < weights.size(); ++j) {
    for (unsigned omega = 0; omega < n; ++omega) {
      res[omega] += weights[j] * densities(j, omega);
    }
  }
  return(res);
}
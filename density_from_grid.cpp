#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector
density_from_grid_cpp(NumericVector data, NumericVector h, NumericVector grid,
                      Nullable<NumericVector> weights_ = R_NilValue,
                      double cutoff = 6.0) {
  int n_data = data.size();
  int G = grid.size();

  NumericVector weights;

  // Handle weights
  if (weights_.isNull()) {
    weights = NumericVector(n_data, 1.0);
  } else {
    weights = weights_.get();

    if (weights.size() != n_data)
      stop("weights and data must have same length");
  }

  double weight_sum = 0.0;
  for (int i = 0; i < n_data; i++) {
    weight_sum += weights[i];
  }

  double norm_const = 1.0 / std::sqrt(2.0 * M_PI);

  NumericVector result(G);

  for (int i = 0; i < n_data; i++) {
    double xi = data[i];
    double acc = 0.0;
    double lower = xi - cutoff * h[i];
    double upper = xi + cutoff * h[i];

    for (int j = 0; j < G; j++) {
      double xj = grid[j];

      // If the data is too far away from the evaulated point, skip
      if (xj < lower)
        continue;
      if (xj > upper)
        break;

      double z = (xi - xj) / h[i];
      result[j] += weights[i] * norm_const * std::exp(-0.5 * z * z) /
                   (weight_sum * h[i]);
    }
  }

  return result;
}

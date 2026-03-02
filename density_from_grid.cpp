#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector
density_from_grid_cpp(NumericVector data, double h, NumericVector grid,
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

  double inv_h = 1.0 / h;
  double norm_const = 1.0 / std::sqrt(2.0 * M_PI);

  NumericVector result(G);

  for (int j = 0; j < G; j++) {
    double xj = grid[j];
    double acc = 0.0;

    double lower = xj - cutoff * h;
    double upper = xj + cutoff * h;

    for (int i = 0; i < n_data; i++) {
      // If the data is too far away from the evaulated point, skip
      if (data[i] < lower)
        continue;
      if (data[i] > upper)
        break;

      double z = (data[i] - xj) * inv_h;
      acc += weights[i] * norm_const * std::exp(-0.5 * z * z);
    }

    result[j] = acc / (weight_sum * h);
  }

  return result;
}

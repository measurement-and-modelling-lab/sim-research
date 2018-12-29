#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "getSample.h"

mat getSample(int n, mat intermediate_P, int seed, double a, double b, double c, double d) {
    // n is the sample size
    // P is the population correlation structure
    // a, b, c, and d are Fleishman coefficients

    // Returns n observation data set with correlation structure P
    // and distribution corresponding to a, b, c, and d

    int nvar = intermediate_P.n_cols;

    // Create normal data with intermediate correlation structure
    arma_rng::set_seed(seed);
    vec mu = zeros(nvar);
    mat X = mvnrnd(mu, intermediate_P, n);
    X = X.t();

    // Scale each marginal by the Fleishman coefficients
    for (uword i = 0; i < nvar; i++) {
	vec Xi = X.col(i);
	X.col(i) = a + b * Xi + c * pow(Xi, 2) + d * pow(Xi, 3);
    }

    return X;
}

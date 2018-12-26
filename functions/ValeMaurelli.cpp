#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "ValeMaurelli.h"

double getICOV(double R, double b, double c, double d) {
    // R is a correlation
    // b, c and d are Fleishman coefficients

    // Returns intermediate correlation for R

    double tol = 0.00000000001;
    double increment = 0.5;
    double rho = 1;
    double eq = rho * (pow(b, 2) + 3 * b * d + 3 * d * b + 9 * pow(d, 2)) +
	pow(rho, 2) * (2 * pow(c, 2)) + pow(rho, 3) * (6 * pow(d, 2));

    while (pow(eq-R, 2) > tol) {

	if (eq > R) {
	    rho = rho - increment;
	} else if (eq < R) {
	    increment = increment / 2;
	    rho = rho + increment;
	} else {
	    continue;
	}

	eq = rho * (pow(b, 2) + 3 * b * d + 3 * d * b + 9 * pow(d, 2)) +
	    pow(rho, 2) * (2 * pow(c, 2)) + pow(rho, 3) * (6 * pow(d, 2));
    }

    return rho;
}

mat ValeMaurelli(int n, mat P, double a, double b, double c, double d,
                 int seed) {
    // n is the sample size
    // P is the population correlation structure
    // a, b, c, and d are Fleishman coefficients

    // Returns n observation data set with correlation structure P
    // and distribution corresponding to a, b, c, and d

    int nvar = P.n_cols;

    // Compute intermediate correlation matrix
    mat ICOR = eye(nvar, nvar);
    for (uword j = 0; j < nvar - 1; j++) {
	for (uword i = j + 1; i < nvar; i++) {
	    if (P(j, i) == 0) {
		continue;
	    } else {
		ICOR(j, i) = getICOV(P(j, i), b, c, d);
		ICOR(i, j) = ICOR(j, i);
	    }
	}
    }

    // Create normal data with intermediate correlation structure
    arma_rng::set_seed(seed);
    vec mu = zeros(nvar);
    mat X = mvnrnd(mu, ICOR, n);
    X = X.t();

    // Scale each marginal by the Fleishman coefficients
    for (uword i = 0; i < nvar; i++) {
	vec Xi = X.col(i);
	X.col(i) = a + b * Xi + c * pow(Xi, 2) + d * pow(Xi, 3);
    }

    return X;
}

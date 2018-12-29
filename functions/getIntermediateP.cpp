#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "getIntermediateP.h"


double getIntermediateRho(double R, double b, double c, double d) {
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

mat getIntermediateP(mat P, double b, double c, double d) {

    int nvar = P.n_cols;

    // Compute intermediate correlation matrix
    mat intermediate_P = eye(nvar, nvar);
    for (uword j = 0; j < nvar - 1; j++) {
	for (uword i = j + 1; i < nvar; i++) {
	    if (P(j, i) == 0) {
		continue;
	    } else {
		intermediate_P(j, i) = getIntermediateRho(P(j, i), b, c, d);
		intermediate_P(i, j) = intermediate_P(j, i);
	    }
	}
    }

    return intermediate_P;
}

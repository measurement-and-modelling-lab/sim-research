#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "getIntermediateP.h"


double getIntermediateRho(double R, mat FTable,int i,int j) {
    // R is a correlation
    // b, c and d are Fleishman coefficients

    // Returns intermediate correlation for R
	double b1 = FTable(i,0);
	double c1 = FTable(i,1);
	double d1 = FTable(i,2);
	double b2 = FTable(j,0);
	double c2 = FTable(j,1);
	double d2 = FTable(j,2);

    double tol = 0.00000000001;
    double increment = 0.5;
    double rho = 1;
    double eq = rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
	rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2);

    while (pow(eq-R, 2) > tol) {

		if (eq > R) {
			rho = rho - increment;
		} else if (eq < R) {
			increment = increment / 2;
			rho = rho + increment;
		} else {
			continue;
		}

		eq = rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
		rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2);
    }

    return rho;
}

mat getIntermediateP(mat P, mat FTable) {

    int nvar = P.n_cols;

    // Compute intermediate correlation matrix
    mat intermediate_P = eye(nvar, nvar);
    for (uword j = 0; j < nvar - 1; j++) {
		for (uword i = j + 1; i < nvar; i++) {
			if (P(j, i) == 0) {
				continue;
			} else {
				intermediate_P(j, i) = getIntermediateRho(P(j, i),FTable,i,j);
				intermediate_P(i, j) = intermediate_P(j, i);
			}
		}
    }

    return intermediate_P;
}

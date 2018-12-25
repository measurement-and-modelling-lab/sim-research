#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "getICOV.h"

double getICOV(double R, double b, double c, double d) {
    // R is a correlation
    // b, c and d are Fleishman coefficients

    // Returns intermediate correlation for R

    double tol = 0.001;
    double increment = 0.01;
    double rho = R + 0.1;
    double eq = rho * (pow(b, 2) + 3 * b * d + 3 * d * b + 9 * pow(d, 2)) +
	pow(rho, 2) * (2 * pow(c, 2)) + pow(rho, 3) * (6 * pow(d, 2));

    while (eq - R > tol) {

	if (eq - R > tol) {
	    increment = increment / 1.01;
	    rho = rho - increment;
	} else if (rho - R < tol) {
	    rho = rho + increment;
	} else {
	    continue;
	}

	eq = rho * (pow(b, 2) + 3 * b * d + 3 * d * b + 9 * pow(d, 2)) +
	    pow(rho, 2) * (2 * pow(c, 2)) + pow(rho, 3) * (6 * pow(d, 2));
    }

    return rho;
}

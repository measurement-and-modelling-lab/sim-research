#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "counsell2015.h"

double counsell2015(mat R, int n, double delta) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the minimum difference between r12 and r13 under the null hypothesis
    // Returns a p value for equivalence test of r12 and r13

    double r12 = R(0, 1);
    double r13 = R(0, 2);
    double r23 = R(1, 2);

    double detR =
	(1 - pow(r12, 2) - pow(r13, 2) - pow(r23, 2)) + (2 * r12 * r13 * r23);
    double s = sqrt(((n - 1) * (1 + r23)) /
		    ((2 * ((n - 1) / (n - 3)) * detR) +
		     ((pow((r12 + r13), 2)) / 4) * (pow((1 - r23), 3))));

    double z1 = (fabs(r12 - r13) - delta) * s;
    double z2 = (-fabs(r12 - r13) - delta) * s;
    double p = normcdf(z1) - normcdf(z2);

    return p;
}

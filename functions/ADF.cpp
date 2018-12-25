#define ARMA_NO_DEBUG
#include <armadillo>
using namespace arma;
#include "ADF.h"

double ADF(mat R, int n, double delta, vec moments) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the difference between r12 and r13 under H0
    // moments is the 4th order moments of the data set

    // Returns p value for test of |r12-r13| <= delta

    double r12 = R(0, 1);
    double r13 = R(0, 2);
    double r23 = R(1, 2);

    // Z ADF
    double gamma_12 = adfCov(0, 1, 0, 1, R, moments);
    double gamma_13 = adfCov(0, 2, 0, 2, R, moments);
    double gamma_12_13 = adfCov(0, 1, 0, 2, R, moments);

    double z1 = sqrt(n) * (fabs(r12 - r13) - delta) *
	(1 / sqrt(gamma_12 + gamma_13 - 2 * gamma_12_13));
    double z2 = sqrt(n) * (-fabs(r12 - r13) - delta) *
	(1 / sqrt(gamma_12 + gamma_13 - 2 * gamma_12_13));
    double p = normcdf(z1) - normcdf(z2);

    return p;
}

double adfCov(int i, int j, int k, int h, mat R, vec moments) {
    // i/j and k/h are the row and column indices of two correlations
    // R is the correlation matrix
    // moments is the 4th order moments for the entire data set

    // Returns asymptotic covariance of rij and rkh

    double term1 = FRHO(i, j, k, h, moments);
    double term2 = 0.25 * R(i, j) * R(k, h) *
	(FRHO(i, i, k, k, moments) + FRHO(j, j, k, k, moments) +
	 FRHO(i, i, h, h, moments) + FRHO(j, j, h, h, moments));
    double term3 =
	0.5 * R(i, j) *
	(FRHO(i, i, k, h, moments) + FRHO(j, j, k, h, moments));
    double term4 =
	0.5 * R(k, h) *
	(FRHO(i, j, k, k, moments) + FRHO(i, j, h, h, moments));
    double cov = term1 + term2 - term3 - term4;

    return cov;
}

int findpos(int i, int j, int k, int h) {
    // i/j and k/h are the row and column indices of two correlations

    // Returns index for kurtosis of rij and rkh in moments vector

    vec indices(4);
    indices(0) = i;
    indices(1) = j;
    indices(2) = k;
    indices(3) = h;

    indices = sort(indices, "descend");

    int a = indices(0) + 1;
    int b = indices(1) + 1;
    int c = indices(2) + 1;
    int d = indices(3) + 1;

    int index = (a - 1) * a * (a + 1) * (a + 2) / 24 +
	        (b - 1) * b * (b + 1) / 6 +
	         c * (c - 1) / 2 +
	         d;

    return index - 1;
}

double FRHO(int i, int j, int k, int h, vec moments) {
    // i/j and k/h are the row and column indices of two correlations
    // moments is the fourth order moments of the data set

    // Returns moments for rij and rkh

    int temp = findpos(i, j, k, h);
    double fpho = moments(temp);

    return fpho;
}

#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "vale-maurelli.h"

double kurtosis(vec X) {
    // X is a vector of scores
    // Returns kurtosis of vector

    double m4 = mean(pow(X - mean(X), 4));
    double kurt = m4 / pow(stddev(X), 4) - 3;
    return kurt;
}

double skewness(vec X) {
    // X is a vector of scores
    // Returns skewness of vector

    double m3 = mean(pow(X - mean(X), 3));
    double skew = m3 / pow(stddev(X), 3);
    return skew;
}

double ksD(vec X) {
    // X is a set of p values
    // Returns uniform Kolmogorov-Smirnov D

    vec expected = sort(X);
    int n = expected.n_rows;
    vec observed = linspace(1, n, n) / n;
    vec d = abs(expected - observed);
    return max(d);
}

vec compute4thOrderMoments(mat X) {
    // X is a matrix of data
    // Returns fourth order moments of X

    int n = X.n_rows;
    int p = X.n_cols;
    rowvec means = mean(X, 0);
    rowvec sds = stddev(X, 0);
    mat zscores = X;

    for (uword i = 0; i < p; i++) {
	zscores.col(i) = zscores.col(i) - means(i);
	zscores.col(i) = zscores.col(i) / sds(i);
    }

    int q = p * (p + 1) * (p + 2) * (p + 3) / 24;
    vec moments(q);
    int a = 0;

    for (uword i = 0; i < p; i++) {
	for (uword j = 0; j <= i; j++) {
	    for (uword k = 0; k <= j; k++) {
		for (uword h = 0; h <= k; h++) {
		    for (uword b = 0; b < n; b++) {
			moments(a) = moments(a) + zscores(b, i) * zscores(b, j) *
			    zscores(b, k) * zscores(b, h);
		    }
		    a++;
		}
	    }
	}
    }

    moments = moments / (n - 1);

    return moments;
}

double ADF(mat R, int n, double delta, vec moments) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the difference between r12 and r13 under H0
    // moments is the 4th order moments of the data set
    // Returns p value for test of |r12-r12| <= delta

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
	0.5 * R(i, j) * (FRHO(i, i, k, h, moments) + FRHO(j, j, k, h, moments));
    double term4 =
	0.5 * R(k, h) * (FRHO(i, j, k, k, moments) + FRHO(i, j, h, h, moments));
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

double counsell(mat R, int n, double delta) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the difference between r12 and r13 under H0
    // Returns p value for test of |r12-r12| <= delta

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

double fisher(double r) {
    // r is a correlation coefficient
    // Returns Fisher transform of r

    double z = 0.5 * (log(1 + r) - log(1 - r));
    return z;
}

double getICOV(double R, double b, double c, double d) {
    // R is a correlation
    // b, c and d are Fleishman coefficients
    // Returns intermediary correlation for R

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

mat ValeMaurelli(int n, mat P, double a, double b, double c, double d,
		 int seed) {
    // n is the sample size
    // P is the population correlation structure
    // a, b, c, and d are Fleishman coefficients
    // Returns n observation data set with correlation structure P
    //  and distribution corresponding to a, b, c, and d

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

    // Create normal data with specified correlation structure
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

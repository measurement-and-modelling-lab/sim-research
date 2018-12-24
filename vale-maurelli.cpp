#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "vale-maurelli.h"

double ksD(vec data) { // Uniform data only
    data = sort(data);
    int n = data.n_rows;
    vec t1error = linspace<vec>(1, n, n) / n;
    vec d = abs(t1error - data);
    return max(d);
}

vec compute4thOrderMoments(mat data) {

    int n = data.n_rows;
    int p = data.n_cols;
    rowvec means = mean(data, 0);
    rowvec sds = stddev(data, 0);
    mat zscores = data;

    for (uword i = 0; i < p; i++) {
	zscores.col(i) = zscores.col(i) - means(i);
	zscores.col(i) = zscores.col(i) / sds(i);
    }

    int q = p * (p + 1) * (p + 2) * (p + 3) / 24;
    vec moments = zeros(q);
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

double ADF(double r12, double r13, double r23, mat R, int n, mat sample,
           double delta, vec moments) {

    // Z ADF
    double gamma_12 = adfCov(0, 1, 0, 1, R, moments);
    double gamma_13 = adfCov(0, 2, 0, 2, R, moments);
    double gamma_12_13 = adfCov(0, 1, 0, 2, R, moments);
    double z1 = sqrt(n) * (fabs(r12 - r13) - delta) *
	(1 / sqrt(gamma_12 + gamma_13 -
		  2 * gamma_12_13));
    double z2 = sqrt(n) * (-fabs(r12 - r13) - delta) *
	(1 / sqrt(gamma_12 + gamma_13 -
		  2 * gamma_12_13));
    double p = normcdf(z1) - normcdf(z2);
    return p;

}

double adfCov(int i, int j, int k, int h, mat R, vec moments) {
    // Called by ComputeWBCorrChiSquare.R
    // a, i and j are the group, row and column number of one correlation from the
    // hypothesis matrix b, k and h are the group, row and column number of
    // another (possibly the same) correlation from the hypothesis matrix R is the
    // list of correlation matrices if using ADF, and the list of OLS matrices if
    // using TSADF moments are the fourth order moments, i.e. values on kurtosis
    // output is the covariance of the two correlations

    double term1 = FRHO(i, j, k, h, moments);
    double term2 = 0.25 * R[i, j] * R[k, h] *
	(FRHO(i, i, k, k, moments) + FRHO(j, j, k, k, moments) +
	 FRHO(i, i, h, h, moments) + FRHO(j, j, h, h, moments));
    double term3 =
	0.5 * R[i, j] * (FRHO(i, i, k, h, moments) + FRHO(j, j, k, h, moments));
    double term4 =
	0.5 * R[k, h] * (FRHO(i, j, k, k, moments) + FRHO(i, j, h, h, moments));
    double cov = term1 + term2 - term3 - term4;

    return cov;
}

int findpos(int i, int j, int k, int h) {
    // Called by FRHO.R
    // i and j are the row and column number of one correlation from the
    // hypothesis matrix k and h are the row and column number of another
    // (possibly the same) correlation from the hypothesis matrix Returns the
    // index number of the kurtosis for the two correlations in the moments vector
    // (M)

    rowvec indices(4);
    indices(0) = i;
    indices(1) = j;
    indices(2) = k;
    indices(3) = h;

    indices = sort(indices, "descend");

    int a = indices[0] + 1;
    int b = indices[1] + 1;
    int c = indices[2] + 1;
    int d = indices[3] + 1;

    int index = (a - 1) * a * (a + 1) * (a + 2) / 24 +
	        (b - 1) * b * (b + 1) / 6 +
	         c * (c - 1) / 2 +
	         d;

    return index - 1;
}

double FRHO(int i, int j, int k, int h, vec M) {
    // Called by adfCov.R
    // j and k are the row and column number of one correlation from the
    // hypothesis matrix h and m are the row and column number of another
    // (possibly the same) correlation from the hypothesis matrix Returns the
    // kurtosis for two the correlations

    int temp = findpos(i, j, k, h);
    double fpho = M[temp];
    return fpho;
}

double counsell(double r12, double r13, double r23, int n, double delta) {
    double detR =
	(1 - pow(r12, 2) - pow(r13, 2) - pow(r23, 2)) + (2 * r12 * r13 * r23);
    double s = sqrt(((n - 1) * (1 + r23)) /
		    ((2 * ((n - 1) / (n - 3)) * detR) +
		     ((pow((r12 + r13), 2)) / 4) * (pow((1 - r23), 3))));
    double p1 = normcdf((fabs(r12 - r13) - delta) * s);
    double p2 = normcdf((-fabs(r12 - r13) - delta) * s);
    return (p1 - p2);
}

double fisher(double r) {
    double z = 0.5 * (log(1 + r) - log(1 - r));
    return z;
}

mat cov2cor(mat S) {
    vec d = S.diag();
    vec d_squared = square(d);
    mat R(3, 3);
    R.zeros();
    for (uword i = 0; i < R.n_cols; i++) {
	for (uword j = 0; j < R.n_cols; j++) {
	    R(j, i) = S(j, i) / (d_squared(i) * d_squared(j));
	}
    }
    return R;
}

double getICOV(double R, double b1, double c1, double d1, double b2, double c2,
               double d2) {

    double tol = 0.001;
    double increment = 0.01;
    double rho = R + 0.1;
    double eq = rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
	rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2);

    while (eq - R > tol) {

	if (eq - R > tol) {
	    increment = increment / 1.01;
	    rho = rho - increment;
	} else if (rho - R < tol) {
	    rho = rho + increment;
	} else {
	    continue;
	}

	eq = rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
	    rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2);
    }
    return rho;
}

mat ValeMaurelli1983(int n, mat COR, double a, double b, double c, double d,
                     int seed) {

    uword nvar = COR.n_cols;

    // Create table of Fleishman coefficients
    mat FTable;
    FTable.zeros(nvar, 4);
    rowvec values(4);
    values(0) = a;
    values(1) = b;
    values(2) = c;
    values(3) = d;
    for (uword i = 0; i < nvar; i++) {
	FTable.row(i) = values;
    }
    // Compute intermediate correlation matrix
    mat ICOR = eye<mat>(nvar, nvar);
    for (uword j = 0; j < nvar - 1; j++) {
	for (uword i = j + 1; i < nvar; i++) {
	    if (COR(j, i) == 0) {
		continue;
	    } else {
		ICOR(j, i) =
		    getICOV(COR(j, i), FTable(j, 1), FTable(j, 2), FTable(j, 3),
			    FTable(i, 1), FTable(i, 2), FTable(i, 3));
		ICOR(i, j) = ICOR(j, i);
	    }
	}
    }

    arma_rng::set_seed(seed);
    vec mu = zeros<vec>(nvar);
    mat Z = mvnrnd(mu, ICOR, n);
    mat X;
    X = Z.t();
    Z = Z.t();
    for (uword i = 0; i < nvar; i++) {
	vec Zi = Z.col(i);
	X.col(i) = FTable(i, 0) + FTable(i, 1) * Zi + FTable(i, 2) * square(Zi) +
	    FTable(i, 3) * pow(Zi, 3);
    }

    return X;
}

mat mvrnonnorm(int n, mat Sigma, double a, double b, double c,
               double d, int seed) {

    uword nvar = Sigma.n_cols;
    mat Z = ValeMaurelli1983(n, cov2cor(Sigma), a, b, c, d, seed);

    // Divide each column by the corresponding diagonal element of Sigma
    for (uword i = 0; i < nvar; i++) {
	double element = Sigma(i, i);
	Z.col(i) = Z.col(i) / linspace<vec>(element, element, n);
    }

    return Z;
}

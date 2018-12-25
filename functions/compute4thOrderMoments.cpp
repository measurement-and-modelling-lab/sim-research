#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "compute4thOrderMoments.h"

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

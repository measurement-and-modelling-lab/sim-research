#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "kolmogorovD.h"

double kolmogorovD(vec X) {
    // X is a set of p values

    // Returns uniform Kolmogorov-Smirnov D

    vec expected = sort(X);
    int n = expected.n_rows;
    vec observed = linspace(1, n, n) / n;
    vec d = abs(expected - observed);
    return max(d);
}

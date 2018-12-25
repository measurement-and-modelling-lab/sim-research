#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "skewness.h"

double skewness(vec X) {
    // X is a vector of scores

    // Returns skewness of vector

    double m3 = mean(pow(X - mean(X), 3));
    double skew = m3 / pow(stddev(X), 3);
    return skew;
}

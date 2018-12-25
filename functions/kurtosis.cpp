#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
using namespace arma;
#include "kurtosis.h"

double kurtosis(vec X) {
    // X is a vector of scores

    // Returns kurtosis of vector

    double m4 = mean(pow(X - mean(X), 4));
    double kurt = m4 / pow(stddev(X), 4) - 3;
    return kurt;
}

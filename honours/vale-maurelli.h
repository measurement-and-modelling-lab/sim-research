#include <armadillo>
using namespace arma;

double counsell(double r12, double r13, double r23, int n, double delta);

double fisher(double r);

mat cov2cor(mat S);

// Retrieve coefficients for skew/kurtosis pair
rowvec fleishman1978(double skewness, double kurtosis);

double getICOV(double R, double b1, double c1, double d1, double b2, double c2,
               double d2);

mat ValeMaurelli1983(int n, mat COR, double a, double b,double c, double d, int seed);

mat mvrnonnorm(int n, double mu, mat Sigma, double a, double b,double c, double d, int seed);

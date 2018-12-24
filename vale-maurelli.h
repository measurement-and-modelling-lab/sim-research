#include <armadillo>
using namespace arma;

double kurtosis(vec X);

double skewness(vec X);

double ksD(vec X);

vec compute4thOrderMoments(mat X);

double ADF(mat R, int n, double delta, vec moments);

double adfCov(int i, int j, int k, int h, mat R, vec moments);

int findpos (int i, int j, int k, int h);

double FRHO (int i, int j, int k, int h, vec M);

double counsell(mat R, int n, double delta);

double fisher(double r);

mat cov2cor(mat S);

double getICOV(double R, double b, double c, double d);

mat ValeMaurelli(int n, mat P, double a, double b,double c, double d, int seed);

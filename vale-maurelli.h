#include <armadillo>
using namespace arma;

double ksD(vec data);

vec compute4thOrderMoments(mat data);

double ADF(double r12, double r13, double r23, mat R, int n, mat sample, double delta,vec moments);

double adfCov(int i, int j, int k, int h, mat R, vec moments);

int findpos (int i, int j, int k, int h);

double FRHO (int i, int j, int k, int h, vec M);

double counsell(double r12, double r13, double r23, int n, double delta);

double fisher(double r);

mat cov2cor(mat S);

double getICOV(double R, double b, double c, double d);

mat ValeMaurelli(int n, mat P, double a, double b,double c, double d,int seed);

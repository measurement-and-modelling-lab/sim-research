#define ARMA_NO_DEBUG // disable bound checks to improve speed
#include <armadillo>
#include <ctime>
using namespace arma;
#include "vale-maurelli.h"

int main(int argc, char const *argv[]) {

  int i, n;
  double rho12, rho13, rho23, a, b, c, d;

  i = atof(argv[1]);
  n = atof(argv[2]);
  rho12 = atof(argv[3]);
  rho23 = atof(argv[4]);
  a = atof(argv[5]);
  b = atof(argv[6]);
  c = atof(argv[7]);
  d = atof(argv[8]);

  double mu = 0;
  double delta = 0.1;

  mat S;
  S.ones(3, 3);
  S(0, 1) = rho12;
  S(0, 2) = rho13;
  S(1, 2) = rho23;
  S = symmatu(S);

  int seed = 12345;

  for (int k = 0; k < 1000; k++) {

    mat sample = mvrnonnorm(n, mu, S, a, b, c, d, seed);

    mat R = cor(sample);
    double r12 = R(0, 1);
    double r13 = R(0, 2);
    double r23 = R(1, 2);

    cout << i
	 << ","
	 << counsell(r12, r13, r23, n, delta)
	 << ","
	 << seed
	 << endl;
  }

  return 0;
}

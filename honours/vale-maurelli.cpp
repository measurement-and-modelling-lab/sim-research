#define ARMA_NO_DEBUG // disable bound checks to improve speed
#include <armadillo>
#include <ctime>
using namespace arma;
#include "vale-maurelli.h"

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

// Retrieve coefficients for skew/kurtosis pair
rowvec fleishman1978(double skewness, double kurtosis) {
  static mat fleishmanTable;
  fleishmanTable.load("coefficients.csv");
  int index = -1;
  for (uword i = 0; i < fleishmanTable.n_rows; i++) {
    if ((fleishmanTable(i, 0) == skewness) &&
        (fleishmanTable(i, 1) == kurtosis)) {
      index = i;
    }
  }

  if (index == -1) {
    return NULL;
  } else {
    double a, b, c, d;
    b = fleishmanTable(index, 2);
    c = fleishmanTable(index, 3);
    d = fleishmanTable(index, 4);
    a = -c;
    rowvec values(4);
    values(0) = a;
    values(1) = b;
    values(2) = c;
    values(3) = d;
    return values;
  }
}

double getICOV(double R, double b1, double c1, double d1, double b2, double c2,
               double d2) {

  double tol = 0.00000001;
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

mat ValeMaurelli1983(int n, mat COR, double a, double b, double c, double d) {

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

  arma_rng::set_seed(1234567); // Should seed with our RNG
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

mat mvrnonnorm(int n, double mu, mat Sigma, double a, double b, double c,
               double d) {

  uword nvar = Sigma.n_cols;
  mat Z = ValeMaurelli1983(n, cov2cor(Sigma), a, b, c, d);
  // Divide each column by the corresponding diagonal element of Sigma
  for (uword i = 0; i < nvar; i++) {
    double element = Sigma(i, i);
    Z.col(i) = Z.col(i) / linspace<vec>(element, element, n);
  }

  // Add mu to every score
  if (mu != 0) {
    Z.for_each([&mu](mat::elem_type &val) { val += mu; });
  }

  return Z;
}

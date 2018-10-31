#include <armadillo>
using namespace arma;

mat cov2cor(mat S) {

  vec d = S.diag();
  vec d_squared = d.transform([](int val) { return (val * val); });
  mat R = S.zeros();
  for (int i = 0; i < R.n_cols; i++) {
    for (int j = 0; j < R.n_cols; j++) {
      R(i, j) = S(i, j) / (d_squared(i) * d_squared(j));
    }
  }
  return R;
}

// Retrieve coefficients for skew/kurtosis pair
rowvec fleishman1978(double skewness, double kurtosis) {
  static mat fleishmanTable;
  fleishmanTable.load("coefficients.csv");
  int index = -1;
  for (int i = 0; i < fleishmanTable.n_rows; i++) {
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
    a = -1 * c;
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
      increment = increment / 1.1;
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

mat ValeMaurelli1983(int n, mat COR, double skewness, double kurtosis) {

  uword nvar = COR.n_cols;

  // Create table of Fleishman coefficients
  mat FTable;
  FTable.zeros(nvar, 4);
  for (int i = 0; i < nvar; i++) {
    FTable.row(i) = fleishman1978(skewness, kurtosis);
  }

  // Compute intermediate correlations between all pairs
  mat ICOR = eye<mat>(nvar, nvar);
  for (int j = 0; j < nvar - 1; j++) {
    for (int i = j + 1; i < nvar - 1; i++) {
      if (COR(i, j) == 0) {
        continue;
      } else {
        ICOR(i, j) =
            getICOV(FTable(i, 1), FTable(i, 2), FTable(i, 3), FTable(j, 1),
                    FTable(j, 2), FTable(j, 3), COR(i, j));
        ICOR(j, i) = ICOR(i, j);
      }
    }
  }

  arma_rng::set_seed(1234567); // Seed with our RNG
  vec mu = zeros<vec>(nvar);
  mat Z = mvnrnd(mu, ICOR, n);

  mat X;
  X = Z.t();
  Z = Z.t();
  for (int i = 0; i < nvar; i++) {
    vec Z1 = Z.col(i);
    vec Z2 = Z1.transform([](double val) { return (val * val); });       // Zi^2
    vec Z3 = Z1.transform([](double val) { return (val * val * val); }); // Zi^3
    X.col(i) = FTable(i, 0) + FTable(i, 1) * Z1 + FTable(i, 2) * Z2 +
               FTable(i, 3) * Z3;
  }

  return X;
}

mat mvrnonnorm(int n, double mu, mat Sigma, double skewness, double kurtosis) {

  int nvar = Sigma.n_cols;
  mat Z = ValeMaurelli1983(n, cov2cor(Sigma), skewness, kurtosis);

  // Divide each column by the corresponding diagonal element of Sigma
  for (int i = 0; i < nvar; i++) {
    double element = Sigma(i, i);
    Z.col(i) = Z.col(i) / linspace<vec>(element, element, n);
  }

  // Add mu to every score
  Z.for_each([&mu](mat::elem_type &val) { val += mu; });

  return Z;
}

int main(int argc, char const *argv[]) {
  // test values
  mat COV;
  COV.eye(5, 5);
  int n = 1000;
  double mu = 0;
  double skewness = 1.75;
  double kurtosis = 3.75;
  mat result = mvrnonnorm(n, mu, COV, skewness, kurtosis);
  cout << cor(result) << endl;

  return 0;
}

#include <armadillo>
using namespace arma;

vec compute4thOrderMoments(mat data) {

  int n = data.n_rows;
  int p = data.n_cols;
  rowvec means = mean(data, 0);
  rowvec sds = stddev(data, 0);
  mat zscores = data;

  for (uword i = 0; i < p; i++) {
    zscores.col(i) = zscores.col(i) - means(i);
    zscores.col(i) = zscores.col(i) / sds(i);
  }

  int q = p * (p + 1) * (p + 2) * (p + 3) / 24;
  vec moments = zeros(q);
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

int main() {

  mat I = eye(5, 5);
  colvec M = zeros(5);
  mat data = mvnrnd(M, I, 50).t();

  cout << compute4thOrderMoments(data) << endl;

  return 0;
}

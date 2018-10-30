#include <iostream>

using namespace std;

float cov2cor (mat S) {

  d = S.diag
  d.squared = d.transform( [](int val) { return (val * val); } );
  R = S.zeros()
    for (int i=0; i<R.n_cols; i++) {
      for (int j=0; j<R.n_cols; j++) {
	R(i,j) = S(i,j) / (d(i) * d(j))
      }
    }

  return(R)

}


int main () {


  cout << R << endl;

}

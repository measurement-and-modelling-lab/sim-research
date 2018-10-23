#include <iostream>

using namespace std;

int main () {

  float b1 = 0.92966052480111;
  float c1 = 0.39949667453766;
  float d1 = -0.03646699281275;
  float b2 = 0.86588620352314;
  float c2 = 0.22102762101262;
  float d2 = 0.02722069915809;
  float R = 0.9;
  float rho = R + 0.1;

  float tol = 0.00000001;
  float increment = 0.01;

  float estimate = rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + rho * rho *(2*c1*c2) + rho * rho * rho *(6*d1*d2);

  while (estimate - R > tol) {
  
    if (estimate - R > tol) {
      increment = increment / 1.1;
      rho = rho - increment;
    } else if (estimate- R < tol) {
      rho = rho + increment;
    } else {
      continue;
    }

    estimate = rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + rho * rho *(2*c1*c2) + rho * rho * rho *(6*d1*d2);

  }

  cout << rho << endl;

}

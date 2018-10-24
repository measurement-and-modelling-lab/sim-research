#include <iostream>

using namespace std;


float getICOV (float R, float b1, float c1, float d1, float b2, float c2, float d2) {

  float tol = 0.00000001;
  float increment = 0.01;
  float rho = R + 0.1;
  float eq = rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + rho * rho *(2*c1*c2) + rho * rho * rho *(6*d1*d2);

  while (eq - R > tol) {
  
    if (eq - R > tol) {
      increment = increment / 1.1;
      rho = rho - increment;
    } else if (rho - R < tol) {
      rho = rho + increment;
    } else {
      continue;
    }

    eq = rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + rho * rho *(2*c1*c2) + rho * rho * rho *(6*d1*d2);

  }

  return rho;
}


int main () {

  float b1 = 0.92966052480111;
  float c1 = 0.39949667453766;
  float d1 = -0.03646699281275;
  float b2 = 0.86588620352314;
  float c2 = 0.22102762101262;
  float d2 = 0.02722069915809;
  float R = 0.9;

  float rho = getICOV(R, b1, c1, d1, b2, c2, d2);

  cout << rho << endl;

}

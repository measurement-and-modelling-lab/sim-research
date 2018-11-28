#define ARMA_NO_DEBUG // disable bound checks to improve speed
#include <armadillo>
#include <ctime>
using namespace arma;
#include "vale-maurelli.h"

int main(int argc, char const *argv[]) {
  int interations = 1000;
  mat results(2*interations,3);
  mat scenarios;
  scenarios.load("conditions.csv");
  vec seeds;
  seeds.load("seeds.csv");
  for (int i = 0; i < scenarios.n_rows; i++) {
  	vec ADFValues(interations);
  	vec counsellValues(interations);
  	int adfCounter =0 ;
  	int counsellCounter = 0;
    for (int k = 0; k < interations; k++) {
      int seedIndex = (i) * interations + k;
      double n = scenarios(i, 2);
      double rho12 = scenarios(i, 3);
      double rho13 = scenarios(i, 4);
      double rho23 = rho13 + scenarios(i, 5);
      double a = scenarios(i, 6);
      double b = scenarios(i, 7);
      double c = scenarios(i, 8);
      double d = scenarios(i, 9);
      double mu = 0;

      mat S;
      S.ones(3, 3);
      S(0, 1) = rho12;
      S(0, 2) = rho13;
      S(1, 2) = rho23;
      S = symmatu(S);

      mat sample = mvrnonnorm(n, mu, S, a, b, c, d,seeds(seedIndex));

      mat R = cor(sample);
      double r12 = R(0, 1);
      double r13 = R(0, 2);
      double r23 = R(1, 2);
      double delta = 0.1;
      vec moments = compute4thOrderMoments(sample);

	  ADFValues(k) = ADF(r12,r13,r23,R,n,sample,delta,moments);
	  if(ADFValues(k) <= .05){
	  	adfCounter++;
	  }
      counsellValues(k) = counsell(r12, r13, r23, n, delta); 
      if( counsellValues(k) <= .05){
	  	counsellCounter++;
	  }

    }
   	double proportionAdf = (double)adfCounter/(double)interations;
    double proportionCounsell = (double)counsellCounter/(double)interations;
    double adfD =  kst(ADFValues);
    double counsellD =  kst(counsellValues);

    // results(2i,0) = "ADF";
    // results(2i,1) = proportionAdf;
    // //results(2i,2) =
    // results(2i+1,0) = "counsell";
    // results(2i+1,1) = proportionCounsell;
    // //results(2i+1,2) =
    cout << "ADF,";
    cout<<proportionAdf;
    cout << ",";
    cout << adfD <<endl;
    cout << "Counsell,";
    cout<<proportionCounsell;
    cout << ",";
    cout << counsellD << endl;



  }


  return 0;
}

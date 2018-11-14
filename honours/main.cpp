#define ARMA_NO_DEBUG // disable bound checks to improve speed
#include <armadillo>
#include <ctime>
using namespace arma;
#include "vale-maurelli.h"

int main(int argc, char const *argv[])
{

	mat scenarios;
	scenarios.load("scenario1.csv");
	for(int i = 1; i<scenarios.n_rows; i++){
		for(int k = 0; k< 1000; k++){
			double a = scenarios(i,6);
			double b = scenarios(i,7);
			double c = scenarios(i,8);
			double d = scenarios(i,9);
			int n = 50;
			double rho12 = scenarios(i,3);
			double rho13 =  scenarios(i,4);
			double rho23 = rho13 + scenarios(i,5);
			double mu =0;
			mat S; 
			S.ones(3,3);
			S(0,1) = rho12;
			S(0,2) = rho13;
			S(1,2) = rho23;
			S = symmatu(S);
			
			mat sample = mvrnonnorm(n,mu,S,a,b,c,d);
			mat R = cor(sample);
			double r12 = R(0,1);
			double r13 = R(0,2);
			double r23 = R(1,2);
			double delta = 0.1;
			double p = counsell(r12,r13,r23,n,delta);
		}
		cout << i << endl;
	}
	return 0;
}
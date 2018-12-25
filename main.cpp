#include <armadillo>
using namespace arma;
#define ARMA_NO_DEBUG // disable bound checks to improve speed

#include "ADF.h"
#include "compute4thOrderMoments.h"
#include "counsell.h"
#include "ksD.h"
#include "kurtosis.h"
#include "skewness.h"
#include "ValeMaurelli.h"

int main(void) {

    // Iterations per condition
    int iterations = 1000;

    mat conditions;
    conditions.load("conditions.csv");

    vec seeds;
    seeds.load("seeds.csv");

    // Header for .csv output
    cout << "condition,statistic,alpha,D" << endl;

    double delta = 0.1;

    for (int i = 0; i < conditions.n_rows; i++) {

    	// The p values for each statistic
    	vec p_adf(iterations);
    	vec p_counsell(iterations);

    	// Count the number of rejections
    	int counter_adf = 0;
    	int counter_counsell = 0;

    	int n = conditions(i, 2);
    	double rho12 = conditions(i, 3);
    	double rho23 = conditions(i, 4);
    	double b = conditions(i, 5);
    	double c = conditions(i, 6);
    	double d = conditions(i, 7);
    	double a = -c;

    	// Population correlation matrix
    	mat P;
    	P.ones(3, 3);
    	P(0, 1) = rho12;
    	P(0, 2) = rho12 - 0.1;
    	P(1, 2) = rho23;
    	P = symmatu(P);

    	for (int j = 0; j < iterations; j++) {

    	    int seed_index = i * iterations + j;
    	    mat sample = ValeMaurelli(n, P, a, b, c, d, seeds(seed_index));

    	    mat R = cor(sample);
	    vec moments = compute4thOrderMoments(sample);

    	    p_adf(j) = ADF(R, n, delta, moments);
    	    if (p_adf(j) <= .05) {
    		counter_adf++;
    	    }

    	    p_counsell(j) = counsell(R, n, delta);
    	    if (p_counsell(j) <= .05) {
    		counter_counsell++;
    	    }
    	}

    	double error_adf = (double)counter_adf / (double)iterations;
    	double D_adf =  ksD(p_adf);

    	double error_counsell = (double)counter_counsell / (double)iterations;
    	double D_counsell = ksD(p_counsell);

    	cout << i << "," << "ADF" << "," << error_adf << "," << D_adf << endl;
    	cout << i << "," << "Counsell" << "," << error_counsell << "," << D_counsell << endl;

    }

    return 0;
}

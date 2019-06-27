#include <armadillo>
using namespace arma;
#define ARMA_NO_DEBUG // disable bound checks to improve speed

#include "serafini2019.h"
#include "compute4thOrderMoments.h"
#include "counsell2015.h"
#include "kolmogorovD.h"
#include "kurtosis.h"
#include "skewness.h"
#include "getIntermediateP.h"
#include "getSample.h"
#include <iomanip> 



rowvec getCoeffs(mat coefficients, double skewness, double kurtosis){
	rowvec values(3);
	for(int i = 0; i < coefficients.n_rows;i++){
		if(coefficients(i,0) == skewness && coefficients(i,1) == kurtosis ){
			values(0) = coefficients(i,2);
			values(1) = coefficients(i,3);
			values(2) = coefficients(i,4);
		}
	}
	return values;
}

int main(void) {

    // The number of iterations per condition
    int iterations = 12000;

    mat conditions;
    conditions.load("conditions.csv");
	mat coefficients;
	coefficients.load("coefficients.csv");

	// Only shed row if the csv has a header row
	coefficients.shed_row(0);

    vec seeds;
    seeds.load("seeds.csv");

    // Print header for .csv
    cout << "skewness_nominal1,kurtosis_nominal1,skewness_nominal2,kurtosis_nominal2,skewness_nominal3,kurtosis_nominal3,rho12,rho23,statistic,alpha,D,skewness_observed,kurtosis_observed" << endl;

    double delta = 0.1;

    for (int i = 1; i < conditions.n_rows; i++) {

    	// The p values for each statistic
    	vec p_serafini(iterations);
    	vec p_counsell(iterations);

	// The sum of the skewness and kurtosis for each marginal of each sample
    	double skewness_observed = 0;
    	double kurtosis_observed = 0;

    	// Counters for the number of rejections
    	int counter_serafini = 0;
    	int counter_counsell = 0;

	// Extract condition parameters
		int skewness_nominal1 = conditions(i, 0);
		int kurtosis_nominal1 = conditions(i, 1);
		int skewness_nominal2 = conditions(i, 2);
		int kurtosis_nominal2 = conditions(i, 3);
		int skewness_nominal3 = conditions(i, 4);
		int kurtosis_nominal3 = conditions(i, 5);
    	int n = conditions(i, 6);
    	double rho12 = conditions(i, 7);
    	double rho23 = conditions(i, 8);
		mat FTable;
		FTable.insert_rows(0,getCoeffs(coefficients,skewness_nominal1,kurtosis_nominal1));
		FTable.insert_rows(1,getCoeffs(coefficients,skewness_nominal2,kurtosis_nominal2));
		FTable.insert_rows(2,getCoeffs(coefficients,skewness_nominal3,kurtosis_nominal3));
		cout.precision(6);

    	// Assemble the population correlation matrix
    	mat P;
    	P.ones(3, 3);
    	P(0, 1) = rho12;
    	P(0, 2) = rho12 - 0.1;
    	P(1, 2) = rho23;
    	P = symmatu(P);


	// Calculate the intermediate correlation matrix (Value & Maurelli, 1983)
	mat intermediate_P = getIntermediateP(P, FTable);
    	for (int j = 0; j < iterations; j++) {

	    	// Generate sample data (Value & Muarelli, 1983)
	    	int seed = seeds(i * iterations + j - 1);
    	    mat sample = getSample(n, intermediate_P, seed, FTable);
    	    mat R = cor(sample);
	    	vec moments = compute4thOrderMoments(sample);
    	    p_serafini(j) = serafini2019(R, n, delta, moments);
    	    if (p_serafini(j) <= .05) {
    			counter_serafini++;
    	    }

    	    p_counsell(j) = counsell2015(R, 	n, delta);
    	    if (p_counsell(j) <= .05) {
    			counter_counsell++;
    	    }

			// Calculate/sum the skewness and kurtosis of each marginal distribution
			for (int k = 0; k < 3; k++) {
				kurtosis_observed = kurtosis_observed + kurtosis(sample.col(k));
				skewness_observed = skewness_observed + skewness(sample.col(k));
			}

	}

	// Calculate the rejection rate for each statistic
    	double alpha_serafini = (double)counter_serafini / iterations;
    	double alpha_counsell = (double)counter_counsell / iterations;

	// Calculate the Kolmogorov-Smirnov D for the p values for each statistic
    	double D_counsell = kolmogorovD(p_counsell);
    	double D_serafini =  kolmogorovD(p_serafini);

	// Calculate the empirical skewness and kurtossi for this condition
	skewness_observed = skewness_observed / (iterations * 3);
	kurtosis_observed = kurtosis_observed / (iterations * 3);

	// Print output (to be piped to a .csv)
        cout << skewness_nominal1  << ","
	     << kurtosis_nominal1  << ","
		 << skewness_nominal2  << ","
	     << kurtosis_nominal2  << ","
		 << skewness_nominal3  << ","
	     << kurtosis_nominal3  << ","
	     << rho12             << ","
	     << rho23             << ","
	     << "serafini"        << ","
	     << alpha_serafini    << ","
	     << D_serafini        << ","
             << skewness_observed << ","
	     << kurtosis_observed << endl; 

        cout << skewness_nominal1  << ","
	     << kurtosis_nominal1  << ","
		 << skewness_nominal2  << ","
	     << kurtosis_nominal2  << ","
		 << skewness_nominal3  << ","
	     << kurtosis_nominal3  << ","
	     << rho12             << ","
	     << rho23             << ","
	     << "counsell"        << ","
	     << alpha_counsell    << ","
	     << D_counsell        << ","
             << skewness_observed << ","
	     << kurtosis_observed << endl; 

    }

    return 0;
}

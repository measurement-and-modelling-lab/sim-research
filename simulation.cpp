#include <armadillo>
#include <thread>
#include <iostream>
#include <time.h>
#include <mutex>
#include <atomic>

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
#include <thread>




// The sum of the skewness and kurtosis for each marginal of each sample
std::atomic<double> skewness_observed;
std::atomic<double> kurtosis_observed;

// Counters for the number of rejections
std::atomic<int> counter_serafini;
std::atomic<int> counter_counsell;

// Extract condition parameters
int n;
double rho12;
double rho23;
mat FTable;
vec seeds;
mat intermediate_P;
int iterations;
double delta = 0.1;
std::mutex mu;

double atomic_addD(std::atomic<double> &f, double d)
{
    double old = f.load(std::memory_order_consume);
    double desired = old + d;
    while (!f.compare_exchange_weak(old, desired, std::memory_order_release, std::memory_order_consume))
    {
        desired = old + d;
    }
    return desired;
}

void slowPart(int start,
              int end,
              int i,
              vec &p_serafini_param,
              vec &p_counsell_param)
{
    for (int j = start; j < end; j++)
    {
        // Generate sample data (Value & Muarelli, 1983)
        mu.lock();
        int seed = seeds(i * iterations + j - 1);
        mat sample = getSample(n, intermediate_P, seed,FTable);
        mu.unlock();
        mat R = cor(sample);
        vec moments = compute4thOrderMoments(sample);
        // Calculate/sum the skewness and kurtosis of each marginal distribution
        for (int k = 0; k < 3; k++)
        {
            while (atomic_addD(kurtosis_observed, kurtosis(sample.col(k))) != kurtosis_observed)
                ;
            while (atomic_addD(skewness_observed, skewness(sample.col(k))) != skewness_observed)
                ;
        }
        p_serafini_param(j) = serafini2019(R, n, delta, moments);
        if (p_serafini_param(j) <= .05)
        {
            counter_serafini.fetch_add(1, std::memory_order_relaxed);
        }

        p_counsell_param(j) = counsell2015(R, n, delta);
        if (p_counsell_param(j) <= .05)
        {
            counter_counsell.fetch_add(1, std::memory_order_relaxed);
            ;
        }

    }
}

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

int main(int argc, char ** argv) {

    // The number of iterations per condition
    iterations = atoi(argv[1]);

    mat conditions;
    conditions.load(argv[2]);
	mat coefficients;
	coefficients.load("coefficients.csv");

	// Only shed row if the csv has a header row
	coefficients.shed_row(0);

    seeds.load(argv[3]);

    // Print header for .csv
    cout << "skewness_nominal1,kurtosis_nominal1,skewness_nominal2,kurtosis_nominal2,skewness_nominal3,kurtosis_nominal3,rho12,rho23,statistic,alpha,D,skewness_observed,kurtosis_observed" << endl;

    delta = 0.1;

    for (int i = 1; i < conditions.n_rows; i++) {

    	// The p values for each statistic
    	vec p_serafini(iterations);
    	vec p_counsell(iterations);

	// The sum of the skewness and kurtosis for each marginal of each sample
    	skewness_observed = 0;
        kurtosis_observed = 0;

    	// Counters for the number of rejections
    	counter_serafini = 0;
    	counter_counsell = 0;
	// Extract condition parameters
		int skewness_nominal1 = conditions(i, 0);
		int kurtosis_nominal1 = conditions(i, 1);
		int skewness_nominal2 = conditions(i, 2);
		int kurtosis_nominal2 = conditions(i, 3);
		int skewness_nominal3 = conditions(i, 4);
		int kurtosis_nominal3 = conditions(i, 5);
    	n = conditions(i, 6);
    	double rho12 = conditions(i, 7);
    	double rho23 = conditions(i, 8);
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

	intermediate_P = getIntermediateP(P, FTable);
    unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
        {
            std::thread range1(
                slowPart,
                0,
                iterations / 4,
                i,
                std::ref(p_serafini),
                std::ref(p_counsell)

            );
            std::thread range2(
                slowPart,
                iterations / 4,
                iterations / 2,
                i,
                std::ref(p_serafini),
                std::ref(p_counsell)

            );
            std::thread range3(
                slowPart,
                iterations / 2,
                iterations / 4 * 3,
                i,
                std::ref(p_serafini),
                std::ref(p_counsell));

            slowPart(iterations / 4 * 3,
                     iterations,
                     i,
                     std::ref(p_serafini),
                     std::ref(p_counsell));

            range1.join();
            range2.join();
            range3.join();
        }
        else
        {
            std::vector<std::thread> threads;
            threads.reserve(num_threads - 1);
            for (int k = 0; k < num_threads - 1; k++)
            {
                threads.emplace_back(std::thread(slowPart,
                                                 iterations / num_threads * k,
                                                 iterations / num_threads * (k + 1),
                                                 i,
                                                 std::ref(p_serafini),
                                                 std::ref(p_counsell)));
            }
            slowPart(
                iterations / (num_threads) * (num_threads - 1),
                iterations,
                i,
                std::ref(p_serafini),
                std::ref(p_counsell));
            for (int k = 0; k < threads.size(); k++)
            {
                threads[k].join();
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

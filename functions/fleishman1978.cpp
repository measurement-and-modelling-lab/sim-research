#define ARMA_NO_DEBUG // Disable bound checks to improve speed
#include <armadillo>
#include "optim.hpp"
#include "fleishman1978.h"
using namespace arma;

//Objective function to minimize
double fleishman_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
    const double b = vals_inp(0);
    const double c = vals_inp(1);
    const double d = vals_inp(2);
    // cast void pointer to arma::vec* to access parameter like a vector
    vec* data = static_cast<vec*>(opt_data);
    double skewness = data->at(0);
    double kurtosis = data->at(1);
    double eq1 = pow(b,2) + 6*b*d + 2*pow(c,2) + 15*pow(d,2) - 1;
    double eq2 = 2*c*(pow(b,2) + 24*b*d + 105*pow(d,2) + 2) - skewness;
    double eq3 = 24*(b*d + pow(c,2)*(1 + pow(b,2) + 28*b*d) + pow(d,2)*(12 + 48*b*d + 141*pow(c,2) + 225*pow(d,2))) - kurtosis;
    vec eq = {eq1,eq2,eq3};
    double obj_val = sum(square(eq));
    return obj_val;
}

rowvec fleishman1978(double skewness,double kurtosis){

    vec opt_params = {skewness,kurtosis};
    vec coeffs = {1,0,0};
    bool success = optim::nm(coeffs,fleishman_fn,&opt_params);
    if(!success){
    std::cout << "ERROR: Failed to compute fleishman coeffcients for skewness=" << skewness << " and kurtosis=" << kurtosis << std::endl;
    }
    rowvec result = {coeffs(0),coeffs(1),coeffs(2)};
    return result;
}


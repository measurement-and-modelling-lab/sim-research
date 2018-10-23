#include <iostream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <armadillo>
using namespace std;
using namespace arma;
// mvrnonnorm <- function(n, mu, Sigma, skewness = NULL, kurtosis = NULL, empirical = FALSE) {
//     p <- 3
//     Z <- ValeMaurelli1983(n = n, COR = cov2cor(Sigma), skewness = skewness, kurtosis = kurtosis)
//     TMP <- scale(Z, center = FALSE, scale = 1/sqrt(diag(Sigma)))[ , , drop = FALSE] ## Divide the ith column of Z by ith diagonal element of Sigma
//         ## int i;
//         ## for (i=0; i<nvar; i++) {
//         ##   uvec indices(1);
//         ##   indices(0) = i;
//         ##   double element = Sigma(i,i);
//         ##   Z = Z.each_col(indices) / linspace<vec>(element,element,nvar); ## Should probably be assigned to TMP
//         ## }
//     X <- sweep(TMP, MARGIN = 2, STATS = mu, FUN = "+") ## mat X = TMP.for_each( [](mat::elem_type& val) { val += mu; } );
//     X
// }
mat ValeMaurelli1983(n = 100L, COR, skewness, kurtosis, debug = FALSE) {
	rowvec fleishman1978(double skewness, double kurtosis) { 
		  //replace this with function for fetching from matrix
			static mat fleishmanTable;
	  		fleishmanTable.load("coefficients.csv");
		  	int index = -1;
		    for (int i = 0; i < fleishmanTable.n_rows; i++){
	  			if((fleishmanTable(i,0) == skewness) && (fleishmanTable(i,1) == kurtosis)){
	  			index = i;
	  			}
	  		}

	  		if(index == -1){
	  			return NULL;
	  		}
	      else{
	  			double a,b,c,d;
	  			b = fleishmanTable(index,2);
	  			c = fleishmanTable(index,3);
	  			d = fleishmanTable(index,4);
	  			a = -1*c;
	  			rowvec values(4);
	        	values(0)=a;
	        	values(1)=b;
	        	values(2)=c;
	        	values(3)=d;
	        return values;
	  		}

	  }

//     getICOV <- function(b1, c1, d1, b2, c2, d2, R) {

//         ## Find value of rho that multiplies with b1, c1, d1, b2, c2, d2 to produce a correlation as close as possible to R
//         objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {

//             rho <- x[1L]
//             eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + rho^2*(2*c1*c2) + rho^3*(6*d1*d2) - R
//             return(eq^2)

//         }

//         out <- nlminb(start=R, objective=objectiveFunction, scale=10, control=list(trace=0), b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)

//         rho <- out$par[1L]

//         return(rho)
//     }

//     nvar <- ncol(COR) ##  int nvar = COR.n_cols;


// //create Fleishman table
	mat FTable = zeroes<mat>(nvar,4);
	int nvar = 10;
	mat FTable = zeros<mat>(nvar,4);
	for (int i=0;i<nvar; i++) {
		rowvec result = fleishman1978(1.25,3);
		FTable.cols(0,3).each_row() = fleishman1978(SK[i],KU[i]);
  	}


// compute intermediate correlations between all pairs
   int i;
    int j;
    mat ICOR = eye<mat>(nvar,nvar);
    for (j=0; j<nvar-1; j++) {
        for (i=j+1; i<nvar-1; i++) {
            if (COR(i,j) == 0) {
              continue;
        	} else {
                   ICOR(i,j) = getICOV(FTable(i,2), FTable(i,3), FTable(i,4), FTable(j,2), FTable(j,3), FTable(j,4), COR(i,j));
                   ICOR(j,i) = ICOR(i,j);
            }
        }
  }

  //     X <- Z <- MASS::mvrnorm(n=n, mu=rep(0,nvar), Sigma=ICOR)
//         ## arma_rng::set_seed(1234567); ## Seed with our RNG
//         ## vec mu = zeros<vec>(nvar);
//         ## Z = mvnrnd(mu, ICOR, nvar);

//     ## transform Z using Fleishman constants
//     for (i in 1:nvar) {
//         X[,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]^2 + FTable[i,4L]*Z[,i]^3
//     }

//     ## for (i=0; i<nvar; i++) {
//     ##     vec Z1 = Z.col(i)
//     ##     vec Z2 = Z1.transform( [](int val) { return (val * val); } ); ## Zi^2
//     ##     vec Z3 = Z1.transform( [](int val) { return (val * val * val); } ); ## Zi^3
//     ##     uvec indices(1);
//     ##     indices(0) = i;
//     ##     X = X.each_col(indices) = FTable(i,1) + FTable(i,2) * Z1 + FTable(i,3) * Z2 + FTable(i,4) * Z3 ## 100% untested
//     ## }

    
//     return(X)
    
}
int main (int argc, char const *argv[]) {

  return 0;
}




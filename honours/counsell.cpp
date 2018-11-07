#include <armadillo>
#include <math.h> 
using namespace arma;
double counsell(double r12, double r13, double r23, double n, double delta){
	double detR = (1-pow(r12,2) - pow(r13,2)-pow(r23,2)) + (2 * r12 * r13 * r23);
	double s = sqrt(((n - 1) * (1 + r23))/((2 * ((n - 1)/(n - 3)) * detR) + ((pow((r12 + r13),2))/4) * (pow((1 - r23),3))));
	double p1 = normcdf((abs(r12 - r13) - delta) * s);
	double p2 = normcdf((-abs(r12 - r13) - delta) * s);
	return (p1-p2);
}
int main(int argc, char const *argv[])
{
	double r12 = .5;
	double r13 = .4;
	double r23 = .6;
	double n =10;
	double delta = 0.1;
	cout << counsell(r12, r13, r23, n, delta) << endl;
	return 0;
}
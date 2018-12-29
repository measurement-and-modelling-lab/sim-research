double zADF(mat R, int n, double delta, vec moments) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the difference between r12 and r13 under H0
    // moments is a vector of kurtosis estimates
    // Returns a p value for test of |r12-r13| <= delta

    double r12 = R(0, 1);
    double r13 = R(0, 2);
    double r23 = R(1, 2);

    // Z ADF

    double gamma_12 = adfCov(0, 1, 0, 1, R, moments);
    double gamma_13 = adfCov(0, 2, 0, 2, R, moments);
    double gamma_12_13 = adfCov(0, 1, 0, 2, R, moments);

    double z1 = sqrt(n) * (fabs(r12 - r13) - delta) *
    	(1 / sqrt(gamma_12 + gamma_13 - 2 * gamma_12_13));
    double z2 = sqrt(n) * (-fabs(r12 - r13) - delta) *
    	(1 / sqrt(gamma_12 + gamma_13 - 2 * gamma_12_13));
    double p = normcdf(z1) - normcdf(z2);

    return p;
}

double zNT(mat R, int n, double delta, vec moments) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the difference between r12 and r13 under H0
    // moments is a vector of kurtosis estimates
    // Returns a p value for test of |r12-r13| <= delta

    double r12 = R(0, 1);
    double r13 = R(0, 2);
    double r23 = R(1, 2);

    double gamma_12 = pow(1 - pow(r12, 2), 2);
    double gamma_13 = pow(1 - pow(r13, 2), 2);
    double gamma_12_13 = r23 * (1 - pow(r13, 2) - pow(r12, 2)) -
	0.5 * (r12 * r13) * (1 - pow(r13, 2) - pow(r12, 2) - pow(r23, 2));

    double z1 = sqrt(n) * (fabs(r12 - r13) - delta) *
    	(1 / sqrt(gamma_12 + gamma_13 - 2 * gamma_12_13));
    double z2 = sqrt(n) * (-fabs(r12 - r13) - delta) *
    	(1 / sqrt(gamma_12 + gamma_13 - 2 * gamma_12_13));
    double p = normcdf(z1) - normcdf(z2);

    return p;
}

double zStarADF(mat R, int n, double delta, vec moments) {
    // R is a correlation matrix
    // n is the sample size of the data set
    // delta is the difference between r12 and r13 under H0
    // moments is a vector of kurtosis estimates
    // Returns a p value for test of |r12-r13| <= delta

    double r12 = R(0, 1);
    double r13 = R(0, 2);
    double r23 = R(1, 2);

    double delta_mod = fabs((0.5 * log((1 + ((r12 + r13) / 2 - delta / 2)) /
				       (1 - ((r12 + r13) / 2 - delta / 2)))) -
			    (0.5 * log((1 + ((r12 + r13) / 2 + delta / 2)) /
				       (1 - ((r12 + r13) / 2 + delta / 2)))));


    double gamma_12 = adfCov(0, 1, 0, 1, R, moments);
    double gamma_13 = adfCov(0, 2, 0, 2, R, moments);
    double gamma_12_13 = adfCov(0, 1, 0, 2, R, moments);
    double psi_prime = gamma_12_13 * (1/sqrt(gamma_13)) * (1/sqrt(gamma_13));

    // double psi_prime = gamma_12_13 /
    // 	   ((1 - pow(r12, 2)) * (1 - pow(r13, 2)));

    double z1 = sqrt(n - 3) *
	        (fabs(fisher(r12) - fisher(r13)) - delta_mod) *
	        (1 / sqrt(2 - 2 * psi_prime));
    double z2 = sqrt(n - 3) *
	        (-fabs(fisher(r12) - fisher(r13)) - delta_mod) *
	        (1 / sqrt(2 - 2 * psi_prime));
    double p = normcdf(z1) - normcdf(z2);

    return p;
}

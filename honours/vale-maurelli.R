## cov2cor
## mvrnorm

mvrnonnorm <- function(n, mu, Sigma, skewness = NULL, kurtosis = NULL, empirical = FALSE) {
    p <- 3
    Z <- ValeMaurelli1983(n = n, COR = cov2cor(Sigma), skewness = skewness, kurtosis = kurtosis)
    TMP <- scale(Z, center = FALSE, scale = 1/sqrt(diag(Sigma)))[ , , drop = FALSE] ## Divide the ith column of Z by ith diagonal element of Sigma
        ## int i;
        ## for (i=0; i<nvar; i++) {
        ##   uvec indices(1);
        ##   indices(0) = i;
        ##   double element = Sigma(i,i);
        ##   Z = Z.each_col(indices) / linspace<vec>(element,element,nvar); ## Should probably be assigned to TMP
        ## }
    X <- sweep(TMP, MARGIN = 2, STATS = mu, FUN = "+") ## mat X = TMP.for_each( [](mat::elem_type& val) { val += mu; } );
    X
}


ValeMaurelli1983 <- function(n = 100L, COR, skewness, kurtosis, debug = FALSE) {

    fleishman1978 <- function(skewness, kurtosis) { ## replace this with function for fetching from matrix

        system.function <- function(x, skewness, kurtosis) {
            b.=x[1L]; c.=x[2L]; d.=x[3L]
            eq1 <- b.^2 + 6*b.*d. + 2*c.^2 + 15*d.^2 - 1
            eq2 <- 2*c.*(b.^2 + 24*b.*d. + 105*d.^2 + 2) - skewness
            eq3 <- 24*(b.*d. + c.^2*(1 + b.^2 + 28*b.*d.) + d.^2*(12 + 48*b.*d. + 141*c.^2 + 225*d.^2)) - kurtosis
            eq <- c(eq1,eq2,eq3)
            sum(eq^2) ## SS
        }

        out <- nlminb(start = c(1, 0, 0), objective = system.function, scale = 10, control = list(trace = 0), skewness = skewness, kurtosis = kurtosis)

        b. <- out$par[1L]
        c. <- out$par[2L]
        d. <- out$par[3L]
        a. <- -c.

        return(c(a.,b.,c.,d.))

    }

    getICOV <- function(b1, c1, d1, b2, c2, d2, R) {

        ## Find value of rho that multiplies with b1, c1, d1, b2, c2, d2 to produce a correlation as close as possible to R
        objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {

            rho <- x[1L]
            eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) + rho^2*(2*c1*c2) + rho^3*(6*d1*d2) - R
            return(eq^2)

        }

        out <- nlminb(start=R, objective=objectiveFunction, scale=10, control=list(trace=0), b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)

        rho <- out$par[1L]

        return(rho)
    }

    nvar <- ncol(COR) ##  int nvar = COR.n_cols;

    ## create Fleishman table
    FTable <- matrix(0, nvar, 4L) ## mat FTable = zeroes<mat>(nvar,4);
    for (i in 1:nvar) { ## Fetch values from coefficients.csv
        FTable[i,] <- fleishman1978(skewness=SK[i], kurtosis=KU[i])
    }

    ## compute intermediate correlations between all pairs
    ICOR <- diag(nvar) ## mat ICOR = eye<mat>(nvar,nvar);
    for (j in 1:(nvar-1L)) {
        for (i in (j+1):nvar) {
            if (COR[i,j] == 0) {
                next
            } else {
                ICOR[i,j] <- getICOV(FTable[i,2], FTable[i,3], FTable[i,4], FTable[j,2], FTable[j,3], FTable[j,4], R=COR[i,j])
                ICOR[j,i] <- ICOR[i,j]
            }
        }
    }

  ##   int i;
  ##   int j;
  ##   mat ICOR = eye<mat>(nvar,nvar);
  ##   for (j=0; j<nvar-1; j++) {
  ##       for (i=j+1; i<nvar-1; i++) {
  ##           if (COR(i,j) == 0) {
  ##             continue;
  ##       	} else {
  ##                  ICOR(i,j) = getICOV(FTable(i,2), FTable(i,3), FTable(i,4), FTable(j,2), FTable(j,3), FTable(j,4), COR(i,j));
  ##                  ICOR(j,i) = ICOR(i,j);
  ##           }
  ##       }
  ## }
    
    

    X <- Z <- MASS::mvrnorm(n=n, mu=rep(0,nvar), Sigma=ICOR)
        ## arma_rng::set_seed(1234567); ## Seed with our RNG
        ## vec mu = zeros<vec>(nvar);
        ## Z = mvnrnd(mu, ICOR, nvar);

    ## transform Z using Fleishman constants
    for (i in 1:nvar) {
        X[,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]^2 + FTable[i,4L]*Z[,i]^3
    }

    ## for (i=0; i<nvar; i++) {
    ##     vec Z1 = Z.col(i)
    ##     vec Z2 = Z1.transform( [](int val) { return (val * val); } ); ## Zi^2
    ##     vec Z3 = Z1.transform( [](int val) { return (val * val * val); } ); ## Zi^3
    ##     uvec indices(1);
    ##     indices(0) = i;
    ##     X = X.each_col(indices) = FTable(i,1) + FTable(i,2) * Z1 + FTable(i,3) * Z2 + FTable(i,4) * Z3 ## 100% untested
    ## }

    
    return(X)
}

fleishman1978 <- function(skewness, kurtosis) {
    system.function <- function(x, skewness, kurtosis) {
        b.=x[1L]; c.=x[2L]; d.=x[3L]
        eq1 <- b.^2 + 6*b.*d. + 2*c.^2 + 15*d.^2 - 1
        eq2 <- 2*c.*(b.^2 + 24*b.*d. + 105*d.^2 + 2) - skewness
        eq3 <- 24*(b.*d. + c.^2*(1 + b.^2 + 28*b.*d.) +
                   d.^2*(12 + 48*b.*d. + 141*c.^2 + 225*d.^2)) - kurtosis
        eq <- c(eq1,eq2,eq3)
        sum(eq^2) ## SS
    }

    out <- nlminb(start = c(1, 0, 0), objective = system.function,
                  scale = 10, control = list(trace = 0),
                  skewness = skewness, kurtosis = kurtosis)
    if(out$convergence != 0) warning("no convergence")
    b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
    c(a.,b.,c.,d.)
}

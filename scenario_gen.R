library(stringr) ## str_split_fixed

skew <- c(0, 0, 0, 1, 0, 1)
kurtosis <- c(-1, 0, 1, 1, 3, 3)
n <- c(125, 250, 500)
rho12 <- c(0.2, 0.4, 0.6, 0.8)
rho23 <- c(0.2, 0.5, 0.8)

table <- expand.grid(paste(skew,kurtosis,sep=","), n, rho12, rho23)
separated <- str_split_fixed(table$Var1, ",", 2)
mode(separated) <- "numeric"
table <- cbind(separated, table[,-1], 0, 0, 0)
colnames(table) <- c("skew", "kurtosis", "n", "rho12", "rho23", "b", "c", "d")

coeff.table <- read.csv("~/git/sim-research/blah/coefficients.csv")
coeff.table <- as.matrix(coeff.table)
for (i in 1:nrow(table)) {
    table[i, 6:8] <- coeff.table[coeff.table[,1] == table[i,"skew"] &
                                 coeff.table[,2] == table[i,"kurtosis"],
                                 3:5, drop=TRUE]
}

write.table(x=table, file="./conditions.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE, append=FALSE)

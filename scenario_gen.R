library(stringr) ## str_split_fixed
source("fleishman1978.R")

skewness <- c(0, 0, 0, 1, 0, 1, 0, 1, 2)
kurtosis <- c(-1, 0, 2, 2, 6, 6, 18, 18, 18)
n <- c(1000)
rho12 <- c(0.2, 0.4, 0.6, 0.8)
rho23 <- c(0.3, 0.5, 0.7)

table <- expand.grid(paste(skewness,kurtosis,sep=","), n, rho12, rho23)
separated <- str_split_fixed(table$Var1, ",", 2)
mode(separated) <- "numeric"
table <- cbind(separated, table[,-1], 0, 0, 0)
colnames(table) <- c("skewness", "kurtosis", "n", "rho12", "rho23", "b", "c", "d")

for (i in 1:nrow(table)) {
    coeff <- fleishman1978(table[i,1], table[i,2])
    table[i, 6:8] <- coeff[-1]
}

write.table(x=table, file="./conditions.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, append=FALSE)

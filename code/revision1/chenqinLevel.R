source("./stat.R")
n <- 300
n1 <- n
n2 <- n
p <- 200
myResult <- vector()
for (i in 1:20000){
    X1 <- matrix(rnorm(n1*p),nrow = n1)
    X2 <- matrix(rnorm(n2*p),nrow = n2)
    myResult[i] <- chenStat(X1,X2,n1,n2)$stat * n1 * n2 /( sqrt(2*p)*(n1+n2) )
}
qqnorm(myResult)
abline(a= 0,b=1)
mean(pnorm(myResult,lower.tail = FALSE)<0.05)

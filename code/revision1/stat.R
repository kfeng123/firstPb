library(ggplot2)
library(MASS)


newModelGenerator <- function(eigValue) {
    # generate V
    V <- svd(matrix(rnorm(p * p, 0, 1), c(p, p)))$v
    normalModelSimulator <- function(n) {
        t(V %*% diag(eigValue) %*% matrix(rnorm(p * n), c(p, n)))
    }
    return(list(V = V,
                normalModelSimulator = normalModelSimulator))
}

# Chen and Qin (2010)
chenStat <- function(X1, X2, n1, n2,...) {
    T1 <- sum(colMeans(X1) ^ 2) * n1 / (n1 - 1) - sum(X1 ^ 2) / n1 / (n1 - 1)
    T2 <- sum(colMeans(X2) ^ 2) * n2 / (n2 - 1) - sum(X2 ^ 2) / n2 / (n2 - 1)
    T3 <- sum(colMeans(X1) * colMeans(X2))
    list(stat=T1 + T2 - 2 * T3)
}
# Srivastava and Du (2008)
sriStat <- function(X1, X2, n1, n2,...) {
    S1 <- colVars(X1)
    S2 <- colVars(X2)
    S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
    list(stat=sum((colMeans(X1) - colMeans(X2)) ^ 2 / S))
}

# Ma et. al. (2015)
maTest <- function(X1,X2,n1,n2,...){
    
    S1 <- var(X1)
    S2 <- var(X2)
    S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
    myEigen <- eigen(S, symmetric = TRUE)
    
    p <- ncol(X1)
    TFast <- sum((colMeans(X1)-colMeans(X2))^2)*(n1*n2/(n1+n2))/p-(myEigen$values[-(1:r)])/p
    
    refDis <- vector()
    for (i in 1:500){
        refDis[i] <- 0
        for (j in 1:r)
            refDis[i] <- refDis[i] + myEigen$values[j]*(rnorm(1)^2)/p
    }
    TFast > quantile(refDis,0.95)
    
}

# New chi squared test
myChiTest <- function(X1,X2,n1,n2,...){
    # Chen's statistic over tau
    tau <- 1 / n1 + 1 / n2
    temp1 <- chenStat(X1,X2,n1,n2,...)$stat/tau
    
    # variance estimator
    S1 <- var(X1)
    S2 <- var(X2)
    S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
    myEigen <- eigen(S, symmetric = TRUE)
    sigmaSqEst <- mean(myEigen$values[-(1:r)])
    
    p <- ncol(X1)
    refDis <- vector()
    for (i in 1:500){
        refDis[i] <- sqrt(2*p)*sigmaSqEst*rnorm(1)
        for (j in 1:r)
            refDis[i] <- refDis[i] + myEigen$values[j]*(rnorm(1)^2-1)
    }
    temp1 > quantile(refDis,0.95)
}


# New test
myStat <- function(X1, X2, n1, n2, ...) {
    r <- list(...)$r
    S1 <- var(X1)
    S2 <- var(X2)
    S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
    myTildeV <- myEigen$vectors[, -(1:r)]
    
    trace1 <- sum(eigen(S1, symmetric = TRUE)$values[-(1:r)])
    trace2 <- sum(eigen(S2, symmetric = TRUE)$values[-(1:r)])
    
    # variance estimator
    myEigen <- eigen(S, symmetric = TRUE)
    sigmaSqEst <- mean(myEigen$values[-(1:r)])
    
    stat <-
        sum((t(myTildeV) %*% (colMeans(X1) - colMeans(X2))) ^ 2) - trace1 / n1 -
        trace2 / n2
    
    # studentized statistic
    tau <- 1 / n1 + 1 / n2
    studentStat <- stat / sigmaSqEst / sqrt(2 * tau ^ 2 * p)
    
    return(list(stat = stat,
                studentStat = studentStat))
}



# do test by studentized statistic and normal theory
doTest = function(X1, X2, n1, n2, p, rmax = 10) {
    my = myStat(X1, X2, n1, n2, rmax)
    return(pnorm(my$studentStat, 0, 1, lower.tail = FALSE))
}

# unequal variance test
doUneqTest <- function(X1, X2, n1, n2, p, rmax = 10) {
    S1 <- var(X1)
    S2 <- var(X2)
    myEigen <- eigen(S1, symmetric = TRUE)
    theTemp = myEigen$values[1:(n1 - 2)] / myEigen$values[2:(n1 - 1)]
    myRhat1 = which.max(theTemp[1:rmax])
    myV1 = myEigen$vectors[, (1:myRhat1)]
    
    myEigen = eigen(S2, symmetric = TRUE)
    theTemp = myEigen$values[1:(n2 - 2)] / myEigen$values[2:(n2 - 1)]
    myRhat2 = which.max(theTemp[1:rmax])
    myV2 <- myEigen$vectors[, (1:myRhat2)]
    
    myV <- svd(cbind(myV1, myV2))$u
    myTildeV <- Null(myV)
    
    
    trace1 = sum(eigen(S1, symmetric = TRUE)$values[(myRhat1 + 1):p])
    trace2 = sum(eigen(S2, symmetric = TRUE)$values[(myRhat2 + 1):p])
    
    # variance estimator
    sigmaSqEst <- 2 * (p - myRhat1 - myRhat2) * (
        (trace1 / (p - myRhat1)) ^ 2 / n1 / (n1 - 1) +
            (trace2 / (p - myRhat2)) ^ 2 / n2 / (n2 - 1) +
            2 * (trace1 / (p - myRhat1)) * (trace2 / (p - myRhat2)) / n1 /
            n2
    )
    
    
    stat <-
        sum((t(myTildeV) %*% (colMeans(X1) - colMeans(X2))) ^ 2) - trace1 / n1 -
        trace2 / n2
    
    studentStat <- stat / sqrt(sigmaSqEst)
    
    return(pnorm(studentStat, 0, 1, lower.tail = FALSE))
}


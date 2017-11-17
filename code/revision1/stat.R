library(ggplot2)
library(MASS)
library(highD2pop)
library(matrixStats)
library(Hotelling)

newModelGenerator <- function(eigValue,theDistribution="normal") {
    p <- length(eigValue)
    # generate V
    V <- svd(matrix(rnorm(p * p, 0, 1), c(p, p)))$v
    if(theDistribution=="normal"){
        modelSimulator <- function(n) {
            t(V %*% diag(eigValue^(1/2)) %*% matrix(rnorm(p * n), c(p, n)))
        }
    }
    if(theDistribution=="chiSquared"){
        modelSimulator <- function(n) {
            t(V %*% diag(eigValue^(1/2)) %*% matrix((rchisq(p * n,df=4)-4)/sqrt(8), c(p, n)))
        }
    }
    if(theDistribution=="t"){
        modelSimulator <- function(n) {
            t(V %*% diag(eigValue^(1/2)) %*% matrix(rt(p * n,df = 4)/sqrt(2), c(p, n)))
        }
    }
    
    return(list(V = V,
                eigValue = eigValue,
                modelSimulator = modelSimulator))
}



# abandoned Chen and Qin (2010)
chenStat <- function(X1, X2, n1, n2,...) {
    T1 <- sum(colMeans(X1) ^ 2) * n1 / (n1 - 1) - sum(X1 ^ 2) / n1 / (n1 - 1)
    T2 <- sum(colMeans(X2) ^ 2) * n2 / (n2 - 1) - sum(X2 ^ 2) / n2 / (n2 - 1)
    T3 <- sum(colMeans(X1) * colMeans(X2))
    list(stat=T1 + T2 - 2 * T3)
}
# Srivastava and Du (2008)
sdStat <- function(X1, X2, n1, n2,...) {
    p <- ncol(X1)
    
    if(is.null(list(...)$S)){
        S1 <- colVars(X1)
        S2 <- colVars(X2)
        S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
    }else{
        S <- list(...)$S
    }
    
    dS <- diag(S)
    R <- cov2cor(S)
    
    mean1 <- colMeans(X1)
    mean2 <- colMeans(X2)
    
    temp <- sum(R^2)
    c <- 1+ temp/p^(3/2)
    
    theNumerator <- n1 * n2 / (n1+n2) * sum( (mean1 - mean2)^2 / dS ) - (n1+n2-2)*p/(n1+n2-4)
    theDenominator <- sqrt( 2*( temp - p^2/(n1+n2-2) ) * c )
    
    list(stat=theNumerator/theDenominator)
}

# Lopes et. al. (2011)
ljwStat <- function(X1,X2,n1,n2,...){
    p <- ncol(X1)
    k <- floor((n1+n2-2)/2)
    Pk <- rnorm(p*k)
    dim(Pk) <- c(p,k)
    (hotelling.test(X1 %*% Pk, X2 %*% Pk)$pval<0.05) + 0
}

# Ma et. al. (2015)
maTest <- function(X1,X2,n1,n2,...){
    r <- list(...)$r
    
   if(is.null(list(...)$myEigen)){
        S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
   }else{
       myEigen <- list(...)$myEigen
   }
    
    p <- ncol(X1)
    TFast <- sum((colMeans(X1)-colMeans(X2))^2)*(n1*n2/(n1+n2))/p-sum(myEigen$values[-(1:r)])/p
    
    #refDis <- vector()
    #for (i in 1:500){
    #    refDis[i] <- 0
    #    for (j in 1:r)
    #        refDis[i] <- refDis[i] + myEigen$values[j]*(rnorm(1)^2)/p
    #}
    refDis <- matrix(rnorm(500*r),nrow = 500)^2 %*% myEigen$values[1:r] / p
    as.numeric(TFast > quantile(refDis,0.95))
    
}

# New chi squared test
myChiTest <- function(X1,X2,n1,n2,...){
    r <- list(...)$r
    p <- ncol(X1)
    
    # Chen's statistic over tau
    tau <- 1 / n1 + 1 / n2
    temp1 <- chenStat(X1,X2,n1,n2,...)$stat/tau
    
    # variance estimator
   if(is.null(list(...)$myEigen)){
        S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
   }else{
       myEigen <- list(...)$myEigen
   }
    oldSigmaSqEst <- mean(myEigen$values[-(1:r)])
    sigmaSqEst <- oldSigmaSqEst*(1+1/(n1+n2-2)*(r+oldSigmaSqEst*sum(1/(myEigen$values[1:r]-oldSigmaSqEst))))
    
    #refDis <- vector()
    #for (i in 1:500){
        #refDis[i] <- sqrt(2*p)*sigmaSqEst*rnorm(1)
        #for (j in 1:r)
            #refDis[i] <- refDis[i] + myEigen$values[j]*(rnorm(1)^2-1)
    #}
    refDis <- (matrix(rnorm(500*r),nrow = 500)^2-1) %*% myEigen$values[1:r] +
                    sqrt(2*p)*sigmaSqEst*rnorm(500)
    
    as.numeric(temp1 > quantile(refDis,0.95))
}


# New test
myStat <- function(X1, X2, n1, n2, ...) {
    r <- list(...)$r
    p <- ncol(X1)
    
    S1 <- var(X1)
    S2 <- var(X2)
   if(is.null(list(...)$myEigen)){
        S <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
   }else{
       myEigen <- list(...)$myEigen
   }
    myTildeV <- myEigen$vectors[, -(1:r)]
    
    trace1 <- sum(eigen(S1, symmetric = TRUE)$values[-(1:r)])
    trace2 <- sum(eigen(S2, symmetric = TRUE)$values[-(1:r)])
    
    # variance estimator
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

# New test version 2
myStat2 <- function(X1, X2, n1, n2, ...) {
    r <- list(...)$r
    p <- ncol(X1)
    
   if(is.null(list(...)$myEigen)){
        S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
   }else{
       myEigen <- list(...)$myEigen
   }
    myTildeV <- myEigen$vectors[, -(1:r)]
    
    # variance estimator
    oldSigmaSqEst <- mean(myEigen$values[-(1:r)])
    newSigmaSqEst <- oldSigmaSqEst*(1+1/(n1+n2-2)*(r+oldSigmaSqEst*sum(1/(myEigen$values[1:r]-oldSigmaSqEst))))
    
    tau <- 1 / n1 + 1 / n2
    stat <-
        sum((t(myTildeV) %*% (colMeans(X1) - colMeans(X2))) ^ 2) - tau*(p-r)*newSigmaSqEst
    
    # studentized statistic
    studentStat <- stat / newSigmaSqEst / sqrt(2 * tau ^ 2 * p)
    
    return(list(stat = stat,
                studentStat = studentStat))
}

# New test final version 
myStatFinal <- function(X1, X2, n1, n2, ...) {
    r <- list(...)$r
    p <- ncol(X1)
    
   if(is.null(list(...)$myEigen)){
        S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
   }else{
       myEigen <- list(...)$myEigen
   }
    myTildeV <- myEigen$vectors[, -(1:r)]
    
    # variance estimator
    oldSigmaSqEst <- mean(myEigen$values[-(1:r)])
    newSigmaSqEst <- (1-r/(n1+n2-2))^(-1) * oldSigmaSqEst
    # eigenvalue estimator
    myEigenEstimator <- myEigen$values[1:r]-((p+n1+n2-r-2)/(n1+n2-2))*newSigmaSqEst
    
    tau <- 1 / n1 + 1 / n2
    stat <-
        sum((t(myTildeV) %*% (colMeans(X1) - colMeans(X2))) ^ 2) -
        tau*(p-r)*newSigmaSqEst -
        tau * sum(p*newSigmaSqEst*myEigenEstimator/((n1+n2)*myEigenEstimator+(n1+n2+p)*newSigmaSqEst))
    
    # studentized statistic
    studentStat <- stat / newSigmaSqEst / sqrt(2 * tau ^ 2 * p)
    
    return(list(stat = stat,
                studentStat = studentStat))
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


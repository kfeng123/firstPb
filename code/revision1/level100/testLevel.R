source("../stat.R")
library(ggplot2)
library(xtable)



simulateLevel <- function(n1,n2,p,r,beta,theDistribution="normal",B=2000){
    theEig <- rep(1,p)
    if(r!=0){
        theEig[1:r] <- rep(p^beta,r)+runif(r,0,1)
    }
    temp <- newModelGenerator(theEig,theDistribution=theDistribution)
    modelSimulator <- temp$modelSimulator
    V <- temp$V[,1:r]
    
    
    myNew1Dis <- NULL
    myNew2Dis <- NULL
    oracleDis <- NULL
    chiDis <- NULL
    fastDis <- NULL
    CQDis <- NULL
    pb <- txtProgressBar(style=3)
    for(i in 1:B){
        # data generation
        X1 <- modelSimulator(n1)
        X2 <- modelSimulator(n2)
        
        # precalculate the eigenvalue decomposition
        S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
        
        # NEW1
        myNew1Dis[i] <- myStat(X1, X2, n1, n2, r=r, myEigen=myEigen)$studentStat
        
        # NEW2
        myNew2Dis[i] <- myStat2(X1, X2, n1, n2, r=r, myEigen=myEigen)$studentStat
        
        # oracle: V and sigma known
        temp <- chenStat(X1%*%Null(V),X2%*%Null(V),n1,n2)$stat
        oracleDis[i] <- n1*n2*temp/(sqrt(2*p)*(n1+n2)*1)
        
        # chi
        chiDis[i] <- myChiTest(X1, X2, n1, n2, r=r, myEigen=myEigen)
        
        # fast
        fastDis[i] <- maTest(X1, X2, n1, n2, r=r, myEigen=myEigen)
        
        # CQ
        CQDis[i] <- ChenQin.test(X1,X2)$ChQ
        
        setTxtProgressBar(pb,i/B)
    }
    close(pb)
    # NEW1
    myNew1Level <- sum(pnorm(myNew1Dis,lower.tail=FALSE)<0.05)/B
    #myNew1Level <- sum(pchisq(myNew1Dis*sqrt(2*(p-r))+p,lower.tail=FALSE)<0.05)/B
    # NEW2
    #myNew2Level <- sum(pnorm(myNew2Dis,lower.tail=FALSE)<0.05)/B
    myNew2Level <- sum(pchisq(myNew2Dis*sqrt(2*(p-r))+p-r,df=p-r,lower.tail=FALSE)<0.05)/B
    # ORACLE
    #oracleLevel <- sum(pnorm(oracleDis,lower.tail=FALSE)<0.05)/B
    oracleLevel <- sum(pchisq(oracleDis*sqrt(2*(p-r))+p-r,df=p-r,lower.tail=FALSE)<0.05)/B
    # chi
    chiLevel <- mean(chiDis)
    # fast
    fastLevel <- mean(fastDis)
    # CQ
    CQLevel <- sum(pnorm(CQDis,lower.tail=FALSE)<0.05)/B
    list(New1=myNew1Level,
         New2=myNew2Level,
         oracle=oracleLevel,
         chi=chiLevel,
         fast=fastLevel,
         CQ=CQLevel
         )
} 

# normal
Out=NULL

for(beta in c(0.5))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta)
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"level.csv")
Out=NULL

for(beta in c(1))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta)
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"level1.csv")

Out=NULL

for(beta in c(2))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta)
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"level2.csv")

# Chi-squared
Out=NULL

for(beta in c(0.5))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta, theDistribution = "chiSquared")
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"levelChi.csv")

Out=NULL

for(beta in c(1))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta, theDistribution = "chiSquared")
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"levelChi1.csv")

Out=NULL

for(beta in c(2))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta, theDistribution = "chiSquared")
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"levelChi2.csv")

# Student's t
Out=NULL

for(beta in c(0.5))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta, theDistribution = "t")
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"levelT.csv")

Out=NULL

for(beta in c(1))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta, theDistribution = "t")
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"levelT1.csv")

Out=NULL

for(beta in c(2))
for(r in 1)
for(p in c(200,400,600,800))
for(n in c(120)){
    n1=n
    n2=n
    level <- simulateLevel(n1,n2,p,r,beta, theDistribution = "t")
    level <- cbind(level)
    dimnames(level)[[2]] <- p
    Out <- cbind(Out,level)
}
write.csv(Out,"levelT2.csv")
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
    
    
    chiDis <- NULL
    fastDis <- NULL
    pb <- txtProgressBar(style=3)
    for(i in 1:B){
        # data generation
        X1 <- modelSimulator(n1)
        X2 <- modelSimulator(n2)
        
        # precalculate the eigenvalue decomposition
        S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
        myEigen <- eigen(S, symmetric = TRUE)
        
        # chi
        chiDis[i] <- myChiTest(X1, X2, n1, n2, r=r, myEigen=myEigen)
        
        # fast
        fastDis[i] <- maTest(X1, X2, n1, n2, r=r, myEigen=myEigen)
        
        
        setTxtProgressBar(pb,i/B)
    }
    close(pb)
    # chi
    chiLevel <- mean(chiDis)
    # fast
    fastLevel <- mean(fastDis)
    
    
    
    
    list(
         chi=chiLevel,
         fast=fastLevel
         )
} 

myOuterFun <- function(n,vecP){
    for(theDistribution in c("normal","chiSquared","t")){
        for(beta in c(0.5,1,2)){
            Out=NULL
            for(p in vecP){
                r<-2
                n1=n
                n2=n
                level <- simulateLevel(n1,n2,p,r,beta,theDistribution=theDistribution)
                level <- cbind(level)
                dimnames(level)[[2]] <- p
                Out <- cbind(Out,level)
            }
            write.csv(Out,paste0(n,theDistribution,beta,".csv"))
        }
    }
    
}
myOuterFun(50,c(200,500,800))
myOuterFun(100,c(200,500,800))
myOuterFun(150,c(200,500,800))











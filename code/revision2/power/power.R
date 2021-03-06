source("../stat.R")
# main simulation
myMainSimulation=function(p,n1,n2,r=2,beta,theDistribution,B=2000){
    
    theEig <- rep(1,p)
    if(r!=0){
        theEig[1:r] <- rep(p^beta,r)+runif(r,0,1)
    }
    temp <- newModelGenerator(theEig,theDistribution=theDistribution)
    modelSimulator <- temp$modelSimulator
    V <- temp$V[,1:r]
    
    tau <- 1/n1+1/n2
    theoryPower <- function(mu1,mu2){
        sum((mu1-mu2)^2)/sqrt(2*tau^2*p)
    }
    
    myNew3Power <-NULL
    chiPower <- NULL
    fastPower <- NULL
    CQPower <- NULL
    SDPower <- NULL
    ljwPower <- NULL
    
    temp1=rep(0,p)
    temp2=rnorm(p,0,1)
    for(hh in c(0,1,2,3,4)){
        myC=sqrt(hh*sqrt(2*tau^2*p)/sum((temp1-temp2)^2))
        mu1=temp1*myC
        mu2=temp2*myC
        
        myNew3Dis <- NULL
        chiDis <- NULL
        fastDis <- NULL
        CQDis <- NULL
        SDDis <- NULL
        ljwDis <- NULL
        pb <- txtProgressBar(style=3)
        for(i in 1:B){
            # data generation
            X1 <- modelSimulator(n1) + outer(rep(1,n1),mu1)
            X2 <- modelSimulator(n2) + outer(rep(1,n2),mu2)
            
            # precalculate the eigenvalue decomposition
            S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
            myEigen <- eigen(S, symmetric = TRUE)
            
            
            # NEW3
            myNew3Dis[i] <- myStatFinal(X1, X2, n1, n2, r=r, myEigen=myEigen)$studentStat
            # chi
            chiDis[i] <- myChiTest(X1, X2, n1, n2, r=r, myEigen=myEigen)
            # fast
            fastDis[i] <- maTest(X1, X2, n1, n2, r=r, myEigen=myEigen)
            # CQ
            CQDis[i] <- ChenQin.test(X1,X2)$ChQ
            # SD
            SDDis[i] <- sdStat(X1,X2, n1, n2, S=S)$stat
            # LJW
            ljwDis[i] <-ljwStat(X1,X2,n1,n2)
            
            setTxtProgressBar(pb,i/B)
        }
        close(pb)
        # NEW3
        myNew3Power <- c(myNew3Power,
              sum(pchisq(myNew3Dis*sqrt(2*(p-r))+p-r,df=p-r,lower.tail=FALSE)<0.05)/B
            )
        # chi
        chiPower <- c(chiPower,
                      mean(chiDis)
                      )
        # fast
        fastPower <- c(fastPower,
                       mean(fastDis)
        )
        # CQ
        CQPower <- c(CQPower,
                     sum(pnorm(CQDis,lower.tail=FALSE)<0.05)/B
        )
        # SD
        SDPower <- c(SDPower,
                     sum(pnorm(SDDis,lower.tail=FALSE)<0.05)/B
        )
        # LJW
        ljwPower <- c(ljwPower,
                      mean(ljwDis)
        )
    }
    
    data.frame(
         SNR=c(0,1,2,3,4),
         NEW=myNew3Power,
         CCQ=chiPower,
         FAST=fastPower,
         CQ=CQPower,
         SD=SDPower,
         LJW=ljwPower
         )
}

for(theDistribution in c("normal","chiSquared","t")){
    for(beta in c(0.5,1,2)){
        p <- 500
        r<-2
        n1 <- 100
        n2 <- 100
        thePower <- myMainSimulation(p,n1,n2,r,beta,theDistribution)
        write.csv(thePower,paste0(1,theDistribution,beta,".csv"),row.names=FALSE)
    }
}
    
for(n in c(50,100,150)){
    for(p in c(200,500,800)){
        n1 <- n
        n2 <- n
        r<-2
        beta <- 1
        theDistribution <- "normal"
        thePower <- myMainSimulation(p,n1,n2,r,beta,theDistribution)
        write.csv(thePower,paste0(2,"n",n,"p",p,".csv"),row.names=FALSE)
    }
}








library(ggplot2)
myPlot=function(uio){
    temp1 <- read.csv(uio)
    myD <- melt(temp1,id.vars="SNR")
    names(myD)<- c("SNR","Method","Power")
    thePlot <- ggplot(data=myD,aes(x=SNR,y=Power,color=Method,linetype=Method))+
        geom_line()+ylim(0,1)+
        theme_bw()+
        theme(
            axis.title.y= element_text(size=rel(2),angle=90),
            axis.title.x= element_text(size=rel(2)),
            axis.text=element_text(size=rel(1.8)),
            legend.position = c(0.05,0.95),
            legend.justification = c("left","top"),
            #legend.box.just ="right",
            legend.margin= margin(3,3,3,3),
            legend.box.background= element_rect(),
            #legend.key= element_rect(colour="black"),
            legend.text= element_text(size=rel(2)),
            legend.title= element_text(size=rel(2.2),face="bold"),
            legend.key.size = unit(1,"cm")
        )
        #scale_colour_brewer(palette = "Set1")
        #xlab(expression( paste("||",mu[1]-mu[2],"||"^2)))
}

for(theDistribution in c("normal","chiSquared","t")){
    for(beta in c(0.5,1,2)){
        p <- 500
        r<-2
        n1 <- 100
        n2 <- 100
        thePlot <- myPlot(paste0(1,theDistribution,beta,".csv"))
        ggsave(paste0("Power",1,theDistribution,beta,".eps"),thePlot)
    }
}
    
for(n in c(50,100,150)){
    for(p in c(200,500,800)){
        n1 <- n
        n2 <- n
        r<-2
        beta <- 1
        theDistribution <- "normal"
        thePlot <- myPlot(paste0(2,"n",n,"p",p,".csv"))
        ggsave(paste0("Power",2,"n",n,"p",p,".eps"),thePlot)
    }
}


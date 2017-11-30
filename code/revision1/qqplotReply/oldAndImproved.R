source("../stat.R")
library(ggplot2)
nnn <-5000

p <- 500
r <- 2
theEig <- rep(1,p)
beta <- 1
if(r!=0){
    theEig[1:r] <- rep(p^beta,r)
}
temp <- newModelGenerator(theEig,theDistribution="normal")
modelSimulator <- temp$modelSimulator


myQQplot <- function(x,y){
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    if (leny < lenx) 
        sx <- approx(1L:lenx, sx, n = leny)$y
    if (leny > lenx) 
        sy <- approx(1L:leny, sy, n = lenx)$y
    d <- data.frame(sx,sy)
    ggplot(d,aes(sx,sy))+
        geom_point(aes(colour=abs(sy-sx)),size=2)+
        scale_colour_gradient(low='royalblue',high='royalblue4',guide='none')+
        annotate("segment",x=-5,xend=5,y=-5,yend=5,colour="red")+
        xlab("Theoretical Quantiles")+
        ylab("Sample Quantiles")+
        scale_x_continuous(breaks=seq(-5,5,2.5))+
        scale_y_continuous(breaks=seq(-5,5,2.5))+
        theme_bw()+
        theme(
            axis.title.y= element_text(size=rel(2),angle=90),
            axis.title.x= element_text(size=rel(2)),
            axis.text= element_text(size=rel(1.8))
        )
}


n <- 50
n1 <- n
n2 <- n
pb <- txtProgressBar(style = 3)
myOld <- vector()
myImproved <- vector()
for (i in 1:nnn){
    # data generation
    X1 <- modelSimulator(n1)
    X2 <- modelSimulator(n2)
    
    myOld[i] <- myStat(X1,X2,n1,n2,r=r)$studentStat
    myImproved[i] <- myStat2(X1,X2,n1,n2,r=r)$studentStat
    setTxtProgressBar(pb, i/nnn)
}
close(pb)
# Normal Q-Q plot
#x <- rnorm(n=100000)
#thePlot <- myQQplot(x,temp1)
#ggsave("QQplotOld1.eps",thePlot)
# chi-Squared Q-Q plot
x <- rchisq(n=100000,df=p-r)
x <- (x-(p-r))/sqrt(2*(p-r))
thePlot <- myQQplot(x,myOld)
thePlot2 <- myQQplot(x,myImproved)
ggsave("QQplotOld50.eps",thePlot)
ggsave("QQplotImproved50.eps",thePlot2)

n <- 100
n1 <- n
n2 <- n
pb <- txtProgressBar(style = 3)
myOld <- vector()
myImproved <- vector()
for (i in 1:nnn){
    X1 <- matrix(rnorm(n1*p),nrow = n1)
    X2 <- matrix(rnorm(n2*p),nrow = n2)
    myOld[i] <- myStat(X1,X2,n1,n2,r=r)$studentStat
    myImproved[i] <- myStat2(X1,X2,n1,n2,r=r)$studentStat
    setTxtProgressBar(pb, i/nnn)
}
close(pb)
# Normal Q-Q plot
#x <- rnorm(n=100000)
#thePlot <- myQQplot(x,temp1)
#ggsave("QQplotOld1.eps",thePlot)
# chi-Squared Q-Q plot
x <- rchisq(n=100000,df=p-r)
x <- (x-(p-r))/sqrt(2*(p-r))
thePlot <- myQQplot(x,myOld)
thePlot2 <- myQQplot(x,myImproved)
ggsave("QQplotOld100.eps",thePlot)
ggsave("QQplotImproved100.eps",thePlot2)
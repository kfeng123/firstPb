source("../stat.R")
library(ggplot2)
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


n <- 300
n1 <- n
n2 <- n
p <- 198
myResult <- vector()
for (i in 1:5000){
    X1 <- matrix(rnorm(n1*p),nrow = n1)
    X2 <- matrix(rnorm(n2*p),nrow = n2)
    myResult[i] <- chenStat(X1,X2,n1,n2)$stat * n1 * n2 /( sqrt(2*p)*(n1+n2) )
}
temp1 <- myResult
# Normal Q-Q plot
x <- rnorm(n=100000)
thePlot <- myQQplot(x,temp1)
ggsave("QQplot1.eps",thePlot)
# chi-Squared Q-Q plot
x <- rchisq(n=100000,df=p)
x <- (x-p)/sqrt(2*p)
thePlot <- myQQplot(x,temp1)
ggsave("QQplotChi1.eps",thePlot)

n <- 50
n1 <- n
n2 <- n
p <- 798
myResult <- vector()
for (i in 1:5000){
    X1 <- matrix(rnorm(n1*p),nrow = n1)
    X2 <- matrix(rnorm(n2*p),nrow = n2)
    myResult[i] <- chenStat(X1,X2,n1,n2)$stat * n1 * n2 /( sqrt(2*p)*(n1+n2) )
}
temp2 <- myResult
# Normal Q-Q plot
x <- rnorm(n=100000)
thePlot <- myQQplot(x,temp2)
ggsave("QQplot2.eps",thePlot)
# chi-Squared Q-Q plot
x <- rchisq(n=100000,df=p)
x <- (x-p)/sqrt(2*p)
thePlot <- myQQplot(x,temp2)
ggsave("QQplotChi2.eps",thePlot)

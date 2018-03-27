set.seed(1)
source("../stat.R")

data1=read.csv("myData.csv")
date=read.csv("myDate.csv",stringsAsFactors = FALSE)
data2=data1[2:nrow(data1),]/data1[1:(nrow(data1)-1),]
data2=rbind(1,data2)
D=data2
D=log(data2)

theMissing=colMeans(is.na(D))
D=D[,theMissing==0]



theWeekDay=format(as.Date(date[,1]),"%u")
X1=D[theWeekDay==1,]
X2=D[theWeekDay!=1,]
X1=as.matrix(X1)
X2=as.matrix(X2)
p=ncol(X1)
n1=nrow(X1)
n2=nrow(X2)

# precalculate the eigenvalue decomposition
S <- ((n1 - 1) * var(X1) + (n2 - 1) * var(X2)) / (n1 + n2 - 2)
myEigen <- eigen(S, symmetric = TRUE)
        
# estimate r
estR <- which.max(myEigen$values[1:50]/myEigen$values[2:51])

tmp <- myStatFinal2(X1, X2, n1, n2, r=estR, myEigen=myEigen)
JJJ <- tmp$studentStat

myPvalue <- pchisq(JJJ*sqrt(2*(p-estR))+p-estR,df=p-estR,lower.tail=FALSE)
chiPvalue <- myChiPvalue(X1,X2,n1,n2,r=estR)
fastPvalue <- maPvalue(X1,X2,n1,n2,r=estR)
cqPvalue <- pnorm(ChenQin.test(X1,X2)$ChQ,lower.tail = FALSE)
sdPvalue <- pnorm(sdStat(X1,X2,n1,n2)$stat,lower.tail=FALSE)
ljwPvalue <- ljwPvalue(X1,X2,n1,n2)


temp <- list(myPvalue,chiPvalue,fastPvalue,cqPvalue,sdPvalue,ljwPvalue)

setwd("C:\\Users\\cspar\\Dropbox\\research files\\2024 group equality of indirect effects\\hayes2018data\\disaster")
D=read.csv("disaster.csv")[,-1]
names(D)=c("X","Y","M","J")

# define categorical moderator G and GX
QJ=quantile(D$J,prob=c(0.25,0.5,0.75),type=6);QJ
D$W <- ifelse(D$J <= QJ[1], 0,
       ifelse(D$J <= QJ[2], 1,
       ifelse(D$J <= QJ[3], 2,3)))
D$G1=ifelse(D$W==1,1,0)
D$G2=ifelse(D$W==2,1,0)
D$G3=ifelse(D$W==3,1,0)
D$G1X=D$G1*D$X
D$G2X=D$G2*D$X
D$G3X=D$G3*D$X
head(D,20)

M.M=lm(M~X*factor(W),D)    # equation (27)
R.M=summary(M.M)
R.M$coefficients

# OLS estimates of model M
a0=R.M$coefficients[1,1]
a1=R.M$coefficients[2,1]
a21=R.M$coefficients[3,1]
a22=R.M$coefficients[4,1]
a23=R.M$coefficients[5,1]
a31=R.M$coefficients[6,1]
a32=R.M$coefficients[7,1]
a33=R.M$coefficients[8,1]
a0;a1;a21;a22;a23;a31;a32;a33

M.Y=lm(Y~X*factor(W)+ M,D)  # equation (28)
R.Y=summary(M.Y)
R.Y$coefficients

# OLS estimates of model Y
cp0=R.Y$coefficients[1,1]
cp1=R.Y$coefficients[2,1]
cp21=R.Y$coefficients[3,1]
cp22=R.Y$coefficients[4,1]
cp23=R.Y$coefficients[5,1]
b1=R.Y$coefficients[6,1]
cp31=R.Y$coefficients[7,1]
cp32=R.Y$coefficients[8,1]
cp33=R.Y$coefficients[9,1]
cp0;cp1;cp21;cp22;cp23;cp31;cp32;cp33;b1

# indirect effects for groups by the product method
IND.0=b1*a1;IND.1=b1*a1+b1*a31;IND.2=b1*a1+b1*a32;IND.3=b1*a1+b1*a33
IND.0;IND.1;IND.2;IND.3

D$Ystar=D$Y-(cp0+cp1*D$X+cp21*D$G1+cp22*D$G2
             +cp23*D$G3+cp31*D$G1X+cp32*D$G2X+cp33*D$G3X)  # equation (30)
M.Yadj=lm(Ystar~X*factor(W),D)
R.Yadj=summary(M.Yadj)

# OLS estimates of model Y*
h0=R.Yadj$coefficients[1,1]
h1=R.Yadj$coefficients[2,1]
h21=R.Yadj$coefficients[3,1]
h22=R.Yadj$coefficients[4,1]
h23=R.Yadj$coefficients[5,1]
h31=R.Yadj$coefficients[6,1]
h32=R.Yadj$coefficients[7,1]
h33=R.Yadj$coefficients[8,1]
h0;h1;h21;h22;h23;h31;h32;h33

# indirect effects by the adjusted response method
INDadj.0=h1;INDadj.1=h1+h31;INDadj.2=h1+h32;INDadj.3=h1+h33
INDadj.0;INDadj.1;INDadj.2;INDadj.3
# indirect effect by the product method
IND.0;IND.1;IND.2;IND.3
# Check the equivalence relations by comparing 
# INDadj.0;INDadj.1;INDadj.2;INDadj.3 with IND.0;IND.1;IND.2;IND.3

H.origin=R.Yadj$coefficients[6:8,1];H.origin # relative conditional indirect effects

M.ANCOVA=lm(Ystar~X*factor(W),D)  # equation (32)
summary(M.ANCOVA)
R.ANCOVA=anova(M.ANCOVA)
F.origin=R.ANCOVA$"Mean Sq"[3]/R.ANCOVA$"Mean Sq"[4]  # equation (34)
F.origin

D$Ystarstar=D$Ystar-(h31*D$G1X+h32*D$G2X+h33*D$G3X) # double-adjusted response in equation (35)

####### bootstrapping ##############
set.seed(123479)
nsize=dim(D)[1]
n.boot=10000              # 10000 number of bootstrap replication
H.boot=matrix(NA,n.boot,3)  # bootstrap estimates of relative indirect effects by groups
Null.F=c()                # F_boot statistic under nulll hypothesis
library(car)

for (i in 1:n.boot){
index=sample(1:nsize,nsize,replace=T)
D.boot=D[index,]

M.M=lm(M~X*factor(W),D.boot)
R.M=summary(M.M)
a0=R.M$coefficients[1,1]
a1=R.M$coefficients[2,1]
a21=R.M$coefficients[3,1]
a22=R.M$coefficients[4,1]
a23=R.M$coefficients[5,1]
a31=R.M$coefficients[6,1]
a32=R.M$coefficients[7,1]
a33=R.M$coefficients[8,1]
a0;a1;a21;a22;a23;a31;a32;a33
 
M.Y=lm(Y~X*factor(W)+ M,D.boot)
R.Y=summary(M.Y)
cp0=R.Y$coefficients[1,1]
cp1=R.Y$coefficients[2,1]
cp21=R.Y$coefficients[3,1]
cp22=R.Y$coefficients[4,1]
cp23=R.Y$coefficients[5,1]
b1=R.Y$coefficients[6,1]
cp31=R.Y$coefficients[7,1]
cp32=R.Y$coefficients[8,1]
cp33=R.Y$coefficients[9,1]
cp0;cp1;cp21;cp22;cp23;cp31;cp32;cp33;b1

# calculation of Y^*
D.boot$Ystar=D.boot$Y-(cp0+cp1*D.boot$X+cp21*D.boot$G1+cp22*D.boot$G2++cp23*D.boot$G3+
                       cp31*D.boot$G1X+cp32*D.boot$G2X+cp33*D.boot$G3X)

M.Yadj=lm(Ystar~X*factor(W),D.boot)
R.Yadj=summary(M.Yadj)
h0=R.Yadj$coefficients[1,1]
h1=R.Yadj$coefficients[2,1]
h21=R.Yadj$coefficients[3,1]
h22=R.Yadj$coefficients[4,1]
h23=R.Yadj$coefficients[5,1]
h31=R.Yadj$coefficients[6,1]
h32=R.Yadj$coefficients[7,1]
h33=R.Yadj$coefficients[8,1]

H.boot[i,1]=h31
H.boot[i,2]=h32
H.boot[i,3]=h33

# generation of pseudo-population
D.boot$Ystarstar.boot=D$Ystarstar[sample(1:nsize,nsize,replace=T)]
Null=lm(Ystarstar.boot~X*factor(W),D.boot)
anova(Null)
Null.F[i]=anova(Null)$"F value"[3]
}

windows(record=T)
hist(H.boot[,1],breaks=200,freq=F,main="",xlab="")
title(expression(paste("(a) Distribution of ",hat(h)[31])))
H1.low=quantile(H.boot[,1],prob=0.025,type=6)
H1.high=quantile(H.boot[,1],prob=0.975,type=6)
segments(H.origin[1],0,H.origin[1],1.7,lty=2,lwd=2)
segments(H1.low,0,H1.low,0.6,lty=2,lwd=2)
segments(H1.high,0,H1.high,0.6,lty=2,lwd=2)
mtext(paste(round(H.origin[1],4)),1,-1,at=H.origin[1],cex=0.8)
mtext(paste(round(H1.low,4)),1,-1,at=H1.low,cex=0.8)
mtext(paste(round(H1.high,4)),1,-1,at=H1.high,cex=0.8)
H1.low;H1.high

hist(H.boot[,2],breaks=200,freq=F,main="",xlab="")
title(expression(paste("(b) Distribution of ",hat(h)[32])))
H2.low=quantile(H.boot[,2],prob=0.025,type=6)
H2.high=quantile(H.boot[,2],prob=0.975,type=6)
segments(H.origin[2],0,H.origin[2],1.6,lty=2,lwd=2)
segments(H2.low,0,H2.low,0.6,lty=2,lwd=2)
segments(H2.high,0,H2.high,0.6,lty=2,lwd=2)
mtext(paste(round(H.origin[2],4)),1,-1,at=H.origin[2],cex=0.8)
mtext(paste(round(H2.low,4)),1,-1,at=H2.low,cex=0.8)
mtext(paste(round(H2.high,4)),1,-1,at=H2.high,cex=0.8)
H2.low;H2.high

hist(H.boot[,3],breaks=200,freq=F,main="",xlab="")
title(expression(paste("(c) Distribution of ",hat(h)[33])))
H3.low=quantile(H.boot[,3],prob=0.025,type=6)
H3.high=quantile(H.boot[,3],prob=0.975,type=6)
segments(H.origin[3],0,H.origin[3],1.6,lty=2,lwd=2)
segments(H3.low,0,H3.low,0.6,lty=2,lwd=2)
segments(H3.high,0,H3.high,0.6,lty=2,lwd=2)
mtext(paste(round(H.origin[3],4)),1,-1,at=H.origin[3],cex=0.8)
mtext(paste(round(H3.low,4)),1,-1,at=H3.low,cex=0.8)
mtext(paste(round(H3.high,4)),1,-1,at=H3.high,cex=0.8)
H3.low;H3.high

hist(Null.F,breaks=200,freq=F,main="",xlab="")
title(expression(paste("(d) Distribution of bootstrap F-statistic")))
segments(F.origin,0,F.origin,0.9,lty=2,lwd=2)
mtext(paste(round(F.origin,4)),1,-1,at=F.origin,cex=0.8)
p.value=1-ecdf(Null.F)(F.origin);p.value
text(F.origin+1, 0.7,paste0("p=",round(p.value,4)))
############  end #########################

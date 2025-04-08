rm(list=ls())
D<-read.csv("https://raw.githubusercontent.com/Changsoon-Park/Omnibus-test-for-indirect-effects/main/protest.csv")[,c(2,5,6)]
head(D)
names(D)=c("X","Y","M")

D$D1=ifelse(D$X==1,1,0)
D$D2=ifelse(D$X==2,1,0)
head(D,20)

M1=lm(M~factor(X),D)   # equation (10)
R1=summary(M1)

# OLS estimates of model M
a0=R1$coefficients[1,1]
a11=R1$coefficients[2,1]
a12=R1$coefficients[3,1]


M2=lm(Y~factor(X)+M,D)   # equation (11)
R2=summary(M2)

# OLS estimates of model Y
cp0=R2$coefficients[1,1]
cp11=R2$coefficients[2,1]
cp12=R2$coefficients[3,1]
b1=R2$coefficients[4,1]
a0;a11;a12;cp0;cp11;cp12;b1

D$Ystar=D$Y-(cp0+cp11*D$D1+cp12*D$D2) # equation (13)
M.Ystar=lm(Ystar~factor(X),D)
R.Ystar=summary(M.Ystar)

# OLS estimates of model Y*
h0=R.Ystar$coefficients[1.1]
h11=R.Ystar$coefficients[2.1]
h12=R.Ystar$coefficients[3.1]

# indirect effects by the adjusted response method
# for reference group, relative indirect effects in equation (14) 
h0;h11;h12  
RIE.h11=h11;RIE.h12=h12
# Indirect effects by the product method
b1*a0;b1*a11;b1*a12
# check equivalence relations by comparing 
# h0;h11;h12 with b1*a0;b1*a11;b1*a12
   
Full=aov(Ystar~factor(X),D)  # equation (15)
summary(Full)
F.origin=anova(Full)$"F value"[1];F.origin  # F_origin

YSS=D$Ystar-(h11*D$D1+h12*D$D2)  # double-adjusted response in equation (18)

######## bootstrapping ##########
nsize=dim(D)[1]
nboot=10000      # 10000 number of bootstrap replication
set.seed(12357)
Null.F=c()       # F_boot statistic under nulll hypothesis
h11.boot=c()  # bootstrap estimate for h11
h12.boot=c()  # bootstrap estimate for h12
h11.low=c()   # bootstrap lower CI estimate for h11
h11.high=c()  # bootstrap upper CI estimate for h11
h12.low=c()   # bootstrap lower CI estimate for h12
h12.high=c()  # bootstrap upper CI estimate for h12

for (i in 1:nboot){
index=sample(1:nsize,nsize,replace=T)
D.boot=D[index,]

M1=lm(M~factor(X),D.boot)
R1=summary(M1)

a0=R1$coefficients[1,1]
a11=R1$coefficients[2,1]
a12=R1$coefficients[3,1]

M2=lm(Y~factor(X)+M,D.boot)
R2=summary(M2)
cp0=R2$coefficients[1,1]
cp11=R2$coefficients[2,1]
cp12=R2$coefficients[3,1]
b1=R2$coefficients[4,1]
a0;a11;a12;cp0;cp11;cp12;b1

D.boot$Ystar=D.boot$Y-(cp0+cp11*D.boot$D1+cp12*D.boot$D2)
M.Ystar=lm(Ystar~factor(X),D.boot)
R.Ystar=summary(M.Ystar)
h0=R.Ystar$coefficients[1.1]
h11=R.Ystar$coefficients[2.1]
h12=R.Ystar$coefficients[3.1]
h0;h11;h12 

h11.boot[i]=h11
h12.boot[i]=h12

# generation of pseudo-population
D.boot$YSS=YSS[sample(1:nsize,nsize,replace=T)]  
Null=aov(YSS~factor(X),D.boot)
Null.F[i]=anova(Null)$"F value"[1];#Null.F
}

windows(record=T)
hist(h11.boot,breaks=200,freq=F,main="",xlab="")
title(expression(paste("(a) Distribution of ",hat(h)[11])))
h11.low=quantile(h11.boot,0.025)
h11.high=quantile(h11.boot,0.975)
segments(RIE.h11,0,RIE.h11,3.5,lty=2,lwd=2)
segments(h11.low,0,h11.low,1.3,lty=2,lwd=2)
segments(h11.high,0,h11.high,1.3,lty=2,lwd=2)
mtext(paste(round(RIE.h11,4)),1,-1,at=RIE.h11,cex=0.8)
mtext(paste(round(h11.low,4)),1,-1,at=h11.low,cex=0.8)
mtext(paste(round(h11.high,4)),1,-1,at=h11.high,cex=0.8)
h11.low;h11.high

hist(h12.boot,breaks=200,freq=F,main="",xlab="")
title(expression(paste("(b) Distribution of ",hat(h)[12])))
h12.low=quantile(h12.boot,0.025)
h12.high=quantile(h12.boot,0.975)
segments(h12,0,h12,3,lty=2,lwd=2)
segments(h12.low,0,h12.low,1.3,lty=2,lwd=2)
segments(h12.high,0,h12.high,1.3,lty=2,lwd=2)
mtext(paste(round(RIE.h12,4)),1,-1,at=RIE.h12,cex=0.8)
mtext(paste(round(h12.low,4)),1,-1,at=h12.low,cex=0.8)
mtext(paste(round(h12.high,4)),1,-1,at=h12.high,cex=0.8)
h12.low;h12.high

hist(Null.F,breaks=200,freq=F,main="",xlab="")
title(expression(paste("(c) Distribution of bootstrap F-statistic")))
segments(F.origin,0,F.origin,0.3,lty=2,lwd=2)
mtext(paste(round(F.origin,4)),1,-1,at=F.origin,cex=0.8)
p.value=1-ecdf(Null.F)(F.origin);p.value
text(F.origin+1, 0.1,paste0("ASL=",round(p.value,4)))
################### end ################################









rm(list=ls())
D<-read.csv("https://raw.githubusercontent.com/Changsoon-Park/Omnibus-test-for-indirect-effects/main/teams.csv")
names(D)=c("X","M","EXP","Y")
head(D)

# define categorical moderator G and GM
QEXP=quantile(D$EXP,prob=c(0.25,0.5,0.75),type=6);QEXP
D$W <- ifelse(D$EXP <= QEXP[1], 0,
       ifelse(D$EXP <= QEXP[2], 1,
       ifelse(D$EXP <= QEXP[3], 2, 3)))
D$G1=ifelse(D$W==1,1,0)
D$G2=ifelse(D$W==2,1,0)
D$G3=ifelse(D$W==3,1,0)
D$G1M=D$G1*D$M
D$G2M=D$G2*D$M
D$G3M=D$G3*D$M
head(D,20)

M.M=lm(M~X,D)         # equation (37)
R.M=summary(M.M)
R.M$coefficients

# OLS estimates of model M
a0=R.M$coefficients[1,1]
a1=R.M$coefficients[2,1]
a0;a1

M.RY=lm(Y~X+G1+G2+G3+M+G1M+G2M+G3M,D) # equation (38)
R.RY=summary(M.RY)
R.RY$coefficients

# OLS estimates of model Y
cp0=R.RY$coefficients[1,1]
cp1=R.RY$coefficients[2,1]
cp21=R.RY$coefficients[3,1]
cp22=R.RY$coefficients[4,1]
cp23=R.RY$coefficients[5,1]
b1=R.RY$coefficients[6.1]
b31=R.RY$coefficients[7.1]
b32=R.RY$coefficients[8.1]
b33=R.RY$coefficients[9.1]
cp0;cp1;cp21;cp22;cp23;b1;b31;b32;b33

IDE.ref=c()
cond.RIDE=c()
h0.T=c()
h1.T=c()
nsize=dim(D)[1]
Y.T=matrix(NA,nsize,3)

for (i in 1:3){
if(i==1){g1=1;g2=0;g3=0;one.1=c(rep(1,nsize),rep(0,2*nsize));X.T1=c(D$X,rep(0,2*nsize))}
if(i==2){g1=0;g2=1;g3=0;one.2=c(rep(0,nsize),rep(1,nsize),rep(0,nsize));X.T2=c(rep(0,nsize),D$X,rep(0,nsize))}
if(i==3){g1=0;g2=0;g3=1;one.3=c(rep(0,2*nsize),rep(1,nsize));X.T3=c(rep(0,2*nsize),D$X)}

D$GC1=D$G1-g1
D$GC2=D$G2-g2
D$GC3=D$G3-g3

D$GC1M=D$GC1*D$M
D$GC2M=D$GC2*D$M
D$GC3M=D$GC3*D$M

M.Y=lm(Y~X+GC1+GC2+GC3+M+GC1M+GC2M+GC3M,D)   # y_{j}in equation (40)
R.Y=summary(M.Y)
R.Y$coefficients

# OLS estimates of model y_{j}
cp0C=R.Y$coefficients[1,1]
cp1C=R.Y$coefficients[2,1]
cp21=R.Y$coefficients[3,1]
cp22=R.Y$coefficients[4,1]
cp23=R.Y$coefficients[5,1]
b1C=R.Y$coefficients[6,1]
b31=R.Y$coefficients[7,1]
b32=R.Y$coefficients[8,1]
b33=R.Y$coefficients[9,1]
cp0C;cp1C;cp21;cp22;cp23;b1C;b31;b32;b33

IDE.ref[i]=b1C*a0  # indirect effect for reference group by product method
cond.RIDE[i]=b1C*a1;cond.RIDE[i]  # conditional relative indirect effect by group i by product method

D$Ystar=D$Y-(cp0C+cp1C*D$X+cp21*D$GC1+cp22*D$GC2+cp23*D$GC3+b31*D$GC1M+b32*D$GC2M+b33*D$GC3M)  # equation (41) 
X.Yadj=lm(Ystar~X,D)
summary(X.Yadj)

# OLS estimates of model (42)
h0=summary(X.Yadj)$coefficients[1,1]
h1=summary(X.Yadj)$coefficients[2,1]
h0;h1  
h0.T[i]=h0   # h_0{j}  in equation (42)
h1.T[i]=h1   # h_1{j}  in equation (42)
Y.T[,i]=D$Ystar   # Ystar_{j} in equation (43)
}

# check equivalence relations for each group
IDE.ref # indirect effect for reference group
h0.T
cond.RIDE  # conditional relative indirect effect for groups by product method
h1.T   # conditional relative indirect effects by adjusted response method
####

Y.T=as.vector(Y.T)   # ystar_T in equation (46)
TOT=data.frame(Y.T,one.1,X.T1,one.2,X.T2,one.3,X.T3)  # use of one vector
REG.TOT=lm(Y.T~-1+.,TOT)          # equation (46)
summary(REG.TOT)$coefficients
anova(REG.TOT)
F.origin=mean(anova(REG.TOT)$"Sum Sq"[c(2,4,6)])/anova(REG.TOT)$"Mean Sq"[7];F.origin  # equation (47)
TOT$Y.T.S=resid(REG.TOT)  # double-adjusted response in equation (47)

####### derivation of F-statistic from three separate models ################
TOT1=TOT[1:nsize,];TOT2=TOT[(nsize+1):(2*nsize),];TOT3=TOT[(2*nsize+1):(3*nsize),]
TOT1
TOT2
TOT3
REG.1=lm(Y.T~X.T1,TOT1)
summary(REG.1)
anova(REG.1)
REG.2=lm(Y.T~X.T2,TOT2)
anova(REG.2)
REG.3=lm(Y.T~X.T3,TOT3)
anova(REG.3)
MSReg.TOT=(anova(REG.1)$"Sum Sq"[1]+anova(REG.2)$"Sum Sq"[1]+anova(REG.3)$"Sum Sq"[1])/3
MSRes.TOT=(anova(REG.1)$"Sum Sq"[2]+anova(REG.2)$"Sum Sq"[2]+anova(REG.3)$"Sum Sq"[2])/174
F.origin=MSReg.TOT/MSRes.TOT;F.origin
#####################################################

####  bootstrapping ##########
set.seed(13479)
nsize=dim(D)[1]
n.boot=10000    # 10000 bootstrao replication

cond.IND.boot=matrix(NA,n.boot,3)
h0.boot=matrix(NA,n.boot,3)  # bootstrap estimate of h_0
h1.boot=matrix(NA,n.boot,3)  # bootstrap estimate of h_1
F.boot=c()                   # F_boot statistic under nulll hypothesis

for (i in 1:n.boot){
index=sample(1:nsize,nsize,replace=T)
D.boot=D[index,]

for (k in 1:3){
if(k==1){g1=1;g2=0;g3=0;one.1=c(rep(1,nsize),rep(0,2*nsize));X.T1=c(D.boot$X,rep(0,2*nsize))}
if(k==2){g1=0;g2=1;g3=0;one.2=c(rep(0,nsize),rep(1,nsize),rep(0,nsize));X.T2=c(rep(0,nsize),D.boot$X,rep(0,nsize))}
if(k==3){g1=0;g2=0;g3=1;one.3=c(rep(0,2*nsize),rep(1,nsize));X.T3=c(rep(0,2*nsize),D.boot$X)}

D.boot$GC1=D.boot$G1-g1
D.boot$GC2=D.boot$G2-g2
D.boot$GC3=D.boot$G3-g3
D.boot$GC1M=D.boot$GC1*D.boot$M
D.boot$GC2M=D.boot$GC2*D.boot$M
D.boot$GC3M=D.boot$GC3*D.boot$M

M.M=lm(M~X,D.boot)
R.M=summary(M.M)
R.M$coefficients
a0=R.M$coefficients[1,1]
a1=R.M$coefficients[2,1]
a0;a1

M.Y=lm(Y~X+GC1+GC2+GC3+M+GC1M+GC2M+GC3M,D.boot)
R.Y=summary(M.Y)
R.Y$coefficients
if(length(R.Y$coefficients[,1])!=9) next

cp0C=R.Y$coefficients[1,1]
cp1=R.Y$coefficients[2,1]
cp21=R.Y$coefficients[3,1]
cp22=R.Y$coefficients[4,1]
cp23=R.Y$coefficients[5,1]
b1C=R.Y$coefficients[6,1]
b31=R.Y$coefficients[7,1]
b32=R.Y$coefficients[8,1]
b33=R.Y$coefficients[9,1]
cp0C;cp1;cp21;cp22;cp23;b1C;b31;b32;b33

cond.IND.boot[i,k]=b1C*a1;cond.IND.boot[i,k]

# calculation of Y^*
D.boot$Ystar=D.boot$Y-(cp0C+cp1*D.boot$X+cp21*D.boot$GC1+cp22*D.boot$GC2+
     cp23*D.boot$GC3+b31*D.boot$GC1M+b32*D.boot$GC2M+b33*D.boot$GC3M)
X.Yadj=lm(Ystar~X,D.boot)
summary(X.Yadj)
h0=summary(X.Yadj)$coefficients[1,1]
h1=summary(X.Yadj)$coefficients[2,1]
h0;h1  
h0.boot[i,k]=h0
h1.boot[i,k]=h1
}

T.index=sample(1:3*nsize,3*nsize,replace=T)
Y.T.S.boot=TOT$Y.T.S[T.index]               # i-th YSS bootstrap sample under null hypothesis (48)
TOT.boot=data.frame(Y.T.S.boot,one.1,X.T1,one.2,X.T2,one.3,X.T3)  # use of one vector
REG.TOT=lm(Y.T.S.boot~-1+.,TOT.boot)
summary(REG.TOT)$coefficients
anova(REG.TOT)
F.boot[i]=mean(anova(REG.TOT)$"Sum Sq"[c(2,4,6)])/anova(REG.TOT)$"Mean Sq"[7]  # i-th bootstrap F-statistic
}

cond.IND.boot
h1.boot=na.omit(h1.boot)
F.boot
p.boot=1-ecdf(F.boot)(F.origin);p.boot


windows(record=T)
hist(h1.boot[,1],breaks=200,freq=F,main="",xlab="",xlim=c(-1.5,1))
title(expression(paste("(a) Distribution of ",hat(h)["1{1}"])))
h11.low=quantile(h1.boot[,1],prob=0.025,type=6)
h11.high=quantile(h1.boot[,1],prob=0.975,type=6)
segments(h1.T[1],0,h1.T[1],4,lty=2,lwd=2)
segments(h11.low,0,h11.low,1,lty=2,lwd=2)
segments(h11.high,0,h11.high,1,lty=2,lwd=2)
mtext(paste(round(h1.T[1],4)),1,-1,at=h1.T[1],cex=0.8)
mtext(paste(round(h11.low,4)),1,-1,at=h11.low,cex=0.8)
mtext(paste(round(h11.high,4)),1,-1,at=h11.high,cex=0.8)
h11.low;h11.high

hist(h1.boot[,2],breaks=200,freq=F,main="",xlab="",xlim=c(-1.5,1))
title(expression(paste("(b) Distribution of ",hat(h)["1{2}"])))
h12.low=quantile(h1.boot[,2],prob=0.025,type=6)
h12.high=quantile(h1.boot[,2],prob=0.975,type=6)
segments(h1.T[2],0,h1.T[2],4,lty=2,lwd=2)
segments(h12.low,0,h12.low,1,lty=2,lwd=2)
segments(h12.high,0,h12.high,1,lty=2,lwd=2)
mtext(paste(round(h1.T[2],4)),1,-1,at=h1.T[2],cex=0.8)
mtext(paste(round(h12.low,4)),1,-1,at=h12.low,cex=0.8)
mtext(paste(round(h12.high,4)),1,-1,at=h12.high,cex=0.8)
h12.low;h12.high

hist(h1.boot[,3],breaks=200,freq=F,main="",xlab="",xlim=c(-1.5,0.5))
title(expression(paste("(c) Bootstrap Distribution of ",hat(h)["1{3}"])))
h13.low=quantile(h1.boot[,3],prob=0.025,type=6)
h13.high=quantile(h1.boot[,3],prob=0.975,type=6)
segments(h1.T[3],0,h1.T[3],4,lty=2,lwd=2)
segments(h13.low,0,h13.low,1,lty=2,lwd=2)
segments(h13.high,0,h13.high,1,lty=2,lwd=2)
mtext(paste(round(h1.T[3],4)),1,-1,at=h1.T[3],cex=0.8)
mtext(paste(round(h13.low,4)),1,-1,at=h13.low,cex=0.8)
mtext(paste(round(h13.high,4)),1,-1,at=h13.high,cex=0.8)
h13.low;h13.high

hist(F.boot,breaks=200,freq=F,main="",xlab="")
title("(d) Distribution of bootstrap F-Statistic")
segments(F.origin,0,F.origin,0.3,lty=2,lwd=2)
mtext(paste(round(F.origin,4)),1,-1,at=F.origin,cex=0.8)
p.value=1-ecdf(F.boot)(F.origin);p.value
text(F.origin+1, 0.2,paste0("ASL=",round(p.value,4)))
############ end ###################################
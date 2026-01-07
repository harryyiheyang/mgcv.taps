require(mgcv)
library(devtools)
document()
devtools::load_all()
#set.seed(100)
#dat <- gamSim(1,n=500,scale=2)
n=30000
PV=matrix(0,1000,5)
for(i in 1:1000){
X=MASS::mvrnorm(n,rep(0,4),matrix(0.25,4,4)+0.75*diag(4))
x1=qbeta(pnorm(X[,1]),1.5,1.5)
x2=qbeta(pnorm(X[,2]),1.5,1.5)
x3=qbeta(pnorm(X[,3]),1.5,1.5)
x4=qbeta(pnorm(X[,4]),1.5,1.5)
f1=x1*2
t2=2*pi*x2
f2=0.4*sin(t2)+0.8*cos(t2)+1.2*sin(t2)^2+1.6*cos(t2)^3+2*sin(t2)^3
t3=2*(x3-0.5)
f3=3*sin(3*t3)+6*exp(-36*t3^2)
f4=0*x4
eta=f1+f2+f3
y=1+eta+rnorm(n,0,1)*sd(eta)
dat=data.frame(x0=x1,x1=x2,x2=x3,x3=x4,y=y)

b<-bam(y~s(x0,bs="AMatern",k=10,m=300,xt=list(getA=polynomial,para=2))+s(x1,bs="cr")+s(x2,bs="cr")+s(x3,bs="cr"),data=dat,family="gaussian",method="REML")
fitb=summary(b)
fit1=taps_score_test(b,test.component=1)
fit2=taps_score_test(b,test.component=1,method="liu")
fit3=taps_score_test(b,test.component=1,method="davies")
fit4=taps_score_test(b,test.component=1,method="farebrother")
fit5=taps_score_test(b,test.component=1,method="imhof")
PV[i,]=c(fit1$smooth.pvalue,fit2$smooth.pvalue,fit3$smooth.pvalue,fit4$smooth.pvalue,fit5$smooth.pvalue)
if(i%%50==0) print(i)
}

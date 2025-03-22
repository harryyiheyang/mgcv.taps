require(mgcv)
library(devtools)
document()
devtools::load_all()
#set.seed(100)
#dat <- gamSim(1,n=500,scale=2)
n=500
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


b<-bam(y~s(x0,bs="Amatern",k=10,m=300,xt=list(getA=polynomial,para=2))+s(x1,bs="Amatern")+s(x2,bs="matern")+s(x3,bs="matern"),data=dat,family="gaussian",method="REML")
fitb=summary(b)
plot(b,pages=1)
summary(b)
mgcv_maps_wald(b)

b1<-bam(y~s(x0,bs="matern")+s(x1,bs="matern")+s(x2,bs="matern")+s(x3,bs="matern"),data=dat,family="gaussian",method="REML")
plot(b1,pages=1)
summary(b1)
mgcv_maps_wald(b1)

b2<-bam(y~s(x0,bs="cosine",xt=list(method="beta",k=250))+s(x1,bs="cosine",xt=list(method="beta",k=250))+s(x2,bs="cosine",xt=list(method="beta",k=250))+s(x3,bs="cosine",xt=list(method="beta",k=250)),data=dat,family="gaussian",method="REML")
plot(b2,pages=1)
summary(b2)
mgcv_maps_wald(b2)

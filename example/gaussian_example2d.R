require(mgcv)
library(scatterplot3d)
library(devtools)
document()
devtools::load_all()
#set.seed(100)
#dat <- gamSim(1,n=500,scale=2)
n=5000
X=MASS::mvrnorm(n,rep(0,4),matrix(0.25,4,4)+0.75*diag(4))
x1=qbeta(pnorm(X[,1]),0.5,0.5)
x2=qbeta(pnorm(X[,2]),0.5,0.5)
x3=qbeta(pnorm(X[,3]),0.5,0.5)
x4=qbeta(pnorm(X[,4]),0.5,0.5)
f1=x1
f2=x2+x1*x2
t3=2*(x3-0.5)
f3=3*sin(3*t3)+6*exp(-36*t3^2)
f4=0*x4
eta=f1+f2+f3
y=1+eta+rnorm(n,0,1)*sd(eta)
dat=data.frame(x0=x1,x1=x2,x2=x3,x3=x4,y=y)

b<-bam(y~s(x0,x1,bs="A2Matern",k=15,m=600)+s(x2,bs="Matern")+s(x3,bs="Matern"),data=dat,family="gaussian",method="REML",gamma=1)
fitb=summary(b)
vis.gam(b, view = c("x0", "x1"))
summary(b)
taps_wald_test(b,test.component = 1)
taps_score_test(b,test.component = 1)

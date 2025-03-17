require(mgcv)
library(devtools)
document()
devtools::load_all()
#set.seed(100)
dat <- gamSim(1,n=1000,scale=1,dist = "normal")
dat$y = dat$y - dat$f0 + 2*dat$x0

b<-bam(y~s(x0,bs="Amatern",k=10,m=300)+s(x1,bs="matern")+s(x2,bs="matern")+s(x3,bs="matern"),data=dat,family="gaussian",select=F)
fitb=summary(b)
plot(b,pages=1)
summary(b)
mgcv_maps_wald(b,approx.method="Gamma")

b1<-bam(y~s(x0,bs="matern")+s(x1,bs="matern")+s(x2,bs="matern")+s(x3,bs="matern"),data=dat,family="gaussian",select=F)
plot(b1,pages=1)
summary(b1)
mgcv_maps_wald(b1)

b2<-bam(y~s(x0,bs="gp")+s(x1,bs="gp")+s(x2,bs="gp")+s(x3,bs="gp"),data=dat,family="gaussian",select=F)
plot(b1,pages=1)
summary(b2)
mgcv_maps_wald(b2)

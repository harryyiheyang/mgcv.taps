require(mgcv)
devtools::load_all()
set.seed(100)
dat <- gamSim(1,n=1000,scale=2)

b<-gam(y~s(x0,bs="cosine",k=10,m=c(2,3))+s(x1,bs="cosine",k=10,m=c(2,3))+
           +          s(x2,bs="cosine",k=10,m=c(2,3))+s(x3,bs="cosine",k=10,m=c(1,9)),data=dat,method="REML")
plot(b,pages=1)

require(mgcv)
devtools::load_all()
#set.seed(100)
dat <- gamSim(1,n=1000,scale=1,dist = "binary")
softthres=function(x,b){
y=x-b
y[y<0]=0
return(y)
}

linear_changepoint=function(x,para){
knot=para
p=length(knot)
G=matrix(0,length(x),p)
for(i in 1:p){
G[,i]=softthres(x,knot[i])
}
return(cbind(1,x,G))
}

b<-bam(y~s(x0,bs="mix_matern",k=10,m=c(100,10),xt=list(getA=linear_changepoint,para=0.5))+s(x1,bs="mix_matern",k=10,m=c(100,10),xt=list(getA=linear_changepoint,para=0.4))
       +s(x2,bs="mix_matern",k=10,m=c(100,10),xt=list(getA=linear_changepoint,para=c(0.24,0.5,0.7)))+s(x3,bs="mix_matern",k=10,m=c(100,10)),data=dat,method="REML",family="binomial",select=F)
plot(b,pages=1)

summary(b)

b1<-bam(y~s(x0,bs="cr")+s(x1,bs="cr")
       +s(x2,bs="cr")+s(x3,bs="cr"),data=dat,method="REML",family="binomial",select=F)
plot(b1,pages=1)
summary(b1)

anova(b1,b)
